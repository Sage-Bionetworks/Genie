from __future__ import absolute_import
from genie import FileTypeFormat, process_functions
import os
import logging
import pandas as pd
import subprocess
logger = logging.getLogger(__name__)

# def createGenePositionsTables():
#   from biomart import BiomartServer
#   import synapseclient
#   syn = synapseclient.login()
#   print(time.time())
#   #server = BiomartServer( "feb2014.archive.ensembl.org/biomart" )
#   server = BiomartServer( "grch37.ensembl.org/biomart" )

#   hsapiens_gene_ensembl = server.datasets['hsapiens_gene_ensembl']
#   #hsapiens_gene_ensembl.show_attributes()

#   chromosomes = range(1,23)
#   chromosomes.extend(["X","Y"])
#   chromosomes = [str(i) for i in chromosomes]
#   response = hsapiens_gene_ensembl.search({
#   'attributes': ['hgnc_symbol',"chromosome_name","start_position","end_position"]
#   })
#   genes = pd.DataFrame()
#   for line in response.iter_lines():
#       line = line.decode('utf-8')
#       row = line.split("\t")
#       if row[0] != '' and row[1] in chromosomes:
#           tempRow = pd.DataFrame([row])
#           genes = genes.append(tempRow)
#   genes.columns = ['hgnc_symbol',"chromosome_name","start_position","end_position"]
#   #genes.to_csv("genepositions.csv",index=False)
#   print(time.time())
#   databaseSynId = "syn11806563"
#   databaseEnt = syn.get(databaseSynId)
#   database = syn.tableQuery("SELECT * FROM %s" % databaseSynId)
#   database = database.asDataFrame()
#   process_functions.updateDatabase(syn, database, genes, databaseSynId, databaseEnt.primaryKey, toDelete=True)
#   return(genes)


def validateSymbol(x, genePositionDf, returnMappedDf=False):
    '''
    The apply function of a DataFrame is called twice on the first row (known
    pandas behavior)
    '''
    valid = True
    matchWithGeneDb = genePositionDf[
        genePositionDf['hgnc_symbol'] == x['Hugo_Symbol']]
    '''
    Make sure there are no overlaps in the submitted gene symbol and the gene
    position database
    '''
    toMap = True
    if not matchWithGeneDb.empty:
        if (matchWithGeneDb['start_position'].values[0] >= x['Start_Position'] and matchWithGeneDb['end_position'].values[0] <= x['End_Position']) or (matchWithGeneDb['start_position'].values[0] <= x['End_Position'] and matchWithGeneDb['end_position'].values[0] >= x['Start_Position']):
            toMap = False
    if toMap:
        chromRows = genePositionDf[genePositionDf['chromosome_name'] == str(x['Chromosome'])]
        startRows = chromRows[chromRows['start_position'] <= x['Start_Position']]

        endRows = startRows[startRows['end_position'] >= x['End_Position']]
        # geneDiff = chromRows['end_position'] - chromRows['start_position']
        bedLength = x['End_Position'] - x['Start_Position']
        # as long as the strand is 90% within the boundary of the
        # Start goes over start boundary, but end is contained in gene
        if len(endRows) == 0:
            if sum(chromRows['end_position'] >= x['End_Position']) > 0:
                overlap = x['End_Position'] - chromRows['start_position']
                # difference =  difference * -1.0
                ratioOverlap = overlap / bedLength
                ratioOverlap = ratioOverlap[ratioOverlap > 0.9]
                ratioOverlap = ratioOverlap[ratioOverlap <= 1]
                if not ratioOverlap.empty:
                    endRows = endRows.append(chromRows.loc[[ratioOverlap.idxmax()]])
            # End goes over end boundary, but start is contained in gene
            if sum(chromRows['start_position'] <= x['Start_Position']) > 0:
                overlap = chromRows['end_position'] - x['Start_Position']
                ratioOverlap = overlap / bedLength
                ratioOverlap = ratioOverlap[ratioOverlap > 0.9]
                ratioOverlap = ratioOverlap[ratioOverlap <= 1]
                if not ratioOverlap.empty:
                    endRows = endRows.append(chromRows.loc[[ratioOverlap.idxmax()]])
            # Start and end go over gene boundary
            check = chromRows[chromRows['start_position'] >= x['Start_Position']]
            check = check[check['end_position'] <= x['End_Position']]
            if not check.empty:
                overlap = chromRows['end_position'] - chromRows['start_position']
                ratioOverlap = overlap / bedLength
                ratioOverlap = ratioOverlap[ratioOverlap > 0.9]
                ratioOverlap = ratioOverlap[ratioOverlap <= 1]
                if not ratioOverlap.empty:
                    endRows = endRows.append(chromRows.loc[[ratioOverlap.idxmax()]])

        if len(endRows) == 0:
            logger.warning("%s cannot be remapped. These rows will have an empty gene symbol" % x['Hugo_Symbol'])
            x['Hugo_Symbol'] = pd.np.nan
            valid = False
        elif len(endRows) > 1:
            if x['Hugo_Symbol'] not in endRows['hgnc_symbol'].tolist():
                # if "MLL4", then the HUGO symbol should be KMT2D and KMT2B
                logger.warning("%s can be mapped to different symbols: %s. Please correct or it will be removed." % (x['Hugo_Symbol'], ", ".join(endRows['hgnc_symbol'])))
                x['Hugo_Symbol'] = pd.np.nan
                valid = False
        else:
            if x['Hugo_Symbol'] != endRows['hgnc_symbol'].values[0]:
                logger.info("%s will be remapped to %s" % (x['Hugo_Symbol'], endRows['hgnc_symbol'].values[0]))
                x['Hugo_Symbol'] = endRows['hgnc_symbol'].values[0]
    if returnMappedDf:
        return(x)
    else:
        return(valid)


class bed(FileTypeFormat):

    _fileType = "bed"

    _process_kwargs = ["newPath", "parentId", "databaseSynId",'seq_assay_id']

    def _get_dataframe(self, filePathList):
        filePath = filePathList[0]
        try:
            beddf = pd.read_csv(filePath, sep="\t", header=None)
        except Exception:
            raise ValueError("Can't read in your bed file. Please make sure the BED file is not binary and does not contain a comment/header line")
        if not str(beddf[0][0]).isdigit() and not str(beddf[0][0]).startswith("chr"):
            raise ValueError("Please make sure your bed file does not contain a comment/header line")
        return(beddf)

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]).startswith("%s-" % self.center) and os.path.basename(filePath[0]).endswith(".bed")

    def createdBEDandGenePanel(self, bed, seq_assay_id,
                               genePanelPath, parentId,
                               createGenePanel=True):
        logger.info("REMAPPING %s" % seq_assay_id)
        # bedname = seq_assay_id + ".bed"
        bed.columns = ["Chromosome", "Start_Position", "End_Position", 
                       "Hugo_Symbol", "includeInPanel"]
        # Validate gene symbols
        # Gene symbols can be split by ; and _ and : and .
        bed['Hugo_Symbol'] = [i.split(";")[0].split("_")[0].split(":")[0].split(".")[0] for i in bed['Hugo_Symbol']]
        bed['Chromosome'] = [str(i).replace("chr","") for i in bed['Chromosome']]
        bed['Start_Position'] = bed['Start_Position'].apply(int)
        bed['End_Position'] = bed['End_Position'].apply(int)

        genePosition = self.syn.tableQuery('SELECT * FROM syn11806563')
        genePositionDf = genePosition.asDataFrame()
        bed['ID'] = bed['Hugo_Symbol']
        bed = bed.apply(lambda x: validateSymbol(x, genePositionDf, returnMappedDf=True), axis=1)

        bed['SEQ_ASSAY_ID'] = seq_assay_id
        bed.to_csv(os.path.join(process_functions.SCRIPT_DIR,"temp.bed"), sep="\t",index=False, header=None)
        command = ['bedtools','intersect','-a',os.path.join(process_functions.SCRIPT_DIR,"temp.bed"),'-b',os.path.join(process_functions.SCRIPT_DIR,'exon.gtf'),'-u','>',os.path.join(process_functions.SCRIPT_DIR,'genie_exons.bed')]
        subprocess.check_call(" ".join(command),shell=True)
        if os.stat(os.path.join(process_functions.SCRIPT_DIR,'genie_exons.bed')).st_size > 0 and createGenePanel:
            temp = pd.read_csv(os.path.join(process_functions.SCRIPT_DIR,'genie_exons.bed'),sep="\t",header=None)
            temp.columns = ["Chromosome","Start_Position","End_Position","Hugo_Symbol","includeInPanel","ID","SEQ_ASSAY_ID"]
            #Only include genes that should be included in the panels
            temp = temp[temp['includeInPanel'] == True]
            # Write gene panel
            allgenes = set(temp['Hugo_Symbol'][~temp['Hugo_Symbol'].isnull()])
            genepanelname = "data_gene_panel_" + seq_assay_id + ".txt"
            metaStableId = "stable_id: %s\n"
            description = "description: %s, Number of Genes - %d\n" % (seq_assay_id, len(allgenes))
            genes = "gene_list:\t%s"
            with open(os.path.join(genePanelPath,genepanelname),"w+") as f:
                f.write(metaStableId % seq_assay_id)
                f.write(description)
                f.write(genes % "\t".join(allgenes))
            process_functions.storeFile(self.syn, os.path.join(genePanelPath,genepanelname), parentId=parentId, center = self.center, fileFormat="bed", dataSubType="metadata",cBioFileFormat="genePanel")
        return(bed)

    def _process(self, gene, seq_assay_id, newPath, parentId, createPanel=True):
        seq_assay_id = seq_assay_id.upper()
        seq_assay_id = seq_assay_id.replace('_','-')

        if len(gene.columns) > 4:
            if not all(gene[4].apply(lambda x: isinstance(x, bool))):
                gene[4] = True
        else:
            gene[4] = True
        bed = gene[[0,1,2,3,4]]
        genePanelPath = os.path.dirname(newPath)
        if not os.path.exists(os.path.join(process_functions.SCRIPT_DIR,"exon.gtf")) or not os.path.exists(os.path.join(process_functions.SCRIPT_DIR,"gene.gtf")):
            command = ['bash',os.path.join(process_functions.SCRIPT_DIR,'createGTF.sh')]
            subprocess.check_call(command)
        bed = self.createdBEDandGenePanel(bed, seq_assay_id, genePanelPath, parentId, createGenePanel =createPanel)
        command = ['bash',os.path.join(process_functions.SCRIPT_DIR,'addFeatureType.sh')]
        subprocess.check_call(command)
        bed = pd.read_csv(os.path.join(process_functions.SCRIPT_DIR,"genie_combined.bed"),sep="\t")
        bed['CENTER'] =self.center
        bed['Chromosome'] = bed['Chromosome'].astype(str)
        return(bed)

    def preprocess(self, filePath):
        seq_assay_id = os.path.basename(filePath).replace(".bed", "")
        seq_assay_id = seq_assay_id.upper().replace("_", "-")
        return({'seq_assay_id': seq_assay_id})

    def process_steps(self, gene, newPath, parentId, databaseSynId, seq_assay_id):
        bed = self._process(gene, seq_assay_id, newPath, parentId)
        process_functions.updateData(self.syn, databaseSynId, bed, seq_assay_id, filterByColumn="SEQ_ASSAY_ID", toDelete=True)
        bed.to_csv(newPath, sep="\t",index=False)
        return(newPath)

    def _validate(self, bed):
        total_error = ""
        warning = ""
        newCols = ["Chromosome","Start_Position","End_Position","Hugo_Symbol"]
        if len(bed.columns) < len(newCols):
            total_error += "Your BED file must at least have four columns in this order: %s.  Make sure there are no headers in your BED file.\n" % ", ".join(newCols)
        else:
            newCols.extend(range(0,len(bed.columns) - len(newCols)))
            bed.columns = newCols
            toValidateSymbol = True
            if not all(bed['Start_Position'].apply(lambda x: isinstance(x, int))):
                total_error += "The Start_Position column must only be integers. Make sure there are no headers in your BED file.\n"
                toValidateSymbol = False
            if not all(bed['End_Position'].apply(lambda x: isinstance(x, int))):
                total_error += "The End_Position column must only be integers. Make sure there are no headers in your BED file.\n"
                toValidateSymbol = False

            logger.info("VALIDATING GENE SYMBOLS")
            if any(bed['Hugo_Symbol'].isnull()):
                total_error += "You cannot submit any null symbols.\n"
            bed = bed[~bed['Hugo_Symbol'].isnull()]
            bed["Hugo_Symbol"] = [str(hugo).split(";")[0].split("_")[0].split(":")[0] for hugo in bed["Hugo_Symbol"]]
            if sum(bed['Hugo_Symbol'] == "+") != 0 or sum(bed['Hugo_Symbol'] == "-") !=0:
                total_error += "Fourth column must be the Hugo_Symbol column, not the strand column\n"

            if toValidateSymbol:
                genePosition = self.syn.tableQuery('SELECT * FROM syn11806563')
                genePositionDf = genePosition.asDataFrame()
                bed = bed.apply(lambda x: validateSymbol(x, genePositionDf, returnMappedDf=True), axis=1)
                if any(bed['Hugo_Symbol'].isnull()):
                    warning += "Any gene names that can't be remapped will be null.\n"
                if all(bed['Hugo_Symbol'].isnull()):
                    total_error += "You have no correct gene symbols. Make sure your gene symbol column (4th column) is formatted like so: SYMBOL(;optionaltext).  Optional text can be semi-colon separated.\n"
        return(total_error, warning)