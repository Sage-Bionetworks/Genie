"""GENIE bed class and functions"""
import os
import logging
import subprocess

import pandas as pd

from .example_filetype_format import FileTypeFormat
from . import process_functions

LOGGER = logging.getLogger(__name__)


# def createGenePositionsTables():
#     '''
#     Create gene position database
#     '''
#     from biomart import BiomartServer
#     import synapseclient
#     import time
#     syn = synapseclient.login()
#     print(time.time())
#     server = BiomartServer("grch37.ensembl.org/biomart")

#     hsapiens_gene_ensembl = server.datasets['hsapiens_gene_ensembl']
#     # hsapiens_gene_ensembl.show_attributes()

#     chromosomes = range(1, 23)
#     chromosomes.extend(["X", "Y"])
#     chromosomes = [str(i) for i in chromosomes]
#     response = hsapiens_gene_ensembl.search({
#         'attributes': [
#             'hgnc_symbol',
#             "chromosome_name",
#             "start_position",
#             "end_position"]
#     })
#     genes = pd.DataFrame()
#     for line in response.iter_lines():
#         line = line.decode('utf-8')
#         row = line.split("\t")
#         if row[0] != '' and row[1] in chromosomes:
#             tempRow = pd.DataFrame([row])
#             genes = genes.append(tempRow)
#     genes.columns = [
#         'hgnc_symbol', "chromosome_name",
#         "start_position", "end_position"]
#     # genes.to_csv("genepositions.csv",index=False)
#     print(time.time())
#     databaseSynId = "syn11806563"
#     databaseEnt = syn.get(databaseSynId)
#     database = syn.tableQuery("SELECT * FROM %s" % databaseSynId)
#     database = database.asDataFrame()
#     process_functions.updateDatabase(
#         syn, database, genes, databaseSynId,
#         databaseEnt.primaryKey, toDelete=True)
#     return(genes)

def create_gtf(dirname):
    """
    Create exon.gtf and gene.gtf from GRCh37 gtf

    Args:
        dirname: Directory where these files should live

    Returns:
        exon_gtf_path: exon GTF
        gene_gtf_path: gene GTF
    """
    exon_gtf_path = os.path.join(dirname, "exon.gtf")
    gene_gtf_path = os.path.join(dirname, "gene.gtf")

    if not os.path.exists(exon_gtf_path) or not os.path.exists(gene_gtf_path):
        download_cmd = ['wget',
                        'http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz',
                        '-P', dirname]
        subprocess.check_call(download_cmd)
        gtfgz_path = os.path.join(dirname, "Homo_sapiens.GRCh37.75.gtf.gz")
        gunzip_cmd = ['gunzip', '-f', gtfgz_path]
        subprocess.check_call(gunzip_cmd)
        gtf_path = os.path.join(dirname, "Homo_sapiens.GRCh37.75.gtf")

        exon_awk_cmd = ['awk', '$3 == "exon" {print}', gtf_path]
        exon_gtf = subprocess.check_output(
            exon_awk_cmd, universal_newlines=True)
        with open(exon_gtf_path, "w") as gtf_file:
            gtf_file.write(exon_gtf)
        gene_awk_cmd = ['awk', '$3 == "gene" {print}', gtf_path]
        gene_gtf = subprocess.check_output(
            gene_awk_cmd, universal_newlines=True)
        with open(gene_gtf_path, "w") as gtf_file:
            gtf_file.write(gene_gtf)
    return(exon_gtf_path, gene_gtf_path)


def _add_feature_type_todf(filepath, featuretype):
    """
    Add Feature_Type to dataframe

    Args:
        filepath: path to bed
        featuretype: exon, intron, or intergenic

    Returns:
        df: empty dataframe or dataframe with appended feature type
    """
    # No need to add anything if the dataframe is empty
    if os.stat(filepath).st_size != 0:
        beddf = pd.read_csv(filepath, sep="\t", header=None)
        beddf.columns = ["Chromosome", "Start_Position", "End_Position",
                         "Hugo_Symbol", "includeInPanel", "clinicalReported",
                         "ID", "SEQ_ASSAY_ID"]
        beddf['Feature_Type'] = featuretype
    else:
        beddf = pd.DataFrame()
    return beddf


def add_feature_type(temp_bed_path, exon_gtf_path, gene_gtf_path):
    """
    Add Feature_Type to bed file (exon, intron, intergenic)

    Args:
        temp_bed_path: BED file without feature type
        exon_gtf_path: exon gtf
        gene_gtf_path: gene gtf

    Returns:
        genie_combined_path: Path to final bed file
    """
    genie_exon_path = os.path.join(process_functions.SCRIPT_DIR,
                                   'genie_exons.bed')
    genie_intron_path = os.path.join(process_functions.SCRIPT_DIR,
                                     'genie_introns.bed')
    genie_intergenic_path = os.path.join(process_functions.SCRIPT_DIR,
                                         'genie_intergenic.bed')
    intron_intergenic_path = os.path.join(process_functions.SCRIPT_DIR,
                                          'intron_intergenic.bed')
    gene_path = os.path.join(process_functions.SCRIPT_DIR, 'gene.bed')
    # GET EXON REGIONS
    # Get intersection between true exon regions and submitted bed regions
    command = ['bedtools', 'intersect', '-a',
               temp_bed_path, '-b', exon_gtf_path, '-wa',
               '|', 'sort', '|', 'uniq', '>', genie_exon_path]
    subprocess.check_call(" ".join(command), shell=True)
    # get intergenic/intron regions
    # Get opposite of intersection between true exon regions and submitted
    # bed regions
    command = ['bedtools', 'intersect', '-a',
               temp_bed_path, '-b', exon_gtf_path, '-wa', '-v'
               '|', 'sort', '|', 'uniq', '>', intron_intergenic_path]
    subprocess.check_call(" ".join(command), shell=True)
    # get gene regions
    # Get intersection between true gene regions and submitted bed regions
    command = ['bedtools', 'intersect', '-a',
               temp_bed_path, '-b', gene_gtf_path, '-wa',
               '|', 'sort', '|', 'uniq', '>', gene_path]
    subprocess.check_call(" ".join(command), shell=True)
    # GET INTRON REGSIONS
    # Difference between the gene regions and exon regions will give
    # intron regions
    command = ['diff', gene_path, genie_exon_path, '|',
               'grep', "'<'", '|', 'sed', "'s/< //'", '>',
               genie_intron_path]
    subprocess.check_call(" ".join(command), shell=True)
    # GET INTERGENIC REGIONS
    # Difference between the intron/intergenic and intron regions will
    # give intergenic regions
    command = ['diff', intron_intergenic_path, genie_intron_path, '|',
               'grep', "'<'", '|', 'sed', "'s/< //'", '>',
               genie_intergenic_path]
    subprocess.check_call(" ".join(command), shell=True)

    genie_exondf = _add_feature_type_todf(genie_exon_path, "exon")
    genie_introndf = _add_feature_type_todf(genie_intron_path, "intron")
    genie_intergenicdf = _add_feature_type_todf(genie_intergenic_path,
                                                "intergenic")
    # Specify the combined df in case there no bed hits at all
    genie_combineddf = pd.DataFrame(columns=["Chromosome", "Start_Position",
                                             "End_Position", "Hugo_Symbol",
                                             "includeInPanel",
                                             "clinicalReported",
                                             "ID", "SEQ_ASSAY_ID",
                                             "Feature_Type"])
    genie_combineddf = genie_combineddf.append(genie_exondf)
    genie_combineddf = genie_combineddf.append(genie_introndf)
    genie_combineddf = genie_combineddf.append(genie_intergenicdf)
    return genie_combineddf


def _check_region_overlap(row, gene_positiondf):
    """
    Check if the submitted bed symbol + region overlaps with
    the actual gene's positions

    Args:
        row: row in bed file (genomic region)
        gene_positiondf: Actual gene position dataframe

    Return:
        True if the region does overlap.
    """
    matching_symbol_ind = gene_positiondf['hgnc_symbol'] == row['Hugo_Symbol']
    match_with_genedb = gene_positiondf[matching_symbol_ind]

    if not match_with_genedb.empty:
        # We are assuming that there is only one matching gene, but
        # there are actually duplicated gene symbols
        start_position = match_with_genedb['start_position'].values[0]
        end_position = match_with_genedb['end_position'].values[0]
        # Submitted bed start position in a gene
        start_in_gene = start_position <= row['Start_Position'] <= end_position
        # Submitted bed end position in a gene
        end_in_gene = start_position <= row['End_Position'] <= end_position
        # Submitted bed region surrounds a gene
        gene_in_region = end_position <= row['End_Position'] and \
                         start_position >= row['Start_Position']
        return start_in_gene or end_in_gene or gene_in_region
    return False


def _get_max_overlap_index(overlap, bed_length, boundary):
    """
    Calculate the ratio of overlap between the submitted bed
    region and gene position dataframe and return the index of
    the max overlapping region

    Args:
        overlap: Possible overlapping region
        bed_length: Length of submitted region
        boundary: specified ratio overlap

    Returns:
        Index of regions with maximum overlap or None
    """
    ratio_overlap = overlap / bed_length
    ratio_overlap = ratio_overlap[ratio_overlap > boundary]
    ratio_overlap = ratio_overlap[ratio_overlap <= 1]
    if not ratio_overlap.empty:
        return ratio_overlap.idxmax()
    return None


def _map_gene_within_boundary(row, gene_positiondf, boundary=0.9):
    """
    Map gene as long as the strand is a percent within the boundary
    of the start goes over start boundary, but end is contained in gene

    Args:
        row: row in bed file (genomic region)
        genePositionDf: Actual gene position dataframe
        boundary: Percent boundary defined

    Return:
        pd.Series: mapped gene
    """
    # First filter just based on whether or not region is within the gene
    # position dataframe
    match_chr = gene_positiondf['chromosome_name'] == str(row['Chromosome'])
    chrom_rows = gene_positiondf[match_chr]
    match_start = chrom_rows['start_position'] <= row['Start_Position']
    start_rows = chrom_rows[match_start]
    end_rows = start_rows[start_rows['end_position'] >= row['End_Position']]

    bed_length = row['End_Position'] - row['Start_Position']
    # as long as the strand is 90% within the boundary of the
    # Start goes over start boundary, but end is contained in gene
    if end_rows.empty:
        if sum(chrom_rows['end_position'] >= row['End_Position']) > 0:
            overlap = row['End_Position'] - chrom_rows['start_position']
            # difference =  difference * -1.0
            max_overlap = _get_max_overlap_index(overlap, bed_length,
                                                 boundary)
            if max_overlap is not None:
                end_rows = end_rows.append(chrom_rows.loc[[max_overlap]])
        # End goes over end boundary, but start is contained in gene
        if sum(chrom_rows['start_position'] <= row['Start_Position']) > 0:
            overlap = chrom_rows['end_position'] - row['Start_Position']
            max_overlap = _get_max_overlap_index(overlap, bed_length,
                                                 boundary)
            if max_overlap is not None:
                end_rows = end_rows.append(chrom_rows.loc[[max_overlap]])
        # Start and end go over gene boundary
        check = chrom_rows[chrom_rows['start_position'] >= row['Start_Position']]
        check = check[check['end_position'] <= row['End_Position']]
        if not check.empty:
            overlap = chrom_rows['end_position'] - chrom_rows['start_position']
            max_overlap = _get_max_overlap_index(overlap, bed_length,
                                                 boundary)
            if max_overlap is not None:
                end_rows = end_rows.append(chrom_rows.loc[[max_overlap]])
    return end_rows


def validate_symbol(row, gene_positiondf, return_mappeddf=False):
    """
    The apply function of a DataFrame is called twice on the first row (known
    pandas behavior)

    Args:
        row: start and end position
        gene_positiondf: Actual gene position dataframe
        return_mappeddf: Return mapped dataframe. Default is False.

    Return:
        bool or Series: if the gene passed in need to be remapped or
                        the remapped gene
    """
    valid = True
    region_overlap = _check_region_overlap(row, gene_positiondf)
    if not region_overlap:
        end_rows = _map_gene_within_boundary(row, gene_positiondf)
        if len(end_rows) == 0:
            LOGGER.warning("{} cannot be remapped. "
                           "These rows will have an empty gene symbol".format(
                                row['Hugo_Symbol']))
            row['Hugo_Symbol'] = pd.np.nan
            valid = False
        elif len(end_rows) > 1:
            if row['Hugo_Symbol'] not in end_rows['hgnc_symbol'].tolist():
                # if "MLL4", then the HUGO symbol should be KMT2D and KMT2B
                LOGGER.warning("{} can be mapped to different symbols: {}. "
                               "Please correct or it will be removed.".format(
                                    row['Hugo_Symbol'],
                                    ", ".join(end_rows['hgnc_symbol'])))
                row['Hugo_Symbol'] = pd.np.nan
                valid = False
        else:
            if row['Hugo_Symbol'] != end_rows['hgnc_symbol'].values[0]:
                LOGGER.info("{} will be remapped to {}".format(
                                row['Hugo_Symbol'],
                                end_rows['hgnc_symbol'].values[0]))
                row['Hugo_Symbol'] = end_rows['hgnc_symbol'].values[0]

    if return_mappeddf:
        return row
    return valid


class bed(FileTypeFormat):
    """GENIE bed format"""
    _fileType = "bed"

    _process_kwargs = ["newPath", "parentId", "databaseSynId", 'seq_assay_id']

    def _get_dataframe(self, filePathList):
        '''
        Bed files don't have a header

        Args:
            filePathList: List of files
        '''
        filePath = filePathList[0]
        try:
            beddf = pd.read_csv(filePath, sep="\t", header=None)
        except Exception:
            raise ValueError(
                "Can't read in your bed file. "
                "Please make sure the BED file is not binary and "
                "does not contain a comment/header line")
        if not str(beddf[0][0]).isdigit() and \
           not str(beddf[0][0]).startswith("chr"):
            raise ValueError(
                "Please make sure your bed file does not "
                "contain a comment/header line")
        return(beddf)

    def _validateFilename(self, filePath):
        '''
        Validates filename
        CENTER-11.bed

        Args:
            filePath: Path to bedfile
        '''

        assert os.path.basename(filePath[0]).startswith("%s-" % self.center)\
            and os.path.basename(filePath[0]).endswith(".bed")

    def create_gene_panel(self, temp_bed_path, seq_assay_id,
                          genePanelPath, parentId,
                          exon_gtf_path,
                          createGenePanel=True):
        '''
        Create bed file and gene panel files from the bed file

        Args:
            temp_bed_path: Temporary bed file
            seq_assay_id: GENIE SEQ_ASSAY_ID
            genePanelPath: Gene panel folder path
            parentId: Synapse id of gene panel folder
            exon_gtf_path: Exon GTF path, created by create_gtf()
            createGenePanel: To create data gene panel files. Default is True

        Returns:
            pd.DataFrame: configured bed dataframe
        '''
        LOGGER.info("CREATING GENE PANEL")
        command = ['bedtools', 'intersect', '-a',
                   temp_bed_path,
                   '-b', exon_gtf_path,
                   '-u']
        # Create GENIE genie_exons.bed for gene panel file
        genie_exon_text = subprocess.check_output(command,
                                                  universal_newlines=True)
        genie_exon_path = os.path.join(process_functions.SCRIPT_DIR,
                                       'genie_exons.bed')
        with open(genie_exon_path, "w") as genie_exon:
            genie_exon.write(genie_exon_text)

        genie_exon_stat = os.stat(genie_exon_path)
        if genie_exon_stat.st_size > 0 and createGenePanel:
            exonsdf = pd.read_csv(genie_exon_path, sep="\t", header=None)
            exonsdf.columns = ["Chromosome", "Start_Position", "End_Position",
                               "Hugo_Symbol", "includeInPanel",
                               "clinicalReported",
                               "ID", "SEQ_ASSAY_ID"]
            # Only include genes that should be included in the panels
            include_exonsdf = exonsdf[exonsdf['includeInPanel']]
            # Write gene panel
            null_genes = include_exonsdf['Hugo_Symbol'].isnull()
            unique_genes = set(include_exonsdf['Hugo_Symbol'][~null_genes])
            gene_panel_text = ("stable_id: {seq_assay_id}\n"
                               "description: {seq_assay_id}, "
                               "Number of Genes - {num_genes}\n"
                               "gene_list:\t{genelist}".format(
                                    seq_assay_id=seq_assay_id,
                                    num_genes=len(unique_genes),
                                    genelist="\t".join(unique_genes)))
            genepanelname = "data_gene_panel_" + seq_assay_id + ".txt"

            with open(os.path.join(genePanelPath, genepanelname), "w+") as f:
                f.write(gene_panel_text)
            process_functions.storeFile(self.syn,
                                        os.path.join(genePanelPath,
                                                     genepanelname),
                                        parentId=parentId,
                                        center=self.center,
                                        fileFormat="bed",
                                        dataSubType="metadata",
                                        cBioFileFormat="genePanel")

    def _process(self, gene, seq_assay_id, newPath,
                 parentId, createPanel=True):
        '''
        Process bed file, add feature type

        Args:
            gene: bed dataframe
            seq_assay_id: GENIE SEQ_ASSAY_ID
            newPath: new GENIE path
            parentId: Synapse id to store the gene panel
            createPanel: Create gene panel

        Returns:
            pd.DataFrame: Conigured bed dataframe
        '''
        seq_assay_id = seq_assay_id.upper()
        seq_assay_id = seq_assay_id.replace('_', '-')

        # Add in 6th column which is the clinicalReported
        if len(gene.columns) > 5:
            if all(gene[5].apply(lambda x: x in [True, False])):
                gene[5] = gene[5].astype(bool)
            else:
                gene[5] = pd.np.nan
        else:
            gene[5] = pd.np.nan
        bed = gene[[0, 1, 2, 3, 4, 5]]
        genePanelPath = os.path.dirname(newPath)
        # Must be .astype(bool) because `1, 0 in [True, False]`
        bed[4] = bed[4].astype(bool)

        exon_gtf_path, gene_gtf_path = create_gtf(process_functions.SCRIPT_DIR)
        LOGGER.info("REMAPPING %s" % seq_assay_id)
        # bedname = seq_assay_id + ".bed"
        bed.columns = ["Chromosome", "Start_Position", "End_Position",
                       "Hugo_Symbol", "includeInPanel", "clinicalReported"]
        # Validate gene symbols
        # Gene symbols can be split by ; and _ and : and .
        bed['Hugo_Symbol'] = [
            i.split(";")[0].split("_")[0].split(":")[0].split(".")[0]
            for i in bed['Hugo_Symbol']]
        # Replace all chr with blank
        bed['Chromosome'] = [
            str(i).replace("chr", "") for i in bed['Chromosome']]
        # Change all start and end to int
        bed['Start_Position'] = bed['Start_Position'].apply(int)
        bed['End_Position'] = bed['End_Position'].apply(int)

        genePosition = self.syn.tableQuery('SELECT * FROM syn11806563')
        genePositionDf = genePosition.asDataFrame()
        bed['ID'] = bed['Hugo_Symbol']
        bed = bed.apply(lambda x: validate_symbol(x, genePositionDf,
                                                  return_mappeddf=True), axis=1)
        bed['SEQ_ASSAY_ID'] = seq_assay_id
        temp_bed_path = os.path.join(process_functions.SCRIPT_DIR, "temp.bed")
        bed.to_csv(temp_bed_path, sep="\t", index=False, header=None)
        self.create_gene_panel(temp_bed_path, seq_assay_id,
                               genePanelPath, parentId, exon_gtf_path,
                               createGenePanel=createPanel)
        final_bed = add_feature_type(temp_bed_path, exon_gtf_path,
                                     gene_gtf_path)

        final_bed['CENTER'] = self.center
        final_bed['Chromosome'] = final_bed['Chromosome'].astype(str)
        return final_bed

    def preprocess(self, filePath):
        '''
        Standardize and grab seq assay id from the bed file path

        Args:
            filePath: bed file path

        Returns:
            dict: GENIE seq assay id
        '''
        seq_assay_id = os.path.basename(filePath).replace(".bed", "")
        seq_assay_id = seq_assay_id.upper().replace("_", "-")
        return {'seq_assay_id': seq_assay_id}

    def process_steps(self, gene, newPath, parentId, databaseSynId,
                      seq_assay_id):
        '''
        Process bed file, update bed database, write bed file to path

        Args:
            gene: Bed dataframe
            newPath: Path to new bed file
            parentId: Synapse id to store gene panel file
            databaseSynId: Synapse id of bed database
            seq_assay_id: GENIE seq assay id

        Returns:
            string: Path to new bed file
        '''
        bed = self._process(gene, seq_assay_id, newPath, parentId)
        process_functions.updateData(self.syn, databaseSynId, bed,
                                     seq_assay_id,
                                     filterByColumn="SEQ_ASSAY_ID",
                                     toDelete=True)
        bed.to_csv(newPath, sep="\t", index=False)
        return(newPath)

    def _validate(self, beddf):
        '''
        Validate bed file

        Args:
            bed: Bed dataframe

        Returns:
            total_error: all the errors
            warning: all the warnings
        '''
        total_error = ""
        warning = ""
        newCols = ["Chromosome", "Start_Position",
                   "End_Position", "Hugo_Symbol",
                   "includeInPanel"]
        if len(beddf.columns) < len(newCols):
            total_error += (
                "BED file: Must at least have five columns in this "
                "order: {}. Make sure there are "
                "no headers.\n".format(", ".join(newCols)))
        else:
            newCols.extend(range(0, len(beddf.columns) - len(newCols)))
            beddf.columns = newCols
            toValidateSymbol = True
            if not all(beddf['Start_Position'].apply(
                    lambda x: isinstance(x, int))):
                total_error += (
                    "BED file: "
                    "The Start_Position column must only be integers. "
                    "Make sure there are no headers.\n")
                toValidateSymbol = False
            if not all(beddf['End_Position'].apply(
                    lambda x: isinstance(x, int))):
                total_error += (
                    "BED file: "
                    "The End_Position column must only be integers. "
                    "Make sure there are no headers.\n")
                toValidateSymbol = False

            LOGGER.info("VALIDATING GENE SYMBOLS")
            if any(beddf['Hugo_Symbol'].isnull()):
                total_error += \
                    "BED file: You cannot submit any null symbols.\n"
            beddf = beddf[~beddf['Hugo_Symbol'].isnull()]
            beddf["Hugo_Symbol"] = [
                str(hugo).split(";")[0].split("_")[0].split(":")[0]
                for hugo in beddf["Hugo_Symbol"]]
            if sum(beddf['Hugo_Symbol'] == "+") != 0 or \
               sum(beddf['Hugo_Symbol'] == "-") != 0:
                total_error += (
                    "BED file: Fourth column must be the Hugo_Symbol column, "
                    "not the strand column\n")

            warn, error = process_functions.check_col_and_values(
                beddf,
                'includeInPanel',
                [True, False],
                filename="BED file",
                required=True)
            warning += warn
            total_error += error

            if toValidateSymbol:
                genePosition = self.syn.tableQuery('SELECT * FROM syn11806563')
                genePositionDf = genePosition.asDataFrame()
                beddf = beddf.apply(
                    lambda x: validate_symbol(
                        x, genePositionDf, return_mappeddf=True), axis=1)

                if any(beddf['Hugo_Symbol'].isnull()):
                    warning += (
                        "BED file: "
                        "Any gene names that can't be "
                        "remapped will be null.\n")
                if all(beddf['Hugo_Symbol'].isnull()):
                    total_error += (
                        "BED file: "
                        "You have no correct gene symbols. "
                        "Make sure your gene symbol column (4th column) "
                        "is formatted like so: SYMBOL(;optionaltext).  "
                        "Optional text can be semi-colon separated.\n")
        return(total_error, warning)
