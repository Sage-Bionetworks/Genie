import os
import logging
import subprocess

import pandas as pd
import synapseclient
try:
    from synapseclient.core.exceptions import SynapseHTTPError
except ModuleNotFoundError:
    from synapseclient.exceptions import SynapseHTTPError

from .example_filetype_format import FileTypeFormat
from . import process_functions

logger = logging.getLogger(__name__)

class maf(FileTypeFormat):
    '''
    MAF file format validation / processing
    '''
    _fileType = "maf"

    _process_kwargs = [
        "processing", "path_to_GENIE", 'databaseToSynIdMappingDf',
        "vcf2mafPath", "veppath", "vepdata", 'reference']

    def _validateFilename(self, filePath):
        '''
        Validates filename.  Should be
        data_mutations_extended_CENTER.txt
        '''
        assert os.path.basename(filePath[0]) == \
            "data_mutations_extended_{}.txt".format(self.center)

    def formatMAF(self, mafDf):
        '''
        Format maf file, shortens the maf file length
        '''
        keepMafColumns = [
            'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',
            'Chromosome', 'Start_Position', 'End_Position', 'Strand',
            'Variant_Classification', 'Variant_Type', 'Reference_Allele',
            'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
            'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
            'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1',
            'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
            'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1',
            'Match_Norm_Validation_Allele2', 'Verification_Status',
            'Validation_Status', 'Mutation_Status', 'Sequencing_Phase',
            'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File',
            'Sequencer', 'HGVSp_Short', 't_ref_count', 't_alt_count',
            'n_ref_count', 'n_alt_count', 'Protein_position', 'Codons',
            'SWISSPROT', 'RefSeq', 't_depth', 'n_depth', 'FILTER',
            'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF',
            'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF',
            'gnomAD_SAS_AF']

        mafDf = mafDf[keepMafColumns]
        mafDf['Center'] = self.center
        mafDf['Tumor_Sample_Barcode'] = [
            process_functions.checkGenieId(i, self.center)
            for i in mafDf['Tumor_Sample_Barcode']]
        mafDf['Sequence_Source'] = float('nan')
        mafDf['Sequencer'] = float('nan')
        mafDf['Validation_Status'][
            mafDf['Validation_Status'].isin(["Unknown", "unknown"])] = ''
        return(mafDf)

    def createFinalMaf(self, mafDf, filePath, maf=False):
        '''
        Creates the maf file at the end by either
        creating the first file or appending to the existing file

        Args:
            mafDf:  Maf dataframe
            filePath: Filepath to append or create
            maf: If maf, always create the file.  Default to False.
        '''
        if not mafDf.empty:
            if os.stat(filePath).st_size == 0 or maf:
                mafset = mafDf.to_csv(sep="\t", index=False)
            else:
                mafset = mafDf.to_csv(sep="\t", index=False, header=None)
            write_or_append = "wb" if maf else "ab"
            with open(filePath, write_or_append) as maffile:
                mafSet = process_functions.removeStringFloat(mafset)
                maffile.write(mafSet.encode("utf-8"))

    def storeProcessedMaf(
            self, filePath, mafSynId, centerMafSynId, isNarrow=False):
        '''
        Stores the processed maf
        There is a isNarrow option, but note that the number of rows
        of the maf file DOES NOT change in this function

        Args:
            filePath: Path to maf file
            mafSynId: database synid
            centerMafSynid: center flat file folder synid
            isNarrow: Is the file a narrow maf. Defaul to False.
        '''
        logger.info('STORING %s' % filePath)
        database = self.syn.get(mafSynId)
        if isNarrow:
            try:
                update_table = synapseclient.Table(
                    database.id, filePath, separator="\t")
                self.syn.store(update_table)
            except SynapseTimeoutError:
                # This error occurs because of waiting for table to index.
                # Don't worry about this.
                pass
        else:
            self.syn.store(
                synapseclient.File(filePath, parentId=centerMafSynId))
        return(filePath)

    def process_steps(
            self, filePath, path_to_GENIE, databaseToSynIdMappingDf,
            vcf2mafPath, veppath, vepdata, processing, reference=None):
        '''
        Processing maf files
        '''
        if processing == self._fileType:
            mafProcessing = "mafSP" if self._fileType == "mafSP" else 'vcf2maf'
            mafSynId = databaseToSynIdMappingDf.Id[
                databaseToSynIdMappingDf['Database'] == mafProcessing][0]
            centerMafSynId = databaseToSynIdMappingDf.Id[
                databaseToSynIdMappingDf['Database'] == "centerMaf"][0]

            logger.info('MAF2MAF %s' % filePath)
            fileName = "data_mutations_extended_%s_MAF.txt" % self.center
            newMafPath = os.path.join(
                path_to_GENIE, self.center, "staging", fileName)
            narrowMafPath = os.path.join(
                path_to_GENIE, self.center, "staging",
                "data_mutations_extended_{}_MAF_narrow.txt".format(
                    self.center))
            narrowMafColumns = [
                col['name'] for col in self.syn.getTableColumns(mafSynId)
                if col['name'] != 'inBED']
            # Strips out windows indentations \r
            command = ['dos2unix', filePath]
            subprocess.check_call(command)
            tempdir = os.path.join(path_to_GENIE, self.center)
            commandCall = ["perl", os.path.join(vcf2mafPath, "maf2maf.pl"),
                           "--input-maf", filePath,
                           "--output-maf", newMafPath,
                           "--vep-fork", '8',
                           "--tmp-dir", tempdir,
                           '--vep-path', veppath,
                           '--vep-data', vepdata,
                           "--custom-enst", os.path.join(
                                vcf2mafPath,
                                "data/isoform_overrides_uniprot")]
            if reference is not None:
                # '--ref-fasta','/root/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa',
                commandCall.extend(["--ref-fasta", reference])
            subprocess.check_call(commandCall)

            process_functions.rmFiles(tempdir, recursive=False)
            open(narrowMafPath, "w").close()
            if os.path.exists(newMafPath):
                # This needs to switch to streaming at some point
                mafDf = pd.read_csv(newMafPath, sep="\t", comment="#")
                mafDf.drop_duplicates(inplace=True)
                mafDf = self.formatMAF(mafDf)
                self.createFinalMaf(mafDf, newMafPath, maf=True)
                narrowMafDf = mafDf[narrowMafColumns]
                self.createFinalMaf(narrowMafDf, narrowMafPath, maf=True)
                # These functions have to be next to each other,
                # because no modifications can happen
                # Store Narrow MAF into db
                if self._fileType == "maf":
                    self.storeProcessedMaf(
                        narrowMafPath, mafSynId, centerMafSynId, isNarrow=True)
                # Store MAF flat file into synapse
                self.storeProcessedMaf(newMafPath, mafSynId, centerMafSynId)
            else:
                logger.error('ERROR PROCESSING %s' % filePath)
                filePath = "NOTPROCESSED"
        else:
            logger.info(
                "Please run with `--process {filetype}` parameter if you want "
                "to reannotate the {filetype} files".format(
                    filetype=self._fileType))
        return(filePath)

    def _validate(self, mutationDF):
        """
        This function validates the mutation file to make sure it
        adheres to the mutation SOP.

        Args:
            mutationDF: mutation dataframe

        Returns:
            Text with all the errors in the mutation file
        """

        first_header = ['CHROMOSOME', 'HUGO_SYMBOL', 'TUMOR_SAMPLE_BARCODE']
        SP = self._fileType == "mafSP"
        if SP:
            correct_column_headers = [
                'CHROMOSOME', 'START_POSITION', 'REFERENCE_ALLELE',
                'TUMOR_SAMPLE_BARCODE', 'TUMOR_SEQ_ALLELE2']
            # T_REF_COUNT + T_ALT_COUNT = T_DEPTH
        else:
            correct_column_headers = [
                'CHROMOSOME', 'START_POSITION', 'REFERENCE_ALLELE',
                'TUMOR_SAMPLE_BARCODE', 'T_ALT_COUNT', 'TUMOR_SEQ_ALLELE2']
            # T_REF_COUNT + T_ALT_COUNT = T_DEPTH
        optional_headers = \
            ['T_REF_COUNT', 'N_DEPTH', 'N_REF_COUNT', 'N_ALT_COUNT']

        mutationDF.columns = [col.upper() for col in mutationDF.columns]

        total_error = ""
        warning = ""

        # CHECK: Everything in correct_column_headers must be in mutation file
        if not all([process_functions.checkColExist(mutationDF, i)
                    for i in correct_column_headers]):
            total_error += (
                "Mutation File: "
                "Must at least have these headers: {}. "
                "If you are writing your maf file with R, please make"
                "sure to specify the 'quote=FALSE' parameter.\n".format(
                    ",".join([i for i in correct_column_headers
                              if i not in mutationDF.columns.values])))
        else:
            # CHECK: First column must be in the first_header list
            if mutationDF.columns[0] not in first_header:
                total_error += ("Mutation File: First column header must be "
                                "one of these: {}.\n".format(
                                    ", ".join(first_header)))
            # No duplicated values
            primary_cols = ['CHROMOSOME', 'START_POSITION',
                            'REFERENCE_ALLELE', 'TUMOR_SAMPLE_BARCODE',
                            'TUMOR_SEQ_ALLELE2']
            if mutationDF.duplicated(primary_cols).any():
                total_error += ("Mutation File: "
                                "Should not have duplicate rows\n")

        check_col = process_functions.checkColExist(mutationDF, "T_DEPTH")
        if not check_col and not SP:
            if not process_functions.checkColExist(mutationDF, "T_REF_COUNT"):
                total_error += (
                    "Mutation File: "
                    "If you are missing T_DEPTH, you must have T_REF_COUNT!\n")

        # CHECK: Must have either TUMOR_SEQ_ALLELE2 column
        if process_functions.checkColExist(mutationDF, "TUMOR_SEQ_ALLELE2"):
            # CHECK: The value "NA" can't be used as a placeholder
            if sum(mutationDF["TUMOR_SEQ_ALLELE2"].fillna('') == "NA") > 0:
                warning += (
                    "Mutation File: "
                    "TUMOR_SEQ_ALLELE2 column contains 'NA' values, "
                    "which cannot be placeholders for blank values.  "
                    "Please put in empty strings for blank values.\n")
            # CHECK: There can't be any null values
            if sum(mutationDF["TUMOR_SEQ_ALLELE2"].isnull()) > 0:
                total_error += (
                    "Mutation File: "
                    "TUMOR_SEQ_ALLELE2 can't have any null values.\n")

        # CHECK: Mutation file would benefit from columns in optional_headers
        if not all([
                process_functions.checkColExist(mutationDF, i)
                for i in optional_headers]) and not SP:
            warning += (
                "Mutation File: "
                "Does not have the column headers that can give extra "
                "information to the processed mutation file: {}.\n".format(
                    ", ".join([
                        i for i in optional_headers
                        if i not in mutationDF.columns.values])))

        if process_functions.checkColExist(mutationDF, "REFERENCE_ALLELE"):
            if sum(mutationDF['REFERENCE_ALLELE'] == "NA") > 0:
                warning += (
                    "Mutation File: "
                    "Your REFERENCE_ALLELE column contains NA values, "
                    "which cannot be placeholders for blank values.  "
                    "Please put in empty strings for blank values.\n")
            # CHECK: mutation file must not have empty reference
            # or variant alleles
            if sum(mutationDF['REFERENCE_ALLELE'].isnull()) > 0:
                total_error += (
                    "Mutation File: "
                    "Cannot have any empty REFERENCE_ALLELE values.\n")

        if process_functions.checkColExist(mutationDF, "CHROMOSOME"):
            # CHECK: Chromosome column can't have any values that start
            # with chr or have any WT values
            invalidValues = [
                str(i).startswith("chr") or str(i) == "WT"
                for i in mutationDF['CHROMOSOME']]
            if sum(invalidValues) > 0:
                total_error += (
                    "Mutation File: "
                    "CHROMOSOME column cannot have any values that "
                    "start with 'chr' or any 'WT' values.\n")

        return(total_error, warning)

    def _get_dataframe(self, filePathList):
        '''
        Get mutation dataframe
        '''
        mutationdf = pd.read_csv(filePathList[0],
                                 sep="\t",
                                 comment="#",
                                 na_values=['-1.#IND', '1.#QNAN', '1.#IND',
                                            '-1.#QNAN', '#N/A N/A', 'NaN',
                                            '#N/A', 'N/A', '#NA', 'NULL',
                                            '-NaN', 'nan', '-nan', ''],
                                 keep_default_na=False,
                                 # This is to check if people write files
                                 # with R, quote=T
                                 quoting=3)
        return mutationdf
