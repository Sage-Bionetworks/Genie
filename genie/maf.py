import os
import logging
import subprocess

import pandas as pd
import synapseclient
from synapseclient.core.exceptions import SynapseHTTPError

from .example_filetype_format import FileTypeFormat
from . import process_functions

logger = logging.getLogger(__name__)

class maf(FileTypeFormat):
    '''
    MAF file format validation / processing
    '''
    _fileType = "maf"

    _process_kwargs = []

    def _validateFilename(self, filePath):
        '''
        Validates filename.  Should be
        data_mutations_extended_CENTER.txt
        '''
        assert os.path.basename(filePath[0]) == \
            "data_mutations_extended_{}.txt".format(self.center)

    def process_steps(self, df):
        """The processing of maf files is specific to GENIE, so
        not included in this function"""
        logger.info("Please run with `--process mutation` parameter "
                    "if you want to reannotate the mutation files")
        return None

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
        """Get mutation dataframe"""
        # Must do this because pandas.read_csv will allow for a file to
        # have more column headers than content.  E.g.
        # A,B,C,D,E
        # 1,2
        # 2,3
        with open(filePathList[0], 'r') as maf_f:
            firstline = maf_f.readline()
            if firstline.startswith("#"):
                firstline = maf_f.readline()
            secondline = maf_f.readline()

        if len(firstline.split("\t")) != len(secondline.split("\t")):
            raise ValueError("Number of fields in a line do not match the "
                             "expected number of columns")

        mutationdf = pd.read_csv(filePathList[0],
                                 sep="\t",
                                 comment="#",
                                 # Keep the value 'NA'
                                 na_values=['-1.#IND', '1.#QNAN', '1.#IND',
                                            '-1.#QNAN', '#N/A N/A', 'NaN',
                                            '#N/A', 'N/A', '#NA', 'NULL',
                                            '-NaN', 'nan', '-nan', ''],
                                 keep_default_na=False,
                                 # This is to check if people write files
                                 # with R, quote=T
                                 quoting=3,
                                 # Retain completely blank lines so that
                                 # validator will cause the file to fail
                                 skip_blank_lines=False)
        return mutationdf
