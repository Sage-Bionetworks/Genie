from io import StringIO
import logging
import os
from typing import List

import pandas as pd

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions, transform, validate

logger = logging.getLogger(__name__)


def _check_tsa1_tsa2(df):
    """If maf file has both TSA1 and TSA2,
    TSA1 must equal REF, or TSA1 must equal TSA2.
    """
    tsa2_col_exist = process_functions.checkColExist(df, "TUMOR_SEQ_ALLELE2")
    tsa1_col_exist = process_functions.checkColExist(df, "TUMOR_SEQ_ALLELE1")
    ref_col_exist = process_functions.checkColExist(df, "REFERENCE_ALLELE")
    error = ""
    if tsa2_col_exist and tsa1_col_exist and ref_col_exist:
        tsa1_eq_ref = all(df["TUMOR_SEQ_ALLELE1"] == df["REFERENCE_ALLELE"])
        tsa1_eq_tsa2 = all(df["TUMOR_SEQ_ALLELE1"] == df["TUMOR_SEQ_ALLELE2"])
        if not (tsa1_eq_ref or tsa1_eq_tsa2):
            error = (
                "maf: Contains both "
                "TUMOR_SEQ_ALLELE1 and TUMOR_SEQ_ALLELE2 columns. "
                "All values in TUMOR_SEQ_ALLELE1 must match all values in "
                "REFERENCE_ALLELE or all values in TUMOR_SEQ_ALLELE2.\n"
            )
    return error


def _check_allele_col(df, col):
    """
    Check the Allele column is correctly formatted.

    Args:
        df: mutation dataframe
        col: Column header name

    Returns:
        error, warning

    """
    col_exist = process_functions.checkColExist(df, col)
    error = ""
    warning = ""
    if col_exist:
        # CHECK: There can't be any null values
        if sum(df[col].isnull()) > 0:
            error = f"maf: {col} can't have any blank or null values.\n"

    return error, warning


class maf(FileTypeFormat):
    """
    MAF file format validation / processing
    """

    _fileType = "maf"

    _process_kwargs = []
    _allele_cols = ["REFERENCE_ALLELE", "TUMOR_SEQ_ALLELE1", "TUMOR_SEQ_ALLELE2"]
    _allowed_comb_alleles = ["A", "T", "C", "G", "N"]
    _allowed_ind_alleles = ["-"]

    def _validateFilename(self, filePath):
        """
        Validates filename.  Should be
        data_mutations_extended_CENTER.txt
        """
        assert os.path.basename(filePath[0]) == "data_mutations_extended_{}.txt".format(
            self.center
        )

    def process_steps(self, df):
        """The processing of maf files is specific to GENIE, so
        not included in this function"""
        logger.info(
            "Please run with `--process mutation` parameter "
            "if you want to reannotate the mutation files"
        )
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

        first_header = ["CHROMOSOME", "HUGO_SYMBOL", "TUMOR_SAMPLE_BARCODE"]
        SP = self._fileType == "mafSP"
        if SP:
            correct_column_headers = [
                "CHROMOSOME",
                "START_POSITION",
                "REFERENCE_ALLELE",
                "TUMOR_SAMPLE_BARCODE",
                "TUMOR_SEQ_ALLELE2",
            ]
            # T_REF_COUNT + T_ALT_COUNT = T_DEPTH
        else:
            correct_column_headers = [
                "CHROMOSOME",
                "START_POSITION",
                "REFERENCE_ALLELE",
                "TUMOR_SAMPLE_BARCODE",
                "T_ALT_COUNT",
                "TUMOR_SEQ_ALLELE2",
            ]
            # T_REF_COUNT + T_ALT_COUNT = T_DEPTH
        optional_headers = ["T_REF_COUNT", "N_DEPTH", "N_REF_COUNT", "N_ALT_COUNT"]

        mutationDF.columns = [col.upper() for col in mutationDF.columns]

        # total_error = ""
        total_error = StringIO()
        warning = StringIO()

        # CHECK: Everything in correct_column_headers must be in mutation file
        if not all(
            [
                process_functions.checkColExist(mutationDF, i)
                for i in correct_column_headers
            ]
        ):
            total_error.write(
                "maf: Must at least have these headers: {}. "
                "If you are writing your maf file with R, please make"
                "sure to specify the 'quote=FALSE' parameter.\n".format(
                    ",".join(
                        [
                            i
                            for i in correct_column_headers
                            if i not in mutationDF.columns.values
                        ]
                    )
                )
            )
        else:
            # CHECK: First column must be in the first_header list
            if mutationDF.columns[0] not in first_header:
                total_error.write(
                    "maf: First column header must be "
                    "one of these: {}.\n".format(", ".join(first_header))
                )
            # No duplicated values
            primary_cols = [
                "CHROMOSOME",
                "START_POSITION",
                "REFERENCE_ALLELE",
                "TUMOR_SAMPLE_BARCODE",
                "TUMOR_SEQ_ALLELE2",
            ]
            # Strip white space if string column
            for col in primary_cols:
                if mutationDF[col].dtype == object:
                    mutationDF[col] = mutationDF[col].str.strip()
            duplicated_idx = mutationDF.duplicated(primary_cols)
            # Find samples with duplicated variants
            duplicated_variants = (
                mutationDF["TUMOR_SAMPLE_BARCODE"][duplicated_idx]
                .unique()
                .astype(str)
                .tolist()
            )

            if duplicated_idx.any():
                total_error.write(
                    "maf: Must not have duplicated variants. "
                    "Samples with duplicated variants: "
                    f"{', '.join(duplicated_variants)}\n"
                )

        t_depth_exists = process_functions.checkColExist(mutationDF, "T_DEPTH")
        t_ref_exists = process_functions.checkColExist(mutationDF, "T_REF_COUNT")
        if not t_depth_exists and not t_ref_exists and not SP:
            total_error.write("maf: If missing T_DEPTH, must have T_REF_COUNT!\n")
        numerical_cols = [
            "T_DEPTH",
            "T_ALT_COUNT",
            "T_REF_COUNT",
            "N_DEPTH",
            "N_REF_COUNT",
            "N_ALT_COUNT",
            "START_POSITION",
            "END_POSITION",
        ]
        actual_numerical_cols = []
        for col in numerical_cols:
            col_exists = process_functions.checkColExist(mutationDF, col)
            if col_exists:
                # Attempt to convert column to float
                try:
                    mutationDF[col] = mutationDF[col].astype(float)
                except ValueError:
                    pass
                if mutationDF[col].dtype not in [int, float]:
                    total_error.write(f"maf: {col} must be a numerical column.\n")
                else:
                    actual_numerical_cols.append(col)

        # CHECK: Must have TUMOR_SEQ_ALLELE2
        error, warn = _check_allele_col(mutationDF, "TUMOR_SEQ_ALLELE2")
        total_error.write(error)
        warning.write(warn)

        # CHECK: Mutation file would benefit from columns in optional_headers
        if (
            not all(
                [
                    process_functions.checkColExist(mutationDF, i)
                    for i in optional_headers
                ]
            )
            and not SP
        ):
            warning.write(
                "maf: Does not have the column headers that can give extra "
                "information to the processed maf: {}.\n".format(
                    ", ".join(
                        [
                            i
                            for i in optional_headers
                            if i not in mutationDF.columns.values
                        ]
                    )
                )
            )

        # CHECK: Must have REFERENCE_ALLELE
        error, warn = _check_allele_col(mutationDF, "REFERENCE_ALLELE")
        total_error.write(error)
        warning.write(warn)

        error, warn = validate._validate_chromosome(
            df=mutationDF, col="CHROMOSOME", fileformat="maf", allow_chr=False
        )
        total_error.write(error)
        warning.write(warn)
        # if process_functions.checkColExist(mutationDF, "CHROMOSOME"):
        #     # CHECK: Chromosome column can't have any values that start
        #     # with chr or have any WT values
        #     invalid_values = [
        #         str(i).startswith("chr") or str(i) == "WT"
        #         for i in mutationDF["CHROMOSOME"]
        #     ]
        #     if sum(invalid_values) > 0:
        #         total_error.write(
        #             "maf: CHROMOSOME column cannot have any values that "
        #             "start with 'chr' or any 'WT' values.\n"
        #         )

        error = _check_tsa1_tsa2(mutationDF)
        total_error.write(error)

        if process_functions.checkColExist(mutationDF, "TUMOR_SAMPLE_BARCODE"):
            error = process_functions.validate_genie_identifier(
                identifiers=mutationDF["TUMOR_SAMPLE_BARCODE"],
                center=self.center,
                filename="maf",
                col="TUMOR_SAMPLE_BARCODE",
            )
            total_error.write(error)

        # only check end position as start position is required col
        if (
            process_functions.checkColExist(mutationDF, "START_POSITION")
            and process_functions.checkColExist(mutationDF, "END_POSITION")
            and set(["START_POSITION", "END_POSITION"]) <= set(actual_numerical_cols)
        ):
            errors, warnings = validate.check_variant_start_and_end_positions(
                input_df=mutationDF,
                start_pos_col="START_POSITION",
                end_pos_col="END_POSITION",
                filename="maf",
            )
            total_error.write(errors)
            warning.write(warnings)

        for allele_col in self._allele_cols:
            if process_functions.checkColExist(mutationDF, allele_col):
                invalid_indices = validate.get_invalid_allele_rows(
                    mutationDF,
                    allele_col,
                    allowed_comb_alleles=self._allowed_comb_alleles,
                    allowed_ind_alleles=self._allowed_ind_alleles,
                    ignore_case=True,
                    allow_na=False,
                )
                errors, warnings = validate.get_allele_validation_message(
                    invalid_indices,
                    invalid_col=allele_col,
                    allowed_comb_alleles=self._allowed_comb_alleles,
                    allowed_ind_alleles=self._allowed_ind_alleles,
                    fileformat=self._fileType,
                )
                total_error.write(errors)
                warning.write(warnings)

        return total_error.getvalue(), warning.getvalue()

    def _cross_validate(self, mutationDF: pd.DataFrame) -> tuple:
        """This function cross-validates the mutation file to make sure it
        adheres to the mutation SOP.

        Args:
            mutationDF (pd.DataFrame): mutation dataframe

        Returns:
            Text with all the errors in the mutation file
        """
        errors = ""
        warnings = ""

        # This section can be removed once we remove the list of lists
        clinical_files = validate.parse_file_info_in_nested_list(
            nested_list=self.ancillary_files, search_str="data_clinical_supp"  # type: ignore[arg-type]
        )
        clinical_file_paths = clinical_files["file_info"]["path"]

        if clinical_files["files"]:
            try:
                clinical_sample_df = process_functions.get_clinical_dataframe(
                    filePathList=clinical_file_paths
                )
                has_file_read_error = False
            except Exception:
                has_file_read_error = True

            if not has_file_read_error:
                if process_functions.checkColExist(clinical_sample_df, "SAMPLE_ID"):
                    errors, warnings = validate.check_values_between_two_df(
                        df1=mutationDF,
                        df1_filename="MAF",
                        df1_id_to_check="TUMOR_SAMPLE_BARCODE",
                        df2=clinical_sample_df,
                        df2_filename="sample clinical",
                        df2_id_to_check="SAMPLE_ID",
                    )
        return errors, warnings

    def _get_dataframe(self, filePathList: List[str]) -> pd.DataFrame:
        """Get mutation dataframe

        1) Starts reading the first line in the file
        2) Skips lines that starts with #
        3) Reads in second line
        4) Checks that first line fields matches second line. Must do this because
        pandas.read_csv will allow for a file to have more column headers than content.
        E.g)  A,B,C,D,E
              1,2
              2,3

        5) We keep the 'NA', 'nan', and 'NaN' as strings in the data because
        these are valid allele values
        then convert the ones in the non-allele columns back to actual NAs

        NOTE: Because allele columns are case-insensitive in maf data, we must
        standardize the case of the columns when checking for the non-allele columns
        to convert the NA strings to NAs

        NOTE: This code allows empty dataframes to pass through
        without errors

        Args:
            filePathList (List[str]): list of filepath(s)

        Raises:
            ValueError: First line fields doesn't match second line fields in file

        Returns:
            pd.DataFrame: mutation data
        """
        with open(filePathList[0], "r") as maf_f:
            firstline = maf_f.readline()
            if firstline.startswith("#"):
                firstline = maf_f.readline()
            secondline = maf_f.readline()

        if len(firstline.split("\t")) != len(secondline.split("\t")):
            raise ValueError(
                "Number of fields in a line do not match the "
                "expected number of columns"
            )

        read_csv_params = {
            "filepath_or_buffer": filePathList[0],
            "sep": "\t",
            "comment": "#",
            "keep_default_na": False,
            "na_values": [
                "-1.#IND",
                "1.#QNAN",
                "1.#IND",
                "-1.#QNAN",
                "#N/A N/A",
                "#N/A",
                "N/A",
                "#NA",
                "NULL",
                "-NaN",
                "-nan",
                "",
            ],
            # This is to check if people write files
            # with R, quote=T
            "quoting": 3,
            # Retain completely blank lines so that
            # validator will cause the file to fail
            "skip_blank_lines": False,
        }

        mutationdf = transform._convert_df_with_mixed_dtypes(read_csv_params)

        mutationdf = transform._convert_values_to_na(
            input_df=mutationdf,
            values_to_replace=["NA", "nan", "NaN"],
            columns_to_convert=[
                col
                for col in mutationdf.columns
                if col.upper() not in self._allele_cols
            ],
        )
        return mutationdf
