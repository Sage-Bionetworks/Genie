import logging
import os
from typing import List

import pandas as pd

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions, transform, validate

logger = logging.getLogger(__name__)


def contains_whitespace(row):
    """Gets the total number of whitespaces from each column of a row"""
    return sum([" " in i for i in row if isinstance(i, str)])


class vcf(FileTypeFormat):
    _fileType = "vcf"

    _process_kwargs = []
    _allele_cols = ["REF"]
    _allowed_comb_alleles = ["A", "T", "C", "G", "N"]
    _allowed_ind_alleles = []

    def _validateFilename(self, filePath):
        basename = os.path.basename(filePath[0])
        startswith_genie = basename.startswith("GENIE-{}-".format(self.center))
        endswith_vcf = basename.endswith(".vcf")
        assert startswith_genie and endswith_vcf

    def _get_dataframe(self, filePathList: List[str]) -> pd.DataFrame:
        """Get mutation dataframe

        1) Looks for the line in the file starting with #CHROM, that will be
        the header line (columns).

        2) When reading in the data, we keep the 'NA', 'nan', and 'NaN'
        as strings in the data because these are valid allele values
        then convert the ones in the non-allele columns back to actual NAs

        Args:
            filePathList (List[str]): list of filepath(s)

        Raises:
            ValueError: when line with #CHROM doesn't exist in file

        Returns:
            pd.DataFrame: mutation data
        """
        headers = None
        filepath = filePathList[0]
        with open(filepath, "r") as vcffile:
            for row in vcffile:
                if row.startswith("#CHROM"):
                    headers = row.replace("\n", "").replace("\r", "").split("\t")
                    break
        if headers is not None:
            vcfdf = pd.read_csv(
                filepath,
                sep="\t",
                comment="#",
                header=None,
                names=headers,
                keep_default_na=False,
                na_values=[
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
            )
        else:
            raise ValueError("Your vcf must start with the header #CHROM")

        vcfdf = transform._convert_values_to_na(
            input_df=vcfdf,
            values_to_replace=["NA", "nan", "NaN"],
            columns_to_convert=[
                col for col in vcfdf.columns if col not in self._allele_cols
            ],
        )
        return vcfdf

    def process_steps(self, df):
        """The processing of vcf files is specific to GENIE, so
        not included in this function"""
        logger.info(
            "Please run with `--process mutation` parameter "
            "if you want to reannotate the mutation files"
        )
        return None

    def _validate(self, vcfdf):
        """
        Validates the content of a vcf file

        Args:
            vcfdf: pandas dataframe containing vcf content

        Returns:
            total_error - error messages
            warning - warning messages
        """
        required_headers = pd.Series(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        )
        total_error = ""
        warning = ""
        if not all(required_headers.isin(vcfdf.columns)):
            total_error += (
                "vcf: Must have these headers: "
                "CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT.\n"
            )
        else:
            # No duplicated values
            primary_cols = ["#CHROM", "POS", "REF", "ALT"]
            if vcfdf.duplicated(primary_cols).any():
                total_error += "vcf: Must not have duplicate variants.\n"

            if vcfdf[["#CHROM", "POS"]].isnull().values.any():
                total_error += (
                    "vcf: May contain rows that are "
                    "space delimited instead of tab delimited.\n"
                )
            if vcfdf["FORMAT"].isnull().values.any():
                total_error += "vcf: Must not have missing values in FORMAT column.\n"
        total_error += self.validate_tumor_and_normal_sample_columns_exist(
            input_df=vcfdf
        )
        # Require that they report variants mapped to
        # either GRCh37 or hg19 without
        # the chr-prefix.
        error, warn = validate._validate_chromosome(
            df=vcfdf, col="#CHROM", fileformat="vcf"
        )
        total_error += error
        warning += warn

        for allele_col in self._allele_cols:
            if process_functions.checkColExist(vcfdf, allele_col):
                invalid_indices = validate.get_invalid_allele_rows(
                    vcfdf,
                    input_col=allele_col,
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
                total_error += errors
                warning += warnings

        # No white spaces
        white_space = vcfdf.apply(lambda x: contains_whitespace(x), axis=1)
        if sum(white_space) > 0:
            warning += "vcf: Should not have any white spaces in any of the columns.\n"

        # I can also recommend a `bcftools query` command that
        # will parse a VCF in a detailed way,
        # and output with warnings or errors if the format is not adhered too
        return total_error, warning

    def validate_tumor_and_normal_sample_columns_exist(
        self, input_df: pd.DataFrame
    ) -> str:
        """Validates that the expected tumor sample column and optional normal
            sample columns are present in the VCF depending on how many
            columns you have present in the VCF and they have no missing values

            Rules:
                - VCFs can only have a max of 11 columns including the 9 required columns
                - For 11 columns VCFs, it is assumed this is a matched tumor normal vcf file
                    which means there should be a tumor sample and normal sample
                    column present
                - For 10 column VCFs, it is assumed this is a single sample vcf file
                    which means there should be a tumor sample column present
                - Anything lower than 10 columns is INVALID because you must have at
                least a tumor sample column on top of the 9 required VCF columns

                - If tumor sample and normal sample columns are present, they must not have
                any missing values.

            Examples:

            VCF with Matched Tumor Normal columns:
            | GENIE-GOLD-1-1-tumor | GENIE-GOLD-1-1-normal |
            | -------------------- | --------------------- |
            |                      |                       |

            VCF with Single Tumor VCF column:
            | TUMOR |
            | ----- |
            |       |

        Args:
            input_df (pd.DataFrame): input vcf data to be validated

        Returns:
            str: error message
        """
        error = ""
        sample_id = None
        normal_id = None
        # vcf can only have max of 11 columns
        if len(input_df.columns) > 11:
            error = (
                "vcf: Should not have more than 11 columns. Only "
                "single sample or matched tumor normal vcf files are accepted.\n"
            )
        # If 11 columns, this is assumed to be a tumor normal vcf
        # so it must have tumor sample and normal sample columns
        elif len(input_df.columns) == 11:
            sample_id = input_df.columns[-2]
            normal_id = input_df.columns[-1]
            error = process_functions.validate_genie_identifier(
                identifiers=pd.Series([sample_id]),
                center=self.center,
                filename="vcf",
                col="tumor sample column",
            )
            error += process_functions.validate_genie_identifier(
                identifiers=pd.Series([normal_id]),
                center=self.center,
                filename="vcf",
                col="normal sample column",
            )
        elif len(input_df.columns) == 10:
            # If 10 columns, it will be assumed to be a single sample vcf.
            # if TUMOR is not the sample column header, then validate
            # the sample column header.
            if "TUMOR" not in input_df.columns:
                sample_id = input_df.columns[-1]
                error = process_functions.validate_genie_identifier(
                    identifiers=pd.Series([sample_id]),
                    center=self.center,
                    filename="vcf",
                    col="tumor sample column",
                )
                if error:
                    error = error.replace("\n", "")
                    error += " if vcf represents a single sample and TUMOR is not the sample column header.\n"
            else:
                sample_id = "TUMOR"
        else:
            # Must have a column called TUMOR or sample column in the header
            error = (
                "vcf: Must have at least 10 columns. "
                "If the vcf represents a single sample, then it's missing a tumor sample column. "
                "If the vcf represents a matched tumor normal, then it's missing both normal sample and tumor sample columns.\n"
            )

        # validate the values in the tumor and/or normal sample columns if present
        if sample_id:
            if input_df[sample_id].isnull().values.any():
                error += f"vcf: Must not have missing values in {sample_id} column.\n"
        if normal_id:
            if input_df[normal_id].isnull().values.any():
                error += f"vcf: Must not have missing values in {normal_id} column.\n"
        return error
