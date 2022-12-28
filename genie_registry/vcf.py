import logging
import os

import pandas as pd

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions, validate

logger = logging.getLogger(__name__)


def contains_whitespace(row):
    """Gets the total number of whitespaces from each column of a row"""
    return sum([" " in i for i in row if isinstance(i, str)])


class vcf(FileTypeFormat):
    _fileType = "vcf"

    _process_kwargs = []

    def _validateFilename(self, filePath):
        basename = os.path.basename(filePath[0])
        startswith_genie = basename.startswith("GENIE-{}-".format(self.center))
        endswith_vcf = basename.endswith(".vcf")
        assert startswith_genie and endswith_vcf

    def _get_dataframe(self, filePathList):
        headers = None
        filepath = filePathList[0]
        with open(filepath, "r") as vcffile:
            for row in vcffile:
                if row.startswith("#CHROM"):
                    headers = row.replace("\n", "").replace("\r", "").split("\t")
                    break
        if headers is not None:
            vcfdf = pd.read_csv(
                filepath, sep="\t", comment="#", header=None, names=headers
            )
        else:
            raise ValueError("Your vcf must start with the header #CHROM")
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
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        )
        total_error = ""
        warning = ""
        if not all(required_headers.isin(vcfdf.columns)):
            total_error += (
                "vcf: Must have these headers: "
                "CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.\n"
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
        # Vcf can only have max of 11 columns
        if len(vcfdf.columns) > 11:
            total_error += (
                "vcf: Should not have more than 11 columns. Only "
                "single sample or matched tumor normal vcf files are accepted.\n"
            )
        elif len(vcfdf.columns) > 8:
            # If there are greater than 8 columns, there must be the FORMAT column
            if "FORMAT" not in vcfdf.columns:
                total_error += "vcf: Must have FORMAT header if sample columns exist.\n"
            # If 11 columns, this is assumed to be a tumor normal vcf
            if len(vcfdf.columns) == 11:
                sample_id = vcfdf.columns[-2]
                normal_id = vcfdf.columns[-1]
                error = process_functions.validate_genie_identifier(
                    identifiers=pd.Series([sample_id]),
                    center=self.center,
                    filename="vcf",
                    col="tumor sample column",
                )
                total_error += error
                error = process_functions.validate_genie_identifier(
                    identifiers=pd.Series([normal_id]),
                    center=self.center,
                    filename="vcf",
                    col="normal sample column",
                )
                total_error += error
            else:
                # Everything else above 8 columns that isn't 11 columns
                # will be assumed to be a single sample vcf.
                # if TUMOR is not the sample column header, then validate
                # the sample column header.
                if "TUMOR" not in vcfdf.columns:
                    sample_id = vcfdf.columns[-1]
                    error = process_functions.validate_genie_identifier(
                        identifiers=pd.Series([sample_id]),
                        center=self.center,
                        filename="vcf",
                        col="tumor sample column",
                    )
                    if error:
                        error = error.replace("\n", "")
                        error += " if vcf represents a single sample and TUMOR is not the sample column header.\n"
                        total_error += error

        # Require that they report variants mapped to
        # either GRCh37 or hg19 without
        # the chr-prefix.
        error, warn = validate._validate_chromosome(
            df=vcfdf, col="#CHROM", fileformat="vcf"
        )
        total_error += error
        warning += warn

        # No white spaces
        white_space = vcfdf.apply(lambda x: contains_whitespace(x), axis=1)
        if sum(white_space) > 0:
            warning += "vcf: Should not have any white spaces in any of the columns.\n"

        # I can also recommend a `bcftools query` command that
        # will parse a VCF in a detailed way,
        # and output with warnings or errors if the format is not adhered too
        return total_error, warning
