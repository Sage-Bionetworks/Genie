from io import StringIO
import logging
import os

from genie.example_filetype_format import FileTypeFormat
from genie import load, process_functions, validate

logger = logging.getLogger(__name__)


class StructuralVariant(FileTypeFormat):
    _fileType = "sv"

    # _validation_kwargs = ["nosymbol_check", "project_id"]

    # VALIDATE FILENAME
    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "data_sv.txt"

    def _process(self, sv_df):
        sv_df.columns = [col.upper() for col in sv_df.columns]
        # Add center column
        center = [sample_id.split("-")[1] for sample_id in sv_df["SAMPLE_ID"]]
        sv_df["CENTER"] = center
        return sv_df

    def process_steps(self, sv_df, newPath, databaseSynId):
        sv_df = self._process(sv_df)
        # TODO: test the col parameter
        load.update_table(
            syn=self.syn,
            databaseSynId=databaseSynId,
            newData=sv_df,
            filterBy=self.center,
            toDelete=True,
            col=sv_df.columns.to_list(),
        )
        sv_df.to_csv(newPath, sep="\t", index=False)
        return newPath

    # def _validate(self, sv_df, nosymbol_check, project_id):
    def _validate(self, sv_df):
        total_error = StringIO()
        total_warning = StringIO()
        sv_df.columns = [col.upper() for col in sv_df.columns]

        have_sample_col = process_functions.checkColExist(sv_df, "SAMPLE_ID")
        if not have_sample_col:
            total_error.write("Structural Variant: Must have SAMPLE_ID column.\n")
        else:
            # if sv_df["SAMPLE_ID"].duplicated().any():
            #     total_error.write(
            #         "Structural Variant: No duplicated SAMPLE_ID allowed.\n"
            #     )
            # TODO: switch to validate_genie_identifier function
            # After GH-444 is merged
            errors = process_functions.validate_genie_identifier(
                identifiers=sv_df["SAMPLE_ID"],
                center=self.center,
                filename="Structural Variant",
                col="SAMPLE_ID",
            )
            total_error.write(errors)

        if sv_df.duplicated().any():
            total_error.write("Structural Variant: No duplicated rows allowed.\n")

        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="SV_STATUS",
            possible_values=["SOMATIC", "GERMLINE"],
            filename="Structural Variant",
            required=True,
        )
        total_warning.write(warn)
        total_error.write(error)

        have_hugo_1 = process_functions.checkColExist(sv_df, "SITE1_HUGO_SYMBOL")
        have_hugo_2 = process_functions.checkColExist(sv_df, "SITE2_HUGO_SYMBOL")
        have_entrez_1 = process_functions.checkColExist(sv_df, "SITE1_ENTREZ_GENE_ID")
        have_entrez_2 = process_functions.checkColExist(sv_df, "SITE2_ENTREZ_GENE_ID")

        if not ((have_hugo_1 or have_entrez_1) and (have_hugo_2 or have_entrez_2)):
            total_error.write(
                "Structural Variant: Either SITE1_HUGO_SYMBOL/SITE1_ENTREZ_GENE_ID "
                "or SITE2_HUGO_SYMBOL/SITE2_ENTREZ_GENE_ID is required.\n"
            )

        # optional_columns = [
        #     "SITE1_REGION_NUMBER",
        #     "SITE2_REGION_NUMBER",
        #     "SITE1_REGION",
        #     "SITE2_REGION",
        #     "SITE1_CHROMOSOME",
        #     "SITE2_CHROMOSOME",
        #     "SITE1_CONTIG",
        #     "SITE2_CONTIG",
        #     "SITE1_POSITION",
        #     "SITE2_POSITION",
        #     "SITE1_DESCRIPTION",
        #     "SITE2_DESCRIPTION",
        #     "SITE2_EFFECT_ON_FRAME",
        #     "NCBI_BUILD",
        #     "CLASS",
        #     "TUMOR_SPLIT_READ_COUNT",
        #     "TUMOR_PAIRED_END_READ_COUNT",
        #     "EVENT_INFO",
        #     "BREAKPOINT_TYPE",
        #     "CONNECTION_TYPE",
        #     "ANNOTATION",
        #     "DNA_SUPPORT",
        #     "RNA_SUPPORT",
        #     "SV_LENGTH",
        #     "NORMAL_READ_COUNT",
        #     "TUMOR_READ_COUNT",
        #     "NORMAL_VARIANT_COUNT",
        #     "TUMOR_VARIANT_COUNT",
        #     "NORMAL_PAIRED_END_READ_COUNT",
        #     "NORMAL_SPLIT_READ_COUNT",
        #     "COMMENTS",
        # ]
        # Check for columns that should be integar columsn
        int_cols = [
            "SITE1_ENTREZ_GENE_ID",
            "SITE2_ENTREZ_GENE_ID",
            "SITE1_REGION_NUMBER",
            "SITE2_REGION_NUMBER",
            "SITE1_POSITION",
            "SITE2_POSITION",
            "TUMOR_SPLIT_READ_COUNT",
            "TUMOR_PAIRED_END_READ_COUNT",
            "SV_LENGTH",
            "NORMAL_READ_COUNT",
            "TUMOR_READ_COUNT",
            "NORMAL_VARIANT_COUNT",
            "TUMOR_VARIANT_COUNT",
            "NORMAL_PAIRED_END_READ_COUNT",
            "NORMAL_SPLIT_READ_COUNT",
        ]
        # Get all columns that are non integers.
        # Allow for NA's
        non_ints = [
            col
            for col in int_cols
            if sv_df.get(col) is not None
            and not sv_df[col].dropna().apply(process_functions.checkInt).all()
        ]
        if len(non_ints) > 0:
            total_error.write(
                "Structural Variant: Only integers allowed in these "
                "column(s): {}.\n".format(", ".join(non_ints))
            )

        region_allow_vals = [
            "5_prime_UTR",
            "3_prime_UTR",
            "Promoter",
            "Exon",
            "Intron",
            "5'UTR",
            "3'UTR",
            "5-UTR",
            "3-UTR",
            "5_Prime_UTR Intron",
            "3_Prime_UTR Intron",
            "Downstream",
            "Upstream",
            "Intergenic",
            "IGR",
        ]
        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="SITE1_REGION",
            possible_values=region_allow_vals,
            filename="Structural Variant",
            na_allowed=True,
            required=False,
        )
        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="SITE2_REGION",
            possible_values=region_allow_vals,
            filename="Structural Variant",
            na_allowed=True,
            required=False,
        )
        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="NCBI_BUILD",
            possible_values=["GRCh37", "GRCh38"],
            filename="Structural Variant",
            na_allowed=True,
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)

        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="BREAKPOINT_TYPE",
            possible_values=["PRECISE", "IMPRECISE"],
            filename="Structural Variant",
            na_allowed=True,
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)

        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="CONNECTION_TYPE",
            possible_values=["3to5", "5to3", "5to5", "3to3"],
            filename="Structural Variant",
            na_allowed=True,
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)

        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="DNA_SUPPORT",
            possible_values=["Yes", "No", "Unknown"],
            filename="Structural Variant",
            na_allowed=True,
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)
        warn, error = process_functions.check_col_and_values(
            df=sv_df,
            col="RNA_SUPPORT",
            possible_values=["Yes", "No", "Unknown"],
            filename="Structural Variant",
            na_allowed=True,
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)
        # check for chromosome columns and don't allow 'chr' for now
        # since in the database there’s nothing with CHR
        chrom_cols = ["SITE1_CHROMOSOME", "SITE2_CHROMOSOME"]
        for chrom_col in chrom_cols:
            error, warn = validate._validate_chromosome(
                df=sv_df,
                col=chrom_col,
                fileformat="Structural Variant",
                allow_chr=False,
                allow_na=True,
            )
            total_error.write(error)

        return total_error.getvalue(), total_warning.getvalue()
