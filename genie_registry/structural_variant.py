from io import StringIO
import logging
import os

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions

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
        process_functions.updateData(
            self.syn,
            databaseSynId,
            sv_df,
            self.center,
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
            if sv_df["SAMPLE_ID"].duplicated().any():
                total_error.write(
                    "Structural Variant: No duplicated SAMPLE_ID allowed.\n"
                )
            # TODO: switch to validate_genie_identifier function
            # After GH-444 is merged
            if not all(sv_df["SAMPLE_ID"].str.startswith(f"GENIE-{self.center}")):
                total_error.write(
                    "Structural Variant: SAMPLE_ID must start with "
                    f"GENIE-{self.center}\n"
                )

        warn, error = process_functions.check_col_and_values(
            sv_df,
            "SV_STATUS",
            ["SOMATIC", "GERMLINE"],
            "Structural Variant",
            required=True,
        )
        total_warning.write(warn)
        total_error.write(error)

        # optional_columns = [
        #     "SITE1_HUGO_SYMBOL",
        #     "SITE2_HUGO_SYMBOL",
        #     "SITE1_ENSEMBL_TRANSCRIPT_ID",
        #     "SITE2_ENSEMBL_TRANSCRIPT_ID",
        #     "SITE1_ENTREZ_GENE_ID",
        #     "SITE2_ENTREZ_GENE_ID",
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
        non_ints = [
            col
            for col in int_cols
            if sv_df.get(col) is not None and sv_df[col].dtype != int
        ]
        if len(non_ints) > 0:
            total_error.write(
                "Structural Variant: Only integers allowed in these "
                "column(s): {}.\n".format(", ".join(non_ints))
            )

        # region_allow_vals = [
        #     "5_PRIME_UTR", "3_PRIME_UTR", "PROMOTER", "EXON", "INTRON"
        # ]
        # warn, error = process_functions.check_col_and_values(
        #     sv_df, "SITE1_REGION", region_allow_vals,
        #     "Structural Variant", required=False
        # )
        # warn, error = process_functions.check_col_and_values(
        #     sv_df, "SITE2_REGION", region_allow_vals,
        #     "Structural Variant", required=False
        # )
        warn, error = process_functions.check_col_and_values(
            sv_df,
            "NCBI_BUILD",
            ["GRCh37", "GRCh38"],
            "Structural Variant",
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)

        warn, error = process_functions.check_col_and_values(
            sv_df,
            "BREAKPOINT_TYPE",
            ["PRECISE", "IMPRECISE"],
            "Structural Variant",
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)

        warn, error = process_functions.check_col_and_values(
            sv_df,
            "CONNECTION_TYPE",
            ["3to5", "5to3", "5to5", "3to3"],
            "Structural Variant",
            required=False,
        )
        # total_warning.write(warn)
        total_error.write(error)

        warn, error = process_functions.check_col_and_values(
            sv_df, "DNA_SUPPORT", ["Yes", "No"], "Structural Variant", required=False
        )
        # total_warning.write(warn)
        total_error.write(error)

        warn, error = process_functions.check_col_and_values(
            sv_df, "RNA_SUPPORT", ["Yes", "No"], "Structural Variant", required=False
        )
        # total_warning.write(warn)
        total_error.write(error)

        return total_error.getvalue(), total_warning.getvalue()
