from unittest.mock import patch

import pandas as pd

from genie_registry.structural_variant import StructuralVariant
from genie import validate


class TestSv:
    def setup_method(self, syn):
        sv_cls = StructuralVariant(syn, "SAGE")
        self.sv_cls = sv_cls

    def test_filetype(self):
        assert self.sv_cls._fileType == "sv"

    def test_processing(self):
        expected_df = pd.DataFrame(
            {
                "SAMPLE_ID": [
                    "GENIE-SAGE-ID1-1",
                    "GENIE-SAGE-ID2-1",
                    "GENIE-SAGE-ID3-1",
                    "GENIE-SAGE-ID4-1",
                    "GENIE-SAGE-ID5-1",
                ],
                "CENTER": ["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
            }
        )
        sv_df = pd.DataFrame(
            {
                "sample_id": [
                    "GENIE-SAGE-ID1-1",
                    "GENIE-SAGE-ID2-1",
                    "GENIE-SAGE-ID3-1",
                    "GENIE-SAGE-ID4-1",
                    "GENIE-SAGE-ID5-1",
                ]
            }
        )
        processed_df = self.sv_cls._process(sv_df)
        assert expected_df.equals(processed_df)

    def test_validation_sample_error(self):
        sv_df = pd.DataFrame(
            {
                "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID1-1", "ID3-1"],
                "SV_STATUS": ["SOMATIC", "SOMATIC", "GERMLINE"],
                "SITE1_HUGO_SYMBOL": ["af", "af", "ff"],
                "SITE2_HUGO_SYMBOL": ["af", "af", "ff"],
            }
        )
        error, warning = self.sv_cls._validate(sv_df)
        assert error == (
            "Structural Variant: SAMPLE_ID must start with GENIE-SAGE\n"
            "Structural Variant: No duplicated rows allowed.\n"
        )
        assert warning == ""

    def test_validation_missing_required_cols(self):
        sv_df = pd.DataFrame(
            {
                "test": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID1-2", "ID3-1"],
                "foo": ["SOMATIC", "SOMATIC", "GERMLINE"],
            }
        )
        error, warning = self.sv_cls._validate(sv_df)
        assert error == (
            "Structural Variant: Must have SAMPLE_ID column.\n"
            "Structural Variant: Must have SV_STATUS column.\n"
            "Structural Variant: Either SITE1_HUGO_SYMBOL/SITE1_ENTREZ_GENE_ID or "
            "SITE2_HUGO_SYMBOL/SITE2_ENTREZ_GENE_ID is required.\n"
        )
        assert warning == ""

    def test_validation_integer_check(self):
        sv_df = pd.DataFrame(
            {
                "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID2-1"],
                "SV_STATUS": ["SOMATIC", "GERMLINE"],
                "SITE1_ENTREZ_GENE_ID": [1, "foo"],
                "SITE2_ENTREZ_GENE_ID": [1, "foo"],
                "SITE1_REGION_NUMBER": [1, "foo"],
                "SITE2_REGION_NUMBER": [1, "foo"],
                "SITE1_POSITION": [1, "foo"],
                "SITE2_POSITION": [1, "foo"],
                "TUMOR_SPLIT_READ_COUNT": [1, "foo"],
                "TUMOR_PAIRED_END_READ_COUNT": [1, "foo"],
                "SV_LENGTH": [1, "foo"],
                "NORMAL_READ_COUNT": [1, "foo"],
                "TUMOR_READ_COUNT": [1, "foo"],
                "NORMAL_VARIANT_COUNT": [1, "foo"],
                "TUMOR_VARIANT_COUNT": [1, "foo"],
                "NORMAL_PAIRED_END_READ_COUNT": [1, "foo"],
                "NORMAL_SPLIT_READ_COUNT": [1, "foo"],
            }
        )
        error, warning = self.sv_cls._validate(sv_df)
        assert error == (
            "Structural Variant: Only integers allowed in these column(s): "
            "SITE1_ENTREZ_GENE_ID, SITE2_ENTREZ_GENE_ID, "
            "SITE1_REGION_NUMBER, SITE2_REGION_NUMBER, SITE1_POSITION, SITE2_POSITION, "
            "TUMOR_SPLIT_READ_COUNT, TUMOR_PAIRED_END_READ_COUNT, SV_LENGTH, "
            "NORMAL_READ_COUNT, TUMOR_READ_COUNT, NORMAL_VARIANT_COUNT, "
            "TUMOR_VARIANT_COUNT, NORMAL_PAIRED_END_READ_COUNT, "
            "NORMAL_SPLIT_READ_COUNT.\n"
        )
        assert warning == ""

    def test_validation_no_errors(self):
        sv_df = pd.DataFrame(
            {
                "sample_id": [
                    "GENIE-SAGE-ID1-1",
                    "GENIE-SAGE-ID2-1",
                    "GENIE-SAGE-ID3-1",
                ],
                "SV_STATUS": ["SOMATIC", "GERMLINE", "GERMLINE"],
                "SITE1_ENTREZ_GENE_ID": [1, 2, 2],
                "SITE2_ENTREZ_GENE_ID": [1, 3, 3],
                "SITE1_REGION_NUMBER": [1, 2, 2],
                "NCBI_BUILD": ["GRCh38", float("nan"), "GRCh37"],
                "BREAKPOINT_TYPE": ["PRECISE", "IMPRECISE", "IMPRECISE"],
                "CONNECTION_TYPE": ["3to5", "5to5", "5to5"],
                "DNA_SUPPORT": ["Yes", "No", "Unknown"],
                "RNA_Support": ["Yes", "No", "Unknown"],
                "SITE1_CHROMOSOME": [1, 22, float("nan")],
                "SITE2_CHROMOSOME": ["X", "2", float("nan")],
                "SITE1_REGION": ["IGR", "Upstream", "5_Prime_UTR Intron"],
                "SITE2_REGION": ["3-UTR", "3_Prime_UTR Intron", "Exon"],
            }
        )
        error, warning = self.sv_cls._validate(sv_df)
        assert error == ""
        assert warning == ""

    def test_validation__validate_chromosome_is_called(self):
        """Tests that _validate_chromosome is called twice
        for the two chromosome columns that sv files may have"""
        sv_df = pd.DataFrame(
            {
                "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID2-1"],
                "SITE1_CHROMOSOME": [1, 22],
            }
        )
        with patch.object(
            validate, "_validate_chromosome", return_value=("", "")
        ) as validation__validate_chromosome_mock:
            self.sv_cls._validate(sv_df)
            assert validation__validate_chromosome_mock.call_count == 2, (
                "_validate_chromosome should be called twice for sv file"
                "since it has two potential chromosome columns to check"
            )
