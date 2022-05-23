from unittest import mock

import pandas as pd

from genie_registry.structural_variant import StructuralVariant


class TestSv:

    def setup_method(self, syn):
        sv_cls = StructuralVariant(syn, "SAGE")
        self.sv_cls = sv_cls

    def test_filetype(self):
        assert self.sv_cls._fileType == "sv"

    def test_processing(self):
        expected_df = pd.DataFrame({
            "SAMPLE_ID": ['GENIE-SAGE-ID1-1', 'GENIE-SAGE-ID2-1', 'GENIE-SAGE-ID3-1',
                          'GENIE-SAGE-ID4-1', 'GENIE-SAGE-ID5-1'],
            "CENTER": ["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"]
        })
        sv_df = pd.DataFrame({
            "sample_id": ['GENIE-SAGE-ID1-1', 'GENIE-SAGE-ID2-1', 'GENIE-SAGE-ID3-1',
                          'GENIE-SAGE-ID4-1', 'GENIE-SAGE-ID5-1']
        })
        processed_df = self.sv_cls._process(sv_df)
        assert expected_df.equals(processed_df)

    def test_validation(self, syn):
        pass
