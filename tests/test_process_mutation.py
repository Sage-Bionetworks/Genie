"""Test process mutation functions"""
import shutil
from unittest.mock import patch

import pandas as pd

from genie import process_mutation


def test_format_maf():
    maf_dict = {}
    maf_dict['Center'] = ["foo", "dsdf", "sdf"]
    maf_dict['Tumor_Sample_Barcode'] = ["GENIE-SAGE-1-3", "1-2", "3-2"]
    maf_dict['Sequence_Source'] = ["3", "e", "sd"]
    maf_dict['Sequencer'] = ["dsf", "sdf", "d"]
    maf_dict['Validation_Status'] = ["Unknown", "unknown", "f"]
    mafdf = pd.DataFrame(maf_dict)

    formatted_mafdf = process_mutation.format_maf(mafdf, center="SAGE")

    expected_maf_dict = {}
    expected_maf_dict['Center'] = ["SAGE", "SAGE", "SAGE"]
    expected_maf_dict['Tumor_Sample_Barcode'] = [
        "GENIE-SAGE-1-3", "GENIE-SAGE-1-2", "GENIE-SAGE-3-2"
    ]
    expected_maf_dict['Sequence_Source'] = [float('nan'), float('nan'),
                                            float('nan')]
    expected_maf_dict['Sequencer'] = [float('nan'), float('nan'),
                                      float('nan')]
    expected_maf_dict['Validation_Status'] = ['', '', "f"]
    expected_mafdf = pd.DataFrame(expected_maf_dict)
    assert expected_mafdf.equals(formatted_mafdf[expected_mafdf.columns])


def test__rename_column_headers():
    """Tests the renaming of column headers"""
    testdf = pd.DataFrame({"foo": ["bar"],
                           "bar": ["baz"]})
    col_map = {"foo": "new_foo",
               "bar": "new_bar"}
    newdf = process_mutation._rename_column_headers(testdf, col_map)
    assert all(newdf.columns == ["new_foo", "new_bar"])


class TestDtype():
    def setup(self):
        self.testdf = pd.DataFrame({"foo": [1],
                                    "bar": ["baz"]})
        self.column_types = {"foo": 'int64',
                             "bar": 'object'}
        self.mutation_path = "/path/test.csv"
        self.input_dir = "/my/dir/here"

    def test__rename_column_headers(self):
        """Tests the renaming of column headers"""
        col_map = {"foo": "new_foo",
                   "bar": "new_bar"}
        newdf = process_mutation._rename_column_headers(self.testdf, col_map)
        assert all(newdf.columns == ["new_foo", "new_bar"])

    def test_determine_dtype(self):
        """Tests determining dtype"""
        with patch.object(pd, "read_csv", return_value=self.testdf):
            col_types = process_mutation.determine_dtype("test.csv")
            assert col_types == self.column_types

    def test__convert_to_str_dtype(self):
        """Tests converting dtypes to str dtypes"""
        new_column_types = process_mutation._convert_to_str_dtype(
            self.column_types, ["foo"]
        )
        assert new_column_types == {"foo": 'object',
                                    "bar": 'object'}

    def test_move_maf_copy(self):
        """Test moving mafs that don't need to rename columns"""
        with patch.object(pd, "read_csv",
                          return_value=self.testdf) as patch_read,\
             patch.object(shutil, "copy") as patch_copy:
            process_mutation.move_maf(self.mutation_path, self.input_dir)
            patch_read.assert_called_once_with(
                self.mutation_path, sep="\t", index_col=0, nrows=1,
                comment="#"
            )
            patch_copy.assert_called_once_with(self.mutation_path,
                                               self.input_dir)

    def test_move_maf_rename(self):
        """Test moving mafs when maf column headers need to be remapped"""
        testdf = pd.DataFrame({"CHROMOSOME": [1]})
        with patch.object(pd, "read_csv", return_value=testdf),\
            patch.object(process_mutation,
                         "determine_dtype",
                         return_value=self.column_types) as patch_determine,\
            patch.object(process_mutation,
                         "_convert_to_str_dtype",
                         return_value=self.column_types) as patch_convert,\
            patch.object(process_mutation,
                         "_rename_column_headers") as patch_rename,\
            patch.object(testdf, "to_csv"):
            process_mutation.move_maf(self.mutation_path, self.input_dir)
            patch_determine.assert_called_once_with(self.mutation_path)
            patch_convert.assert_called_once_with(
                self.column_types, process_mutation.KNOWN_STRING_COLS
            )
            patch_rename.assert_called_once_with(
                testdf, col_map=process_mutation.MAF_COL_MAPPING
            )
