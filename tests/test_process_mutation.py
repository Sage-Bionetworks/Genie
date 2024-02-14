"""Test process mutation functions"""

import os
from collections import namedtuple
from pandas.testing import assert_frame_equal
import shutil
import subprocess
import tempfile
from unittest.mock import patch, call

import pandas as pd
import pytest
import synapseclient

from genie import extract, load, process_mutation


def test_format_maf():
    maf_dict = {}
    maf_dict["Center"] = ["foo", "dsdf", "sdf"]
    maf_dict["Tumor_Sample_Barcode"] = ["GENIE-SAGE-1-3", "1-2", "3-2"]
    maf_dict["Sequence_Source"] = ["3", "e", "sd"]
    maf_dict["Sequencer"] = ["dsf", "sdf", "d"]
    maf_dict["Validation_Status"] = ["Unknown", "unknown", "f"]
    mafdf = pd.DataFrame(maf_dict)

    formatted_mafdf = process_mutation.format_maf(mafdf, center="SAGE")

    expected_maf_dict = {}
    expected_maf_dict["Center"] = ["SAGE", "SAGE", "SAGE"]
    expected_maf_dict["Tumor_Sample_Barcode"] = [
        "GENIE-SAGE-1-3",
        "GENIE-SAGE-1-2",
        "GENIE-SAGE-3-2",
    ]
    expected_maf_dict["Sequence_Source"] = [float("nan"), float("nan"), float("nan")]
    expected_maf_dict["Sequencer"] = [float("nan"), float("nan"), float("nan")]
    expected_maf_dict["Validation_Status"] = ["", "", "f"]
    expected_mafdf = pd.DataFrame(expected_maf_dict)
    assert expected_mafdf.equals(formatted_mafdf[expected_mafdf.columns])


class TestDtype:
    def setup_method(self):
        self.testdf = pd.DataFrame({"foo": [1], "bar": ["baz"]})
        self.column_types = {"foo": "int64", "bar": "object"}
        self.mutation_path = "/path/test.csv"
        self.input_dir = "/my/dir/here"
        self.final_maf_path = "/my/dir/here/test.csv"

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
        assert new_column_types == {"foo": "object", "bar": "object"}

    def test_move_maf_rename(self):
        """Test moving mafs when maf column headers need to be remapped"""
        testdf = pd.DataFrame({"CHROMOSOME": [1]})
        with patch.object(pd, "read_csv", return_value=testdf), patch.object(
            process_mutation, "determine_dtype", return_value=self.column_types
        ) as patch_determine, patch.object(
            process_mutation, "_convert_to_str_dtype", return_value=self.column_types
        ) as patch_convert, patch.object(
            testdf, "rename"
        ) as patch_rename, patch(
            "builtins.open"
        ) as patch_open:
            moved_maf = process_mutation.move_and_configure_maf(
                self.mutation_path, self.input_dir
            )
            patch_determine.assert_called_once_with(self.mutation_path)
            patch_convert.assert_called_once_with(
                self.column_types, process_mutation.KNOWN_STRING_COLS
            )
            patch_rename.assert_called_once_with(
                columns=process_mutation.MAF_COL_MAPPING
            )
            patch_open.assert_called_once_with(self.final_maf_path, "w")
            assert moved_maf == self.final_maf_path

    def test_move_mutation_vcf(self):
        """Test moving vcfs"""
        with patch.object(shutil, "copy") as patch_copy:
            process_mutation.move_mutation("/path/to/my.vcf", self.input_dir)
            patch_copy.assert_called_once_with("/path/to/my.vcf", self.input_dir)

    def test_move_mutation_maf(self):
        """Test moving maf files"""
        with patch.object(process_mutation, "move_and_configure_maf") as patch_move:
            process_mutation.move_mutation(self.mutation_path, self.input_dir)
            patch_move.assert_called_once_with(self.mutation_path, self.input_dir)


@pytest.fixture
def annotation_paths():
    Filepaths = namedtuple(
        "Filepaths",
        [
            "input_files_dir",
            "output_files_dir",
            "error_dir",
            "merged_maf_path",
            "narrow_maf_path",
            "full_maf_path",
            "full_error_report_path",
        ],
    )
    yield Filepaths(
        input_files_dir="input/dir",
        output_files_dir="input/dir",
        error_dir="input/dir/SAGE_error_reports",
        merged_maf_path="input/dir/data_mutations_extended_SAGE.txt",
        narrow_maf_path="input/SAGE/staging/data_mutations_extended_SAGE_MAF_narrow.txt",
        full_maf_path="input/SAGE/staging/data_mutations_extended_SAGE.txt",
        full_error_report_path="input/SAGE/staging/failed_annotations_report.txt",
    )


def test_process_mutation_workflow(syn, genie_config, annotation_paths):
    """Integration test to make sure workflow runs"""
    validfiles = pd.DataFrame(
        {"fileType": ["vcf", "maf"], "path": ["path/to/vcf", "path/to/maf"]}
    )
    genie_annotation_pkg = genie_config["genie_annotation_pkg"]
    syn_get_calls = [
        call(
            "syn22053204",
            ifcollision="overwrite.local",
            downloadLocation=genie_annotation_pkg,
            # version=1,  # TODO: This should pull from a config file in the future
        ),
        call(
            "syn22084320",
            ifcollision="overwrite.local",
            downloadLocation=genie_annotation_pkg,
            # version=13,  # TODO: This should pull from a config file in the future
        ),
    ]
    center = "SAGE"
    workdir = "working/dir/path"
    maf_table_id = "syn22493903"
    sample_error_report = pd.DataFrame({"col1": [1, 2, 3], "col2": [2, 3, 4]})
    with patch.object(
        process_mutation, "create_annotation_paths", return_value=annotation_paths
    ) as patch_annotation_paths, patch.object(syn, "get") as patch_synget, patch.object(
        process_mutation, "annotate_mutation"
    ) as patch_annotation, patch.object(
        process_mutation, "split_and_store_maf"
    ) as patch_split, patch.object(
        process_mutation,
        "concat_annotation_error_reports",
        return_value=sample_error_report,
    ) as patch_concat_error, patch.object(
        process_mutation, "check_annotation_error_reports"
    ) as patch_check_error, patch.object(
        process_mutation, "store_annotation_error_reports"
    ) as patch_store_error:
        maf = process_mutation.process_mutation_workflow(
            syn, center, validfiles, genie_config, workdir
        )
        patch_annotation_paths.assert_called_once_with(center=center, workdir=workdir)
        patch_synget.assert_has_calls(syn_get_calls)
        patch_annotation.assert_called_once_with(
            annotation_paths=annotation_paths,
            center=center,
            mutation_files=["path/to/vcf", "path/to/maf"],
            genie_annotation_pkg=genie_annotation_pkg,
        )
        patch_split.assert_called_once_with(
            syn=syn,
            center=center,
            maf_tableid=maf_table_id,
            annotation_paths=annotation_paths,
            flatfiles_synid="syn12279903",
        )
        patch_concat_error.assert_called_once_with(
            center=center,
            input_dir=annotation_paths.error_dir,
        )
        patch_check_error.assert_called_once_with(
            syn=syn,
            maf_table_synid=maf_table_id,
            full_error_report=sample_error_report,
            center=center,
        )
        patch_store_error.assert_called_once_with(
            full_error_report=sample_error_report,
            full_error_report_path=annotation_paths.full_error_report_path,
            syn=syn,
            errors_folder_synid="syn53239079",
        )
        assert maf == annotation_paths.merged_maf_path


def test_that_create_annotation_paths_returns_expected_paths(annotation_paths):
    center = "SAGE"
    input_dir = "input/dir"
    workdir = "test/dir"

    with patch.object(tempfile, "mkdtemp", return_value=input_dir) as patch_mktemp:
        test_paths = process_mutation.create_annotation_paths(
            center=center,
            workdir=workdir,
        )
        mktemp_calls = [call(dir=workdir)] * 2
        patch_mktemp.assert_has_calls(mktemp_calls)
        assert test_paths.input_files_dir == annotation_paths.input_files_dir
        assert test_paths.output_files_dir == annotation_paths.output_files_dir
        assert test_paths.error_dir == annotation_paths.error_dir
        assert test_paths.merged_maf_path == annotation_paths.merged_maf_path


class TestAnnotationErrorReports:
    @classmethod
    def setup_class(cls):
        cls.source_dir = "source_test_directory"
        os.makedirs(cls.source_dir, exist_ok=True)

        # Create sample files in the source directory
        with open(os.path.join(cls.source_dir, "file1.tsv"), "w") as f1:
            f1.write("col1\tcol2\tcol3\nval1\tval2\tval3\n")
        with open(os.path.join(cls.source_dir, "file2.tsv"), "w") as f2:
            f2.write("col1\tcol2\tcol3\nval4\tval5\tval6\n")

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.source_dir)

    @pytest.fixture
    def test_error_report(self):
        yield pd.DataFrame(
            {
                "col1": ["val1", "val4"],
                "col2": ["val2", "val5"],
                "col3": ["val3", "val6"],
                "Center": ["SAGE", "SAGE"],
            }
        )

    def test_check_annotation_error_reports_passes_if_match(
        self, syn, test_error_report
    ):
        maf_table_synid = "synZZZZ"
        with patch.object(
            extract, "get_syntabledf", return_value=test_error_report
        ) as patch_get_syntabledf:
            process_mutation.check_annotation_error_reports(
                syn=syn,
                maf_table_synid="synZZZZ",
                full_error_report=test_error_report,
                center="SAGE",
            )
            patch_get_syntabledf.assert_called_once_with(
                syn=syn,
                query_string=f"SELECT * FROM {maf_table_synid} "
                "WHERE Center = 'SAGE' AND "
                "Annotation_Status = 'FAILED'",
            )

    def test_check_annotation_error_reports_raises_warning_if_no_match(
        self, syn, test_error_report, caplog
    ):
        with patch.object(
            extract,
            "get_syntabledf",
            return_value=pd.DataFrame({"test": [1], "test2": [2]}),
        ):
            process_mutation.check_annotation_error_reports(
                syn=syn,
                maf_table_synid="synZZZZ",
                full_error_report=test_error_report,
                center="SAGE",
            )
            assert (
                "Genome nexus's failed annotations error report rows doesn't match"
                "maf table's failed annotations for SAGE" in caplog.text
            )
            # check that at least one is a warning
            assert any(record.levelname == "WARNING" for record in caplog.records)

    def test_concat_annotation_error_reports_returns_expected(self, test_error_report):
        compiled_report = process_mutation.concat_annotation_error_reports(
            input_dir="source_test_directory",
            center="SAGE",
        )
        assert_frame_equal(
            compiled_report.sort_values(by="col1").reset_index(drop=True),
            test_error_report.sort_values(by="col1").reset_index(drop=True),
        )

    def test_store_annotation_error_reports(self, syn, test_error_report):
        errors_folder_synid = "syn11111"
        full_error_report_path = "test.tsv"
        with patch.object(load, "store_file", return_value=None) as patch_store:
            process_mutation.store_annotation_error_reports(
                full_error_report=test_error_report,
                full_error_report_path=full_error_report_path,
                syn=syn,
                errors_folder_synid=errors_folder_synid,
            )
            patch_store.assert_called_once_with(
                syn=syn,
                filepath=full_error_report_path,
                parentid=errors_folder_synid,
            )


def test_annotate_mutation(annotation_paths):
    """Integration test, test that annotate mutation is called currectly"""
    center = "SAGE"
    mutation_files = ["path/to/vcf"]
    genie_annotation_pkg = "annotation/pkg/path"

    with patch.object(process_mutation, "move_mutation") as patch_move, patch.object(
        subprocess, "check_call"
    ) as patch_call:
        maf_path = process_mutation.annotate_mutation(
            annotation_paths=annotation_paths,
            center=center,
            mutation_files=mutation_files,
            genie_annotation_pkg=genie_annotation_pkg,
        )
        patch_move.assert_called_once_with(
            "path/to/vcf", annotation_paths.input_files_dir
        )
        patch_call.assert_called_once_with(
            [
                "bash",
                "annotation/pkg/path/annotation_suite_wrapper.sh",
                f"-i={annotation_paths.input_files_dir}",
                f"-o={annotation_paths.input_files_dir}",
                f"-e={annotation_paths.error_dir}",
                f"-m={annotation_paths.merged_maf_path}",
                f"-c={center}",
                "-s=WXS",
                f"-p={genie_annotation_pkg}",
            ]
        )


def test_append_or_createdf_append():
    """Test appending dataframe"""
    test_df = pd.DataFrame({"test": ["testme"]})
    with tempfile.NamedTemporaryFile() as temp_file, patch.object(
        test_df, "to_csv"
    ) as patch_tocsv:
        temp_file.write(b"Hello world!")
        temp_file.seek(0)
        process_mutation.append_or_createdf(test_df, temp_file.name)
        patch_tocsv.assert_called_once_with(
            temp_file.name, sep="\t", mode="a", index=False, header=None
        )


def test_append_or_createdf_create_none_exist_path():
    """Test creating dataframe when filepath passed in doesn't exist"""
    test_df = pd.DataFrame({"test": ["testme"]})
    with patch.object(test_df, "to_csv") as patch_tocsv:
        process_mutation.append_or_createdf(test_df, "DNE")
        patch_tocsv.assert_called_once_with("DNE", sep="\t", index=False)


def test_append_or_createdf_create_file_0size():
    """Test creating dataframe when file passed in is 0 in size"""
    test_df = pd.DataFrame({"test": ["testme"]})
    with tempfile.NamedTemporaryFile() as temp_file, patch.object(
        test_df, "to_csv"
    ) as patch_tocsv:
        temp_file.seek(0)
        process_mutation.append_or_createdf(test_df, temp_file.name)
        patch_tocsv.assert_called_once_with(temp_file.name, sep="\t", index=False)


def test_split_and_store_maf(syn, annotation_paths):
    """Integration test, check splitting and storing of maf functions are
    called"""
    # getTableColumns
    # pd.read_csv(
    # format_maf
    # append_or_createdf
    # store_narrow_maf
    # store_full_maf
    columns = [{"name": "inBED"}, {"name": "colA"}, {"name": "colB"}]
    exampledf = pd.DataFrame(
        {"colA": ["test", "foo"], "colB": ["bar", "baz"], "colC": [2, 3]}
    )
    annotated_maf_path = "maf/path"
    center = "SAGE"
    full_maf_path = "workdir/path/SAGE/staging/data_mutations_extended_SAGE.txt"
    narrow_maf_path = (
        "workdir/path/SAGE/staging/data_mutations_extended_SAGE_MAF_narrow.txt"
    )

    with patch.object(
        syn, "getTableColumns", return_value=columns
    ) as patch_getcols, patch.object(
        pd, "read_csv", return_value=[exampledf]
    ) as patch_readcsv, patch.object(
        process_mutation, "format_maf", return_value=exampledf
    ) as patch_format, patch.object(
        process_mutation, "append_or_createdf"
    ) as patch_append, patch.object(
        load, "store_table"
    ) as patch_store_table, patch.object(
        load, "store_file"
    ) as patch_store_file:
        process_mutation.split_and_store_maf(
            syn=syn,
            center=center,
            maf_tableid="sy12345",
            annotation_paths=annotation_paths,
            flatfiles_synid="syn2345",
        )
        patch_getcols.assert_called_once_with("sy12345")
        patch_readcsv.assert_called_once_with(
            annotation_paths.merged_maf_path, sep="\t", chunksize=100000, comment="#"
        )
        patch_format.assert_called_once_with(exampledf, center)

        assert patch_append.call_count == 2

        patch_store_table.assert_called_once_with(
            syn=syn, filepath=annotation_paths.narrow_maf_path, tableid="sy12345"
        )
        patch_store_file.assert_called_once_with(
            syn=syn, filepath=annotation_paths.full_maf_path, parentid="syn2345"
        )
