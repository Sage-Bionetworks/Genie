"""Test process mutation functions"""
from distutils.command.build import build
import shutil
import subprocess
import tempfile
from unittest.mock import create_autospec, patch, call

import pandas as pd
import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError

from genie import process_mutation

SYN = create_autospec(synapseclient.Synapse)


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
    def setup(self):
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


def test_process_mutation_workflow(genie_config):
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
        ),
        call(
            "syn22084320",
            ifcollision="overwrite.local",
            downloadLocation=genie_annotation_pkg,
        ),
    ]
    center = "SAGE"
    workdir = "working/dir/path"
    maf_path = "path/to/maf"
    with patch.object(SYN, "get") as patch_synget, patch.object(
        process_mutation, "annotate_mutation", return_value=maf_path
    ) as patch_annotation, patch.object(
        process_mutation, "split_and_store_maf"
    ) as patch_split:

        maf = process_mutation.process_mutation_workflow(
            SYN, center, validfiles, genie_config, workdir
        )
        patch_synget.assert_has_calls(syn_get_calls)
        patch_annotation.assert_called_once_with(
            center=center,
            mutation_files=["path/to/vcf", "path/to/maf"],
            genie_annotation_pkg=genie_annotation_pkg,
            workdir=workdir,
        )
        patch_split.assert_called_once_with(
            syn=SYN,
            center=center,
            maf_tableid="syn22493903",
            annotated_maf_path=maf_path,
            flatfiles_synid="syn12279903",
            workdir=workdir,
        )
        assert maf == maf_path


def test_annotate_mutation():
    """Integration test, test that annotate mutation is called currectly"""
    center = "SAGE"
    mutation_files = ["path/to/vcf"]
    genie_annotation_pkg = "annotation/pkg/path"
    workdir = "working/dir/path"
    mktemp_calls = [call(dir=workdir)] * 2
    input_dir = "input/dir"
    with patch.object(
        tempfile, "mkdtemp", return_value=input_dir
    ) as patch_mktemp, patch.object(
        process_mutation, "move_mutation"
    ) as patch_move, patch.object(
        subprocess, "check_call"
    ) as patch_call:
        maf_path = process_mutation.annotate_mutation(
            center=center,
            mutation_files=mutation_files,
            genie_annotation_pkg=genie_annotation_pkg,
            workdir=workdir,
        )
        patch_mktemp.assert_has_calls(mktemp_calls)
        patch_move.assert_called_once_with("path/to/vcf", input_dir)
        patch_call.assert_called_once_with(
            [
                "bash",
                "annotation/pkg/path/annotation_suite_wrapper.sh",
                f"-i={input_dir}",
                f"-o={input_dir}",
                f"-m=input/dir/data_mutations_extended_{center}.txt",
                f"-c={center}",
                "-s=WXS",
                f"-p={genie_annotation_pkg}",
            ]
        )
        assert maf_path == f"input/dir/data_mutations_extended_{center}.txt"


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


def test_store_full_maf():
    """Test storing of full maf"""
    with patch.object(SYN, "store") as patch_store:
        process_mutation.store_full_maf(SYN, "full/path", "syn1234")
        patch_store.assert_called_once_with(
            synapseclient.File("full/path", parentId="syn1234")
        )


def test_store_narrow_maf():
    """Test storing of narrow maf"""
    with patch.object(SYN, "store") as patch_store:
        process_mutation.store_narrow_maf(SYN, "full/path", "syn1234")
        patch_store.assert_called_once()


def test_store_narrow_maf_test_error():
    """Test storing of narrow maf catches and passes error"""
    with patch.object(SYN, "store", side_effect=SynapseTimeoutError) as patch_store:
        process_mutation.store_narrow_maf(SYN, "full/path", "syn1234")
        patch_store.assert_called_once()


def test_split_and_store_maf():
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
        SYN, "getTableColumns", return_value=columns
    ) as patch_getcols, patch.object(
        pd, "read_csv", return_value=[exampledf]
    ) as patch_readcsv, patch.object(
        process_mutation, "format_maf", return_value=exampledf
    ) as patch_format, patch.object(
        process_mutation, "append_or_createdf"
    ) as patch_append, patch.object(
        process_mutation, "store_narrow_maf"
    ) as patch_narrow, patch.object(
        process_mutation, "store_full_maf"
    ) as patch_full:
        process_mutation.split_and_store_maf(
            syn=SYN,
            center=center,
            maf_tableid="sy12345",
            annotated_maf_path=annotated_maf_path,
            flatfiles_synid="syn2345",
            workdir="workdir/path",
        )
        patch_getcols.assert_called_once_with("sy12345")
        patch_readcsv.assert_called_once_with(
            annotated_maf_path, sep="\t", chunksize=100000, comment="#"
        )
        patch_format.assert_called_once_with(exampledf, center)

        assert patch_append.call_count == 2

        patch_narrow.assert_called_once_with(SYN, narrow_maf_path, "sy12345")
        patch_full.assert_called_once_with(SYN, full_maf_path, "syn2345")
