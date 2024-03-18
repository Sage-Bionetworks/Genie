from unittest.mock import patch, mock_open

import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError

from genie import load, __version__
import pandas as pd
import tempfile
import os


def test_store_file(syn):
    """Test storing of file"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234")
        patch_store.assert_called_once_with(
            synapseclient.File("full/path", parentId="syn1234", versionComment=None),
            used=None,
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_file_version_comment(syn):
    """Test storing of file"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234", version_comment="test")
        patch_store.assert_called_once_with(
            synapseclient.File("full/path", parentId="syn1234", versionComment="test"),
            used=None,
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_file_version_provenance(syn):
    """Test storing of file"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234", used="temp")
        patch_store.assert_called_once_with(
            synapseclient.File("full/path", parentId="syn1234", versionComment=None),
            used="temp",
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_file_annotations(syn):
    """Test storing of file with annotations"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234", annotations={"test": "test"})
        patch_store.assert_called_once_with(
            synapseclient.File(
                "full/path",
                parentId="syn1234",
                versionComment=None,
                annotations={"test": "test"},
            ),
            used=None,
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_table(syn):
    """Test storing of table"""
    with patch.object(syn, "store") as patch_store:
        load.store_table(syn, "full/path", "syn1234")
        patch_store.assert_called_once()


def test_store_table_error(syn):
    """Test storing of table catches and passes error"""
    with patch.object(syn, "store", side_effect=SynapseTimeoutError) as patch_store:
        load.store_table(syn, "full/path", "syn1234")
        patch_store.assert_called_once()


def test__update_table_non_empty_dataframe(syn):
    """Test _update_table function with both new_dataset and database as non_empty dataframe"""
    database_synid = "syn123"
    primary_key_cols = ["test", "foo"]
    to_delete = (False,)
    new_dataset = pd.DataFrame(
        {
            "test": ["test1", "test2", "test3", "test4"],
            "foo": [1, 2, 3, 4],
            "baz": [float("nan"), float("nan"), float("nan"), 3.2],
        }
    )
    database = pd.DataFrame(
        {
            "test": ["test1", "test2", "test3"],
            "foo": [1, 2, 3],
            "baz": [float("nan"), float("nan"), float("nan")],
        }
    )
    expected_results = ["ROW_ID,ROW_VERSION,test,foo,baz\n", ",,test4,4,3.2\n"]
    with patch("os.unlink") as mock_unlink, patch(
        "tempfile.NamedTemporaryFile"
    ) as mock_tempfile:
        with patch("builtins.open", mock_open()) as mock_file_open:
            # set the tempfile name
            mock_tempfile.return_value.name = "test.csv"
            load._update_table(
                syn, database, new_dataset, database_synid, primary_key_cols, to_delete
            )
            mock_file_open.assert_called_once_with("test.csv", "w")
            mock_file_handle = mock_file_open()
            write_calls = mock_file_handle.write.call_args_list
            results = [call_args[0][0] for call_args in write_calls]
            assert results == expected_results
            mock_unlink.assert_called_once_with("test.csv")


def test__update_table_empty_dataframe(syn):
    """Test _update_table function with empty new_dataset"""
    database_synid = "syn123"
    primary_key_cols = ["test", "foo"]
    to_delete = False
    new_dataset = pd.DataFrame(columns=["test", "foo", "baz"])
    database = pd.DataFrame(
        {
            "test": ["test1", "test2", "test3"],
            "foo": [1, 2, 3],
            "baz": [float("nan"), float("nan"), float("nan")],
        }
    )
    expected_results = ["ROW_ID,ROW_VERSION,test,foo,baz\n"]
    with patch("os.unlink") as mock_unlink, patch(
        "tempfile.NamedTemporaryFile"
    ) as mock_tempfile:
        with patch("builtins.open", mock_open()) as mock_file_open:
            # set the tempfile name
            mock_tempfile.return_value.name = "test.csv"
            load._update_table(
                syn, database, new_dataset, database_synid, primary_key_cols, to_delete
            )
            mock_file_open.assert_called_once_with("test.csv", "w")
            mock_file_handle = mock_file_open()
            write_calls = mock_file_handle.write.call_args_list
            results = [call_args[0][0] for call_args in write_calls]
            assert results == expected_results
            mock_unlink.assert_called_once_with("test.csv")


def test__update_table_empty_dataframes(syn):
    """Test _update_table function with empty new_dataset and database"""
    database_synid = "syn123"
    primary_key_cols = ["test", "foo"]
    to_delete = False
    new_dataset = pd.DataFrame(columns=["test", "foo", "baz"])
    database = pd.DataFrame(columns=["test", "foo", "baz"])
    expected_results = ["ROW_ID,ROW_VERSION,test,foo,baz\n"]
    with patch("os.unlink") as mock_unlink, patch(
        "tempfile.NamedTemporaryFile"
    ) as mock_tempfile:
        with patch("builtins.open", mock_open()) as mock_file_open:
            # set the tempfile name
            mock_tempfile.return_value.name = "test.csv"
            load._update_table(
                syn, database, new_dataset, database_synid, primary_key_cols, to_delete
            )
            mock_file_open.assert_called_once_with("test.csv", "w")
            mock_file_handle = mock_file_open()
            write_calls = mock_file_handle.write.call_args_list
            results = [call_args[0][0] for call_args in write_calls]
            assert results == expected_results
            mock_unlink.assert_called_once_with("test.csv")
