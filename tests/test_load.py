from unittest.mock import Mock, patch

import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError

from genie import load, __version__


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


@pytest.mark.parametrize(
    "cols_subset, to_delete, subsetted_data",
    [
        (None, True, pd.DataFrame({"CENTER": ["test-center"], "data": [123]})),
        (None, False, pd.DataFrame({"CENTER": ["test-center"], "data": [123]})),
        (["CENTER"], True, pd.DataFrame({"CENTER": ["test-center"]})),
        (["CENTER", "col_extra"], False, pd.DataFrame({"CENTER": ["test-center"]})),
    ],
    ids=[
        "to_delete_is_true",
        "to_delete_is_false",
        "col_is_not_none",
        "col_has_column_not_in_db",
    ],
)
def test_that_update_table_has_expected_calls(
    syn, cols_subset, to_delete, subsetted_data
):
    test_table_synid = "synZZZZ"
    test_data = pd.DataFrame({"CENTER": ["test-center"], "data": [123]})
    test_new_data = pd.DataFrame({"CENTER": ["test-center"], "data": [123]})

    mock_database_ent = Mock()
    mock_database_ent.primaryKey = "PRIMARY_KEY"
    mock_database = Mock()

    with patch.object(syn, "get", return_value=mock_database_ent), patch.object(
        syn, "tableQuery", return_value=mock_database
    ) as patch_table_query, patch.object(
        mock_database, "asDataFrame", return_value=test_data
    ), patch.object(
        load, "_update_table"
    ) as patch__update_table:
        load.update_table(
            syn,
            databaseSynId=test_table_synid,
            newData=test_new_data,
            filterBy="test-center",
            filterByColumn="CENTER",
            col=cols_subset,
            toDelete=to_delete,
        )
        patch_table_query.assert_called_with(
            f"SELECT * FROM {test_table_synid} where CENTER ='test-center'"
        )

        # use this method to be able to compare dataframe arg values directly
        called_kwargs = patch__update_table.call_args.kwargs
        assert called_kwargs["syn"] == syn
        assert_frame_equal(called_kwargs["database"], subsetted_data)
        assert_frame_equal(called_kwargs["new_dataset"], test_new_data)
        assert called_kwargs["database_synid"] == test_table_synid
        assert called_kwargs["primary_key_cols"] == "PRIMARY_KEY"
        assert called_kwargs["to_delete"] == to_delete
