from unittest.mock import Mock, patch, call

import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError
from synapseclient.models import Table
from genie import load, __version__, process_functions


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


def test__get_col_order():
    database = pd.DataFrame(columns=["a", "b", "c"])
    orig_database_cols = database.columns
    expected = ["ROW_ID", "ROW_VERSION", "a", "b", "c"]
    results = load._get_col_order(orig_database_cols)
    assert results == expected


@pytest.mark.parametrize(
    "ori_dataset,new_dataset",
    [
        (
            pd.DataFrame(columns=["a", "b", "c"]),
            pd.DataFrame(columns=["b", "c", "a"]),
        ),
    ],
)
def test__reorder_new_dataset(ori_dataset, new_dataset):
    """Test if new dataset is re-ordered as the column order of original dataset.
    No need to check different number of columns since we are using commmon columns as
    specified in load.update_table()

    """
    orig_database_cols = ori_dataset.columns
    reorder_new_dataset = load._reorder_new_dataset(orig_database_cols, new_dataset)
    pd.testing.assert_index_equal(reorder_new_dataset.columns, orig_database_cols)


@pytest.mark.parametrize(
    "dataset,expected",
    [
        (
            pd.DataFrame(
                {
                    "test": ["test1", "test2", "test3", "test4"],
                    "foo": [1.0, 2.0, float("nan"), float("nan")],
                    "baz": ["b", float("nan"), "c", float("nan")],
                }
            ),
            pd.DataFrame(
                {
                    "test": ["test1", "test2", "test3", "test4"],
                    "foo": [1.0, 2.0, float("nan"), float("nan")],
                    "baz": ["b", float("nan"), "c", float("nan")],
                    "UNIQUE_KEY": ["test1 1.0 b", "test2 2.0 ", "test3  c", "test4  "],
                }
            ),
        ),
        (
            pd.DataFrame(columns=["test", "foo", "baz"]),
            pd.DataFrame(columns=["test", "foo", "baz", "UNIQUE_KEY"]),
        ),
    ],
    ids=["non_empty dataframe", "empty dataframe"],
)
def test__generate_primary_key(dataset, expected):
    results = load._generate_primary_key(
        dataset=dataset,
        primary_key_cols=["test", "foo", "baz"],
        primary_key="UNIQUE_KEY",
    )
    pd.testing.assert_frame_equal(results, expected, check_dtype=False)


@pytest.mark.parametrize(
    "to_delete,expected_to_delete_rows",
    [
        (True, pd.DataFrame({0: ["2", "3"], 1: ["3", "5"]})),
        (False, pd.DataFrame()),
    ],
    ids=["to_delted is True", "to_delted is False"],
)
def test_check_database_changes(to_delete, expected_to_delete_rows):
    with (
        patch.object(
            process_functions,
            "_append_rows",
            return_value=pd.DataFrame({"test": ["test4"], "foo": [4], "baz": [3.2]}),
        ) as patch_append_rows,
        patch.object(
            process_functions,
            "_update_rows",
            return_value=pd.DataFrame(
                {
                    "test": ["test", "test2"],
                    "foo": [1, 3],
                    "baz": ["", 5],
                    "ROW_ID": ["1", "2"],
                    "ROW_VERSION": ["3", "3"],
                }
            ),
        ) as patch_update,
        patch.object(
            process_functions,
            "_delete_rows",
            return_value=pd.DataFrame({0: ["2", "3"], 1: ["3", "5"]}),
        ) as patch_delete,
        patch.object(load, "_generate_primary_key", return_value="table"),
        patch.object(
            load,
            "_get_col_order",
            return_value=["ROW_ID", "ROW_VERSION", "test", "foo", "baz"],
        ) as col_order,
    ):
        database = pd.DataFrame(columns=["test", "foo", "baz"])
        new_dataset = pd.DataFrame(columns=["test", "foo", "baz"])
        primary_key_cols = ["test", "foo", "baz"]
        primary_key = "UNIQUE_KEY"
        allupdates = pd.DataFrame(columns=col_order.return_value)
        expected_allupdates = pd.concat(
            [allupdates, patch_append_rows.return_value, patch_update.return_value],
            sort=False,
        )

        #  # check if to_delete is False
        results = load.check_database_changes(
            database, new_dataset, primary_key_cols, to_delete
        )
        if to_delete:
            patch_delete.assert_called_once_with("table", "table", primary_key)
        else:
            patch_delete.assert_not_called()
        pd.testing.assert_frame_equal(
            results["allupdates"].sort_index(axis=1).reset_index(drop=True),
            expected_allupdates.sort_index(axis=1).reset_index(drop=True),
        )
        pd.testing.assert_frame_equal(
            results["to_delete_rows"], expected_to_delete_rows
        )


def test_store_database_with_rows_to_update_and_delete(syn):
    """Test store_database function"""
    database_synid = "syn123"
    all_updates = pd.DataFrame(
        {"test": ["test", "test2"], "foo": [1, 3], "baz": ["", 5]}
    )
    to_delete_rows = pd.DataFrame({0: ["3"], 1: ["5"]})
    mock_database_ent = Mock()
    mock_database_ent.name = "test_table"
    # get the table entity
    with (
        patch.object(syn, "get", return_value=mock_database_ent) as patch_get,
        patch.object(Table, "store_rows") as patch_store_rows,
        patch.object(Table, "delete_rows") as patch_delete_rows,
        patch.object(load.logger, "info") as patch_info,
    ):
        load.store_database(syn, database_synid, all_updates, to_delete_rows)
        patch_get.assert_called_once_with(database_synid)
        patch_store_rows.assert_called_once_with(values=all_updates)
        patch_delete_rows.assert_called_once_with(df=to_delete_rows)
        patch_info.assert_has_calls(
            [
                call(
                    f"Upserting {len(all_updates)} rows from {mock_database_ent.name} table"
                ),
                call(
                    f"Deleting {len(to_delete_rows)} rows from {mock_database_ent.name} table"
                ),
            ]
        )


def test_store_database_only_rows_to_update(syn):
    """Test store_database function"""
    database_synid = "syn123"
    all_updates = pd.DataFrame(
        {"test": ["test", "test2"], "foo": [1, 3], "baz": ["", 5]}
    )
    to_delete_rows = pd.DataFrame()
    mock_database_ent = Mock()
    mock_database_ent.name = "test_table"
    # get the table entity
    with (
        patch.object(syn, "get", return_value=mock_database_ent) as patch_get,
        patch.object(Table, "store_rows") as patch_store_rows,
        patch.object(Table, "delete_rows") as patch_delete_rows,
        patch.object(load.logger, "info") as patch_info,
    ):
        load.store_database(syn, database_synid, all_updates, to_delete_rows)
        patch_get.assert_called_once_with(database_synid)
        patch_store_rows.assert_called_once_with(values=all_updates)
        patch_delete_rows.assert_not_called()
        patch_info.assert_called_once_with(
            f"Upserting {len(all_updates)} rows from {mock_database_ent.name} table"
        )


def test_store_database_only_rows_to_delete(syn):
    """Test store_database function"""
    database_synid = "syn123"
    all_updates = pd.DataFrame()
    to_delete_rows = pd.DataFrame({0: ["3"], 1: ["5"]})
    mock_database_ent = Mock()
    mock_database_ent.name = "test_table"
    # get the table entity
    with (
        patch.object(syn, "get", return_value=mock_database_ent) as patch_get,
        patch.object(Table, "store_rows") as patch_store_rows,
        patch.object(Table, "delete_rows") as patch_delete_rows,
        patch.object(load.logger, "info") as patch_info,
    ):
        load.store_database(syn, database_synid, all_updates, to_delete_rows)
        patch_get.assert_called_once_with(database_synid)
        patch_store_rows.assert_not_called()
        patch_delete_rows.assert_called_once_with(df=to_delete_rows)
        patch_info.assert_called_once_with(
            f"Deleting {len(to_delete_rows)} rows from {mock_database_ent.name} table"
        )


def test_store_database_no_rows_to_update_and_delete(syn):
    """Test store_database function"""
    database_synid = "syn123"
    all_updates = pd.DataFrame()
    to_delete_rows = pd.DataFrame()
    mock_database_ent = Mock()
    mock_database_ent.name = "test_table"
    # get the table entity
    with (
        patch.object(syn, "get", return_value=mock_database_ent) as patch_get,
        patch.object(Table, "store_rows") as patch_store_rows,
        patch.object(Table, "delete_rows") as patch_delete_rows,
        patch.object(load.logger, "info") as patch_info,
    ):
        load.store_database(syn, database_synid, all_updates, to_delete_rows)
        patch_get.assert_called_once_with(database_synid)
        patch_store_rows.assert_not_called()
        patch_delete_rows.assert_not_called()
        patch_info.assert_not_called()


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
    test_table_synid = "synZZZ"
    test_data = pd.DataFrame(
        {
            "ROW_ID": ["1"],
            "ROW_VERSION": ["1"],
            "CENTER": ["test-center"],
            "data": [123],
        }
    )
    test_new_data = pd.DataFrame(
        {
            "ROW_ID": ["1"],
            "ROW_VERSION": ["1"],
            "CENTER": ["test-center"],
            "data": [123],
        }
    )

    mock_database_ent = Mock()
    mock_database_ent.primaryKey = "PRIMARY_KEY"
    mock_database = Mock()

    with (
        patch.object(syn, "get", return_value=mock_database_ent),
        patch(
            "synapseclient.models.mixins.table_components.QueryMixin.query_async",
            return_value=test_data,
        ) as patch_table_query,
        patch.object(load, "_update_table") as patch__update_table,
    ):
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
