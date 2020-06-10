from unittest import mock
from unittest.mock import patch

import pandas as pd
import pytest
import synapseclient

import genie.process_functions

syn = mock.create_autospec(synapseclient.Synapse)

DATABASE_DF = pd.DataFrame({
    'UNIQUE_KEY': ['test1', 'test2', 'test3'],
    "test": ['test1', 'test2', 'test3'],
    "foo": [1, 2, 3],
    "baz": [float('nan'), float('nan'), float('nan')]})
DATABASE_DF.index = ['1_3', '2_3', '3_5']
ENTITY = synapseclient.Project("foo", annotations={"dbMapping": ["syn1234"]})

@pytest.mark.parametrize("input_str,output", [
        ("1.0\t", "1\t"),
        ("1.0\n", "1\n"),
        ("1.5\t", "1.5\t"),
        ("1\t", "1\t"),
        ("0\t", "0\t"),
        ("'a'\t'b'\n1.0\t2.0\n", "'a'\t'b'\n1\t2\n"),
    ])
def test_removeStringFloat(input_str, output):
    """Remove string float - will always assume that there is a \n
    at the end.  This is because if a value was 2.01, we dont want to
    remove the .0 from this."""
    assert genie.process_functions.removeStringFloat(input_str) == output


def test_valid__check_valid_df():
    genie.process_functions._check_valid_df(DATABASE_DF, "test")


def test_invalid__check_valid_df():
    with pytest.raises(ValueError, match="Must pass in pandas dataframe"):
        genie.process_functions._check_valid_df("foo", "test")
    with pytest.raises(
            ValueError,
            match="'error' column must exist in dataframe"):
        genie.process_functions._check_valid_df(DATABASE_DF, "error")


def test__get_left_diff_df():
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test1', 'test2', 'test3', 'test4'],
        "test": ['test1', 'test2', 'test3', 'test4'],
        "foo": [1, 2, 3, 4],
        "baz": [float('nan'), float('nan'), float('nan'), 3.2]})
    get_diff = genie.process_functions._get_left_diff_df(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    expecteddf = new_datadf.loc[[3]]
    assert get_diff.equals(expecteddf[get_diff.columns])


def test_norows_get_left_diff_df():
    append_rows = genie.process_functions._get_left_diff_df(
        DATABASE_DF, DATABASE_DF, 'UNIQUE_KEY')
    assert append_rows.empty


def test_first_validation_get_left_diff_df():
    '''
    This checks to make sure that validation is called
    - In a situation where someone comments out the validation
    line, this will cause an error
    '''
    with pytest.raises(
            ValueError,
            match="'FOO' column must exist in dataframe"):
        genie.process_functions._get_left_diff_df(
            DATABASE_DF, DATABASE_DF, 'FOO')


def test_second_validation_get_left_diff_df():
    '''
    This checks to make sure that validation is called
    - In a situation where someone comments out the validation
    line, this will cause an error
    '''
    testing = DATABASE_DF.copy()
    testing['FOO'] = float('nan')
    with pytest.raises(
            ValueError,
            match="'FOO' column must exist in dataframe"):
        genie.process_functions._get_left_diff_df(
            testing, DATABASE_DF, 'FOO')


def test_first_validation_get_left_union_df():
    '''
    This checks to make sure that validation is called
    - In a situation where someone comments out the 1st validation
    line, this will cause an error
    '''
    with pytest.raises(
            ValueError,
            match="'FOO' column must exist in dataframe"):
        genie.process_functions._get_left_union_df(
            DATABASE_DF, DATABASE_DF, 'FOO')


def test_second_validation_get_left_union_df():
    '''
    This checks to make sure that validation is called
    - In a situation where someone comments out the 2nd validation
    line, this will cause an error
    '''
    testing = DATABASE_DF.copy()
    testing['FOO'] = float('nan')
    with pytest.raises(
            ValueError,
            match="'FOO' column must exist in dataframe"):
        genie.process_functions._get_left_union_df(
            testing, DATABASE_DF, 'FOO')


def test_append__append_rows():
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test1', 'test2', 'test3', 'test4'],
        "test": ['test1', 'test2', 'test3', 'test4'],
        "foo": [1, 2, 3, 4],
        "baz": [float('nan'), float('nan'), float('nan'), 3.2]})
    expecteddf = pd.DataFrame({
        'test': ['test4'],
        'foo': [4],
        'baz': [3.2]})
    append_rows = genie.process_functions._append_rows(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    append_rows.fillna('', inplace=True)
    expecteddf.fillna('', inplace=True)
    assert append_rows.equals(expecteddf[append_rows.columns])


def test___create_update_rowsdf():
    differentrows = [True, True, False]
    database = pd.DataFrame({
        "test": ['test', 'test2', 'test3'],
        "foo": [1, 3, 3],
        "baz": [float('nan'), 5, float('nan')]},
        index=['test1', 'test5', 'test4'])
    new_datadf = pd.DataFrame({
        "test": ['test1', 'test4', 'test3'],
        "foo": [2, 3, 3],
        "baz": [3, 5, float('nan')]},
        index=['test1', 'test5', 'test4'])

    to_update_rowsdf = genie.process_functions._create_update_rowsdf(
        database, new_datadf, DATABASE_DF.index, differentrows)
    expecteddf = pd.DataFrame({
        "test": ['test1', 'test4'],
        "foo": [2, 3],
        "baz": [3.0, 5.0],
        "ROW_ID": ["1", "2"],
        "ROW_VERSION": ["3", "3"]})
    assert to_update_rowsdf.equals(expecteddf[to_update_rowsdf.columns])


def test_none__create_update_rowsdf():
    differentrows = [False, False, False]
    database = pd.DataFrame({
        "test": ['test', 'test2', 'test3'],
        "foo": [1, 3, 3],
        "baz": [float('nan'), 5, float('nan')]},
        index=['test1', 'test5', 'test4'])
    new_datadf = pd.DataFrame({
        "test": ['test1', 'test4', 'test3'],
        "foo": [2, 3, 3],
        "baz": [3, 5, float('nan')]},
        index=['test1', 'test5', 'test4'])

    to_update_rowsdf = genie.process_functions._create_update_rowsdf(
        database, new_datadf, DATABASE_DF.index, differentrows)
    assert to_update_rowsdf.empty


def test___get_left_union_df():
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test1', 'test5', 'test4'],
        "test": ['test', 'test2', 'test3'],
        "foo": [1, 3, 3],
        "baz": [float('nan'), 5, float('nan')]})
    left_union = genie.process_functions._get_left_union_df(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    expecteddf = pd.DataFrame({
        'UNIQUE_KEY': ['test1'],
        'test': ['test'],
        'foo': [1],
        'baz': [float('nan')]})
    assert left_union.equals(expecteddf[left_union.columns])


def test_none__get_left_union_df():
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test7', 'test5', 'test4'],
        "test": ['test', 'test2', 'test3'],
        "foo": [1, 3, 3],
        "baz": [float('nan'), 5, float('nan')]})
    left_union = genie.process_functions._get_left_union_df(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    assert left_union.empty


def test_update__update_rows():
    '''
    Tests index comparison for updating rows
    '''
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test1', 'test2', 'test3'],
        "test": ['test', 'test2', 'test3'],
        "foo": [1, 3, 3],
        "baz": [float('nan'), 5, float('nan')]})

    expecteddf = pd.DataFrame({
        "test": ['test', 'test2'],
        "foo": [1, 3],
        "baz": ['', 5],
        'ROW_ID': ['1', '2'],
        'ROW_VERSION': ['3', '3']})
    update_rows = genie.process_functions._update_rows(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    assert update_rows.equals(expecteddf[update_rows.columns])


def test_maintaintype__update_rows():
    '''
    Test pandas behavior.  Integer -> Float if NA exists
    '''
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test1', 'test2', 'test3'],
        "test": ['test1', 'test2', 'test3'],
        "foo": [1, 3, 3],
        "baz": [float('nan'), 5, float('nan')]})
    # Test that the datatype passed into from new_datadf gets preserved
    expecteddf = pd.DataFrame({
        "test": ['test2'],
        "foo": [3],
        "baz": [5],
        'ROW_ID': ['2'],
        'ROW_VERSION': ['3']})
    expecteddf = expecteddf.astype({'baz': object})
    update_rows = genie.process_functions._update_rows(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    assert update_rows.equals(expecteddf[update_rows.columns])


def test_noupdate__update_rows():
    '''
    Tests the index comparison to get no updates
    '''
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test4'],
        "test": ['test'],
        "foo": [1],
        "baz": [float('nan')]})
    update_rows = genie.process_functions._update_rows(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    assert update_rows.empty


def test_delete__delete_rows():
    new_datadf = pd.DataFrame({
        'UNIQUE_KEY': ['test1'],
        "test": ['test1'],
        "foo": [1],
        "baz": [float('nan')]})
    expecteddf = pd.DataFrame({
        0: ['2', '3'],
        1: ['3', '5']})
    delete_rows = genie.process_functions._delete_rows(
        new_datadf, DATABASE_DF, 'UNIQUE_KEY')
    assert delete_rows.equals(expecteddf)


def test_norows__delete_rows():
    delete_rows = genie.process_functions._delete_rows(
        DATABASE_DF, DATABASE_DF, 'UNIQUE_KEY')
    assert delete_rows.empty


class argparser:
    def asDataFrame(self):
        database_dict = {"Database": ["centerMapping"],
                         "Id": ["syn123"]}
        databasetosynid_mappingdf = pd.DataFrame(database_dict)
        return(databasetosynid_mappingdf)


def test_get_synid_database_mappingdf():
    '''
    Test getting database mapping config
    no flags
    staging flag
    test flag
    '''
    arg = argparser()
    with patch.object(syn, "get", return_value=ENTITY), \
         patch.object(genie.process_functions, "get_syntabledf",
                      return_value=arg.asDataFrame()) as patch_gettabledf:
        df = genie.process_functions.get_synid_database_mappingdf(
            syn, project_id=None)
        patch_gettabledf.assert_called_once_with(
            syn, "SELECT * FROM {}".format(ENTITY.dbMapping[0]))
        assert df.equals(arg.asDataFrame())


def test_get_syntabledf():
    '''
    Test helper function that queries synapse tables and returns dataframes
    '''
    arg = argparser()
    with patch.object(syn, "tableQuery",
                      return_value=arg) as patch_syn_tablequery:
        querystring = "select * from foo"
        df = genie.process_functions.get_syntabledf(syn, querystring)
        patch_syn_tablequery.assert_called_once_with(querystring)
        assert df.equals(arg.asDataFrame())
