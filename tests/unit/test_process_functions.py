import pandas as pd
import genie
import pytest

DATABASE_DF = pd.DataFrame({
    'UNIQUE_KEY': ['test1', 'test2', 'test3'],
    "test": ['test1', 'test2', 'test3'],
    "foo": [1, 2, 3],
    "baz": [float('nan'), float('nan'), float('nan')]})
DATABASE_DF.index = ['1_3', '2_3', '3_5']


def test_valid__check_valid_df():
    genie.process_functions._check_valid_df(DATABASE_DF, "test")


def test_invalid__check_valid_df():
    with pytest.raises(ValueError):
        genie.process_functions._check_valid_df("foo", "test")
    with pytest.raises(ValueError):
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
