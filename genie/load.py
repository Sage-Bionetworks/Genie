"""
This module contains all the functions that stores data
to Synapse
"""

import logging
import os
import time
import tempfile
from typing import Dict, List, Optional, Union

import pandas as pd
import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError

from . import __version__, extract, process_functions

logger = logging.getLogger(__name__)


# TODO Edit docstring
def store_file(
    syn: synapseclient.Synapse,
    filepath: str,
    parentid: str,
    name: Optional[str] = None,
    annotations: Optional[Dict] = None,
    used: Optional[Union[List[str], str]] = None,
    version_comment: Optional[str] = None,
) -> synapseclient.File:
    """Stores file into Synapse

    Args:
        syn (synapseclient.Synapse): Synapse connection
        filepath (str): Path to file
        parentid (str): Synapse Id of a folder or project
        name (str, optional): Name of entity. Defaults to basename of your file path.
        annotations (Dict, optional): Synapse annotations to the File Entity. Defaults to None.
        used (List[str], optional): Entities used to generate file. Defaults to None.
        version_comment (str, optional): File Entity version comment. Defaults to None.

    Returns:
        synapseclient.File: Synapse File entity
    """
    file_ent = synapseclient.File(
        filepath, parentId=parentid, versionComment=version_comment
    )
    if name is not None:
        file_ent.name = name
    if annotations is not None:
        file_ent.annotations = annotations
    file_ent = syn.store(
        file_ent,
        used=used,
        executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
    )
    return file_ent


def store_files(
    syn: synapseclient.Synapse, filepaths: List[str], parentid: str
) -> List[synapseclient.File]:
    """Stores a list of files

    Args:
        syn (synapseclient.Synapse): Synapse connection
        filepaths (List[str]): List of filepaths
        parentid (str): Synapse Id of a folder or project

    Returns:
        List[synapseclient.File]: List of Synaps File entities
    """
    file_entities = [
        store_file(syn=syn, filepath=path, parentid=parentid) for path in filepaths
    ]
    return file_entities


def store_table(syn: synapseclient.Synapse, filepath: str, tableid: str):
    """Stores a tsv to a Synapse Table.  This can append, update, and delete
    rows in a Synapse table depending on how the tsv file is formatted.

    Args:
        syn (synapseclient.Synapse): Synapse connection
        filepath (str): Path to a TSV
        tableid (str): Synapse Id of a Synapse Table

    """
    try:
        update_table = synapseclient.Table(tableid, filepath, separator="\t")
        syn.store(update_table)
    except SynapseTimeoutError:
        # This error occurs because of waiting for table to index.
        # Don't worry about this error
        pass


def update_process_trackingdf(
    syn: synapseclient.Synapse,
    process_trackerdb_synid: str,
    center: str,
    process_type: str,
    start: bool = True,
):
    """
    Updates the processing tracking database

    Args:
        syn (synapseclient.Synapse): Synapse connection
        process_trackerdb_synid: Synapse Id of Process Tracking Table
        center: GENIE center (ie. SAGE)
        process_type: processing type (dbToStage or public)
        start: Start or end of processing.  Default is True for start
    """
    logger.info("UPDATE PROCESS TRACKING TABLE")
    column = "timeStartProcessing" if start else "timeEndProcessing"
    query_str = (
        f"SELECT {column} FROM {process_trackerdb_synid} where center = '{center}' "
        f"and processingType = '{process_type}'"
    )
    process_trackerdf = extract.get_syntabledf(syn=syn, query_string=query_str)
    process_trackerdf[column].iloc[0] = str(int(time.time() * 1000))
    syn.store(synapseclient.Table(process_trackerdb_synid, process_trackerdf))


def update_table(
    syn: synapseclient.Synapse,
    databaseSynId: str,
    newData: pd.DataFrame,
    filterBy: str,
    filterByColumn: str = "CENTER",
    col: Optional[List[str]] = None,
    toDelete: bool = False,
):
    """Update Synapse table given a new dataframe

    Args:
        syn (Synapse): Synapse connection
        databaseSynId (str): Synapse Id of Synapse Table
        newData (pd.DataFrame): New data in a dataframe
        filterBy (str): Value to filter new data by
        filterByColumn (str, optional): Column to filter values by. Defaults to "CENTER".
        col (List[str], optional): List of columns to ingest. Defaults to None.
        toDelete (bool, optional): Delete rows given the primary key. Defaults to False.
    """
    databaseEnt = syn.get(databaseSynId)
    database = syn.tableQuery(
        f"SELECT * FROM {databaseSynId} where {filterByColumn} ='{filterBy}'"
    )
    database = database.asDataFrame()
    db_cols = set(database.columns)
    if col is not None:
        new_data_cols = set(col)
        # Make sure columns from file exists in database columns
        use_cols = db_cols.intersection(new_data_cols)
        # No need to fail, because there is bound to be at least one
        # column that will exist in the database
        database = database[list(use_cols)]
    else:
        newData = newData[database.columns]
    _update_table(
        syn, database, newData, databaseSynId, databaseEnt.primaryKey, toDelete
    )


def _update_table(
    syn: synapseclient.Synapse,
    database: pd.DataFrame,
    new_dataset: pd.DataFrame,
    database_synid: str,
    primary_key_cols: List[str],
    to_delete: bool = False,
):
    """
    A helper function to compare new dataset with existing data,
    and store any changes that need to be made to the database
    """
    changes = check_database_changes(database, new_dataset, primary_key_cols, to_delete)
    store_database(
        syn,
        database_synid,
        changes["col_order"],
        changes["allupdates"],
        changes["to_delete_rows"],
    )


def _get_col_order(orig_database_cols: pd.Index) -> List[str]:
    """
    Get column order

    Args:
        orig_database_cols (pd.Index): A list of column names of the original database

    Returns:
        The list of re-ordered column names
    """
    col_order = ["ROW_ID", "ROW_VERSION"]
    col_order.extend(orig_database_cols.tolist())
    return col_order


def _reorder_new_dataset(
    orig_database_cols: pd.Index, new_dataset: pd.DataFrame
) -> pd.DataFrame:
    """
    Reorder new dataset based on the original database

    Args:
        orig_database_cols (pd.Index): A list of column names of the original database
        new_dataset(pd.DataFrame): New Data

    Returns:
        The re-ordered new dataset
    """
    # Columns must be in the same order as the original data
    new_dataset = new_dataset[orig_database_cols]
    return new_dataset


def _generate_primary_key(
    dataset: pd.DataFrame, primary_key_cols: List[str], primary_key: str
) -> pd.DataFrame:
    """
    Generate primary key column a dataframe

    Args:
        dataset(pd.DataFrame): A dataframe
        new_dataset: The re-ordered new dataset
        primary_key_cols (list): Column(s) that make up the primary key
        primary_key: The column name of the primary_key
    Returns:
        The dataframe with primary_key column added
    """
    # replace NAs with emtpy string
    dataset = dataset.fillna("")
    # generate primary key column for original database
    dataset[primary_key_cols] = dataset[primary_key_cols].applymap(str)
    if dataset.empty:
        dataset[primary_key] = ""
    else:
        dataset[primary_key] = dataset[primary_key_cols].apply(
            lambda x: " ".join(x), axis=1
        )
    return dataset


def check_database_changes(
    database: pd.DataFrame,
    new_dataset: pd.DataFrame,
    primary_key_cols: List[str],
    to_delete: bool = False,
) -> Dict[pd.DataFrame, List[str]]:
    """
    Check changes that need to be made, i.e. append/update/delete rows to the database
    based on its comparison with new data

    Args:
        database (pd.DataFrame): Original Data
        new_dataset (pd.DataFrame): New Data
        primary_key_cols (list): Column(s) that make up the primary key
        to_delete (bool, optional): Delete rows. Defaults to False
    """
    # get a list of column names of the original database
    orig_database_cols = database.columns
    # get the final column order
    col_order = _get_col_order(orig_database_cols)
    # reorder new_dataset
    new_dataset = _reorder_new_dataset(orig_database_cols, new_dataset)
    # set the primary_key name
    primary_key = "UNIQUE_KEY"
    # generate primary_key column for dataset comparison
    ori_data = _generate_primary_key(database, primary_key_cols, primary_key)
    new_data = _generate_primary_key(new_dataset, primary_key_cols, primary_key)
    # output dictionary
    changes = {"col_order": col_order, "allupdates": None, "to_delete_rows": None}
    # get rows to be appened or updated
    allupdates = pd.DataFrame(columns=col_order)
    to_append_rows = process_functions._append_rows(new_data, ori_data, primary_key)
    to_update_rows = process_functions._update_rows(new_data, ori_data, primary_key)
    allupdates = pd.concat([allupdates, to_append_rows, to_update_rows], sort=False)
    changes["allupdates"] = allupdates
    # get rows to be deleted
    if to_delete:
        to_delete_rows = process_functions._delete_rows(new_data, ori_data, primary_key)
    else:
        to_delete_rows = pd.DataFrame()
    changes["to_delete_rows"] = to_delete_rows
    return changes


def store_database(
    syn: synapseclient.Synapse,
    database_synid: str,
    col_order: List[str],
    allupdates: pd.DataFrame,
    to_delete_rows: pd.DataFrame,
):
    """
    Store changes to the database

    Args:
        syn (synapseclient.Synaps): Synapse object
        database_synid (str): Synapse Id of the Synapse table
        col_order (List[str]): The ordered column names to be saved
        allupdates (pd.DataFrame): rows to be appended and/or updated
        to_deleted_rows (pd.DataFrame): rows to be deleted

    Returns:
        None
    """
    storedatabase = False
    update_all_file = tempfile.NamedTemporaryFile(
        dir=process_functions.SCRIPT_DIR, delete=False
    )
    with open(update_all_file.name, "w") as updatefile:
        # Must write out the headers in case there are no appends or updates
        updatefile.write(",".join(col_order) + "\n")
        if not allupdates.empty:
            """
            This is done because of pandas typing.
            An integer column with one NA/blank value
            will be cast as a double.
            """
            updatefile.write(
                allupdates[col_order]
                .to_csv(index=False, header=None)
                .replace(".0,", ",")
                .replace(".0\n", "\n")
            )
            storedatabase = True
        if not to_delete_rows.empty:
            updatefile.write(
                to_delete_rows.to_csv(index=False, header=None)
                .replace(".0,", ",")
                .replace(".0\n", "\n")
            )
            storedatabase = True
    if storedatabase:
        syn.store(synapseclient.Table(database_synid, update_all_file.name))
    # Delete the update file
    os.unlink(update_all_file.name)
