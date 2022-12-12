"""
This module contains all the functions that stores data
to Synapse
"""
import logging
import os
import time
import tempfile
from typing import Dict, List

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
    name: str = None,
    annotations: Dict = None,
    used: List[str] = None,
    version_comment: str = None,
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
    col: List[str] = None,
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
    Updates synapse tables by a row identifier with another
    dataset that has the same number and order of columns

    Args:
        syn (synapseclient.Synaps): Synapse object
        database (pd.DataFrame): Original Data
        new_dataset (pd.DataFrame): New Data
        database_synid (str): Synapse Id of the Synapse table
        primary_key_cols (list): Column(s) that make up the primary key
        to_delete (bool, optional): Delete rows. Defaults to False
    """
    primary_key = "UNIQUE_KEY"
    database = database.fillna("")
    orig_database_cols = database.columns
    col_order = ["ROW_ID", "ROW_VERSION"]
    col_order.extend(orig_database_cols.tolist())
    new_dataset = new_dataset.fillna("")
    # Columns must be in the same order
    new_dataset = new_dataset[orig_database_cols]
    database[primary_key_cols] = database[primary_key_cols].applymap(str)
    database[primary_key] = database[primary_key_cols].apply(
        lambda x: " ".join(x), axis=1
    )

    new_dataset[primary_key_cols] = new_dataset[primary_key_cols].applymap(str)
    new_dataset[primary_key] = new_dataset[primary_key_cols].apply(
        lambda x: " ".join(x), axis=1
    )

    allupdates = pd.DataFrame(columns=col_order)
    to_append_rows = process_functions._append_rows(new_dataset, database, primary_key)
    to_update_rows = process_functions._update_rows(new_dataset, database, primary_key)
    if to_delete:
        to_delete_rows = process_functions._delete_rows(
            new_dataset, database, primary_key
        )
    else:
        to_delete_rows = pd.DataFrame()
    allupdates = pd.concat([allupdates, to_append_rows, to_update_rows], sort=False)
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
