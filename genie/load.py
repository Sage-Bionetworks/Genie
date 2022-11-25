"""
This module contains all the functions that stores data
to Synapse
"""
import logging
import time
from typing import Dict, List

import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError

from . import __version__, extract

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
    start: bool = True
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
