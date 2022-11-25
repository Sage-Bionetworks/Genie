"""
This module contains all the functions that stores data
to Synapse
"""
import logging
from typing import Dict, List

import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError

logger = logging.getLogger(__name__)


# TODO make sure to remove from database_to_staging
# def store_file(
#     syn,
#     filePath,
#     genieVersion="database",
#     name=None,
#     parent=None,
#     fileFormat=None,
#     cBioFileFormat=None,
#     tag_or_commit=None,
#     used=None,
# ):
#     """
#     Convenience function to store files

#     Args:
#         syn: Synapse object
#         filePath: path to file to store
#         genieVersion: Version of genie release
#         name: Name of entity
#         fileFormat: GENIE file format
#         cBioFileFormat: cBioPortal file format
#         staging: Staging GENIE release.  Default to False
#         caseLists: Case lists are stored elsewhere
#         tag_or_commit: Github tag or commit
#         used: List of entities used in creation of file
#     """
#     logger.info("STORING FILE: {}".format(os.path.basename(filePath)))
#     if name is None:
#         name = os.path.basename(filePath)
#     ent = synapseclient.File(
#         filePath, name=name, parent=parent, versionComment=genieVersion
#     )
#     if fileFormat is not None:
#         ent.fileFormat = fileFormat
#     if cBioFileFormat is not None:
#         ent.cBioFileFormat = cBioFileFormat
#     if tag_or_commit is None:
#         tag_or_commit = f"v{__version__}"
#     ent = syn.store(
#         ent,
#         executed=f"https://github.com/Sage-Bionetworks/Genie/tree/{tag_or_commit}",
#         used=used,
#     )
#     return ent


def store_file(
    syn: synapseclient.Synapse, filepath: str, parentid: str,
    annotations: Dict = None, used: List[str] = None, executed: List[str] = None
) -> synapseclient.File:
    """Stores file into Synapse

    Args:
        syn (synapseclient.Synapse): Synapse connection
        filepath (str): Path to file
        parentid (str): Synapse Id of a folder or project

    Returns:
        synapseclient.File: Synapse File entity
    """
    file_ent = synapseclient.File(filepath, parentId=parentid)
    if annotations is not None:
        file_ent.annotations = annotations
    file_ent = syn.store(file_ent, used=used, executed=executed)
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
