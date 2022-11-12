"""
This module contains all the functions to extract data
from Synapse
"""

from collections import defaultdict
import datetime
import logging
import os
import time
from typing import List

import synapseclient  # lgtm [py/import-and-import-from]
from synapseclient import Synapse
from synapseclient.core.utils import to_unix_epoch_time
import synapseutils
import pandas as pd

from . import process_functions, process_mutation, toRetract, validate

logger = logging.getLogger(__name__)


def get_center_input_files(
        syn: synapseclient.Synapse, synid: str,
        center: str, process: str = "main", downloadFile: bool = True
    ) -> List[List[synapseclient.Entity]]:
    """Walks through each center's input directory
    to get a list of tuples of center files

    Args:
        syn (synapseclient.Synapse): Synapse connection
        synid (str): Synapse Id of a folder
        center (str): GENIE center name
        process (str, optional): Process type include "main", "mutation".
                                 Defaults to "main".
        downloadFile (bool, optional): Downloads the file. Defaults to True.

    Returns:
        list: List of Synapse entities
    """
    logger.info("GETTING {center} INPUT FILES".format(center=center))
    clinical_pair_name = [
        "data_clinical_supp_sample_{center}.txt".format(center=center),
        "data_clinical_supp_patient_{center}.txt".format(center=center),
    ]

    center_files = synapseutils.walk(syn, synid)
    clinicalpair_entities = []
    prepared_center_file_list = []

    for _, _, entities in center_files:
        for name, ent_synid in entities:
            # This is to remove vcfs from being validated during main
            # processing. Often there are too many vcf files, and it is
            # not necessary for them to be run everytime.
            if name.endswith(".vcf") and process != "mutation":
                continue

            ent = syn.get(ent_synid, downloadFile=downloadFile)

            # HACK: Clinical file can come as two files.
            # The two files need to be merged together which is
            # why there is this format
            if name in clinical_pair_name:
                clinicalpair_entities.append(ent)
                continue

            prepared_center_file_list.append([ent])

    if clinicalpair_entities:
        prepared_center_file_list.append(clinicalpair_entities)

    return prepared_center_file_list


# TODO: Add to transform.py
def _map_name_to_filetype(name: str) -> str:
    """Maps file name to filetype

    Args:
        name (str): File name

    Returns:
        str: filetype
    """
    # By default, set the filetype to be the filename
    filetype = name
    if not name.startswith("meta"):
        if name.startswith("data_clinical_sample"):
            filetype = "sample"
        elif name.endswith("fusions.txt"):
            filetype = "fusion"
        elif name.endswith("CNA.txt"):
            filetype = "cna"
        elif name.endswith(".seg"):
            filetype = "seg"
        elif name == "data_sv.txt":
            filetype = "sv"
    return filetype


def get_file_mapping(syn: synapseclient.Synapse, synid: str) -> dict:
    """Get mapping between Synapse entity name and Synapse ids
    of all entities in a folder

    Args:
        syn (synapseclient.Synapse): Synapse connection
        synid (str): Synapse Id of folder

    Returns:
        dict: mapping between Synapse Entity name and Id
    """
    files = syn.getChildren(synid)
    file_mapping = {_map_name_to_filetype(name=metadata["name"]): metadata["id"]
                    for metadata in files}
    return file_mapping


def find_caselistid(syn: synapseclient.Synapse, parentid: str):
    """
    Search for case_lists folder based on parentId given

    Args:
        syn: Synapse object
        parentid: Synapse Id of Folder or Project

    Returns:
        string: Synapse id of case list
    """
    file_mapping = syn.findEntityId(entity_name, parent=parentid)
    # if case_lists doesn't exist
    if file_mapping.get("case_lists") is None:
        caselist_folder = synapseclient.Folder(name="case_lists", parent=parentid)
        caselistid = syn.store(caselist_folder).id
    else:
        caselistid = file_mapping.get("case_lists")
    return caselistid