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
    syn: synapseclient.Synapse,
    synid: str,
    center: str,
    process: str = "main",
    downloadFile: bool = True,
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
    file_mapping = {
        _map_name_to_filetype(name=metadata["name"]): metadata["id"]
        for metadata in files
    }
    return file_mapping


def get_public_to_consortium_synid_mapping(
    syn: synapseclient.Synapse, release_synid: str
) -> dict:
    """
    Gets the mapping between potential public release names and
    the consortium release folder

    Args:
        syn (synapseclient.Synapse): Synapse connection
        release_synid (str): Release folder fileview

    Returns:
        dict: Mapping between potential public release and consortium
              release synapse id
    """
    # This dict contains the mapping between public release name and
    # consortium release folder
    public_to_consortium_map = dict()
    # release_files = synapseutils.walk(syn, releaseSynId)
    # TODO: fix the database to mapping table
    consortium_release_folders = syn.tableQuery(
        f"SELECT name, id FROM {release_synid} WHERE "
        "name NOT LIKE 'Release %' "
        "and name NOT LIKE '%-public' "
        "and name NOT IN ('case_lists', 'potential_artifacts')"
        "ORDER BY name"
    )
    consortium_release_folders_df = consortium_release_folders.asDataFrame()
    # Get major release version
    consortium_release_folders_df["major_release"] = [
        release.split(".")[0] for release in consortium_release_folders_df["name"]
    ]
    # only keep the latest consortium release for the public release
    consortium_release_folders_df.drop_duplicates(
        "major_release", keep="last", inplace=True
    )

    for _, release_info in consortium_release_folders_df.iterrows():
        major_release = release_info["major_release"]
        # add support for potential patch releases
        for num in [0, 1, 2, 3]:
            # This has to exist because the the first three GENIE releases
            # used semantic versioning
            if release_info["major_release"] in ["0", "1", "2"]:
                public_release_name = f"{int(major_release) + 1}.{num}.0"
                public_to_consortium_map[public_release_name] = release_info["id"]
            else:
                public_release_name = f"{major_release}.{num}-public"
                public_to_consortium_map[public_release_name] = release_info["id"]
    return public_to_consortium_map
