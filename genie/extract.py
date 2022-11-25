"""
This module contains all the functions that extract data
from Synapse
"""

import logging
from typing import List

import synapseclient
import synapseutils
import pandas as pd

from genie import process_functions

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


def get_syntabledf(syn: synapseclient.Synapse, query_string: str) -> pd.DataFrame:
    """
    Get dataframe from Synapse Table query

    Args:
        syn (synapseclient.Synapse): Synapse connection
        query_string (str): Table query

    Returns:
        pd.DataFrame: Query results in a dataframe
    """
    table = syn.tableQuery(query_string)
    tabledf = table.asDataFrame()
    return tabledf


def _get_synid_database_mappingdf(syn, project_id):
    """
    Get database to synapse id mapping dataframe

    Args:
        syn: Synapse object
        project_id: Synapse Project ID with a 'dbMapping' annotation.

    Returns:
        database to synapse id mapping dataframe
    """

    project = syn.get(project_id)
    database_mapping_synid = project.annotations["dbMapping"][0]
    database_map_query = f"SELECT * FROM {database_mapping_synid}"
    mappingdf = get_syntabledf(syn, database_map_query)
    return mappingdf


def getDatabaseSynId(syn, tableName, project_id=None, databaseToSynIdMappingDf=None):
    """
    Get database synapse id from database to synapse id mapping table

    Args:
        syn: Synapse object
        project_id: Synapse Project ID with a database mapping table.
        tableName: Name of synapse table
        databaseToSynIdMappingDf: Avoid calling rest call to download table
                                  if the mapping table is already downloaded

    Returns:
        str:  Synapse id of wanted database
    """
    if databaseToSynIdMappingDf is None:
        databaseToSynIdMappingDf = _get_synid_database_mappingdf(
            syn, project_id=project_id
        )

    synId = process_functions.lookup_dataframe_value(
        databaseToSynIdMappingDf, "Id", f'Database == "{tableName}"'
    )
    return synId


def _get_database_mapping_config(syn: synapseclient.Synapse, synid: str) -> dict:
    """Gets Synapse database to Table mapping in dict

    Args:
        syn (synapseclient.Synapse): Synapse connection
        synid (str): Synapse id of database mapping table

    Returns:
        dict: {'databasename': 'synid'}
    """
    configdf = get_syntabledf(syn=syn, query_string=f"SELECT * FROM {synid}")
    configdf.index = configdf["Database"]
    config_dict = configdf.to_dict()
    return config_dict["Id"]


def get_genie_config(
    syn: synapseclient.Synapse,
    project_id: str,
) -> dict:
    """Get configurations needed for the GENIE codebase

    Args:
        syn (synapseclient.Synapse): Synapse connection
        project_id (str): Synapse project id

    Returns:
        dict: GENIE table type/name to Synapse Id
    """
    # Get the Synapse Project where data is stored
    # Should have annotations to find the table lookup
    project = syn.get(project_id)

    # Get project GENIE configurations
    database_to_synid_mapping_synid = project.annotations.get("dbMapping", "")
    genie_config = _get_database_mapping_config(
        syn=syn, synid=database_to_synid_mapping_synid[0]
    )
    # Fill in GENIE center configurations
    center_mapping_id = genie_config["centerMapping"]
    center_mapping_df = get_syntabledf(
        syn=syn, query_string=f"SELECT * FROM {center_mapping_id} where release is true"
    )
    center_mapping_df.index = center_mapping_df.center
    # Add center configurations including input/staging synapse ids
    genie_config["center_config"] = center_mapping_df.to_dict("index")

    genie_config["ethnicity_mapping"] = "syn7434242"
    genie_config["race_mapping"] = "syn7434236"
    genie_config["sex_mapping"] = "syn7434222"
    genie_config["sampletype_mapping"] = "syn7434273"

    return genie_config


# TODO: Remove oncotree_link parameter from this function
def _get_oncotreelink(
    syn: synapseclient.Synapse, genie_config: dict, oncotree_link: str = None
) -> str:
    """
    Gets oncotree link unless a link is specified by the user

    Args:
        syn (synapseclient.Synapse): Synapse connection
        genie_config (dict): database name to synid mapping
        oncotree_link (str): link to oncotree. Default is None

    Returns:
        str: oncotree link
    """
    if oncotree_link is None:
        onco_link_ent = syn.get(genie_config["oncotreeLink"])
        oncotree_link = onco_link_ent.externalURL
    return oncotree_link
