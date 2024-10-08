#!/usr/bin/env python3
from collections import defaultdict
import datetime
import logging
import os
import time
from typing import List, Optional

import synapseclient
from synapseclient import Synapse
from synapseclient.core.utils import to_unix_epoch_time
import pandas as pd

from genie import (
    extract,
    load,
    process_functions,
    process_mutation,
    toRetract,
    validate,
)

logger = logging.getLogger(__name__)


DUPLICATED_FILE_ERROR = (
    "Duplicated filename! Files should be uploaded as new versions "
    "and the entire dataset should be uploaded."
)

"""
TODO:
Could potentially get all the inforamation of the file entity right here
To avoid the syn.get rest call later which doesn't actually download the file
"""


# TODO: add to transform.py
def entity_date_to_timestamp(entity_date_time):
    """Convert Synapse object date/time string (from modifiedOn or createdOn properties) to a timestamp."""

    date_and_time = entity_date_time.split(".")[0]
    date_time_obj = datetime.datetime.strptime(date_and_time, "%Y-%m-%dT%H:%M:%S")
    return to_unix_epoch_time(date_time_obj)


# TODO: Add to validation.py
def check_existing_file_status(validation_status_table, error_tracker_table, entities):
    """
    This function checks input files against the existing validation and error
    tracking dataframe

    Args:
        validation_status_table: Validation status Synapse Table query result
        error_tracker_table: Error tracking Synapse Table query result
        entities: list of center input entites

    Returns:
        dict: Input file status
            status_list: file validation status
            error_list: Errors of the files if they exist,
            to_validate: Boolean value for whether of not an input
                         file needs to be validated
    """
    if len(entities) > 2:
        raise ValueError("There should never be more than 2 files being validated.")

    statuses = []
    errors = []

    validation_statusdf = validation_status_table.asDataFrame()
    error_trackerdf = error_tracker_table.asDataFrame()
    # This should be outside fo the forloop so that it doesn't
    # get reset
    to_validate = False
    for ent in entities:
        # Get the current status and errors from the tables.
        current_status = validation_statusdf[validation_statusdf["id"] == ent.id]
        current_error = error_trackerdf[error_trackerdf["id"] == ent.id]

        if current_status.empty:
            to_validate = True
        else:
            # This to_validate is here, because the following is a
            # sequential check of whether files need to be validated
            statuses.append(current_status["status"].values[0])
            if current_error.empty:
                to_validate = current_status["status"].values[0] == "INVALID"
            else:
                errors.append(current_error["errors"].values[0])
            # Add Name check here (must add name of the entity as a column)
            if (
                current_status["md5"].values[0] != ent.md5
                or current_status["name"].values[0] != ent.name
            ):
                to_validate = True
            else:
                status_str = "{filename} ({id}) FILE STATUS IS: {filestatus}"
                logger.info(
                    status_str.format(
                        filename=ent.name,
                        id=ent.id,
                        filestatus=current_status["status"].values[0],
                    )
                )

    return {"status_list": statuses, "error_list": errors, "to_validate": to_validate}


# TODO: Add to validation.py
def _send_validation_error_email(syn, user, message_objs):
    """
    Sends validation error email

    Args:
        syn: Synapse object
        user: username to send message to
        message_objs: list of dicts with 'filenames' and 'messages' to send
    """

    username = syn.getUserProfile(user)["userName"]

    errors = ""
    for message_obj in message_objs:
        file_names = ", ".join(message_obj["filenames"])
        error_message = message_obj["messages"]
        errors += f"Filenames: {file_names}, Errors:\n {error_message}\n\n"

    email_message = (
        f"Dear {username},\n\n"
        "You have invalid files! "
        f"Here are the reasons why:\n\n{errors}"
    )

    date_now = datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    syn.sendMessage(
        userIds=[user],
        messageSubject=f"GENIE Validation Error - {date_now}",
        messageBody=email_message,
    )


# TODO: Add to validation.py
def _get_status_and_error_list(valid, message, entities):
    """
    Helper function to return the status and error list of the
    files based on validation result.

    Args:
        valid: Boolean value of results of validation
        message: Validation message
        entities: List of Synapse Entities

    Returns:
        tuple: input_status_list - status of input files list,
               invalid_errors_list - error list
    """
    if valid:
        input_status_list = [{"entity": ent, "status": "VALIDATED"} for ent in entities]
        invalid_errors_list = []
    else:
        input_status_list = [{"entity": ent, "status": "INVALID"} for ent in entities]
        invalid_errors_list = [{"entity": ent, "errors": message} for ent in entities]
    return input_status_list, invalid_errors_list


# def get_ancillary_files(
#     syn: synapseclient.Synapse,
#     synid: str,
#     project_id: str,
#     center: str,
#     process: str = "main",
#     downloadFile: bool = True,
#     genie_config: list = None,
#     format_registry: list = None,
# ) -> Dict[str, Dict[str, object]]:
#     """Walks through each center's input directory
#     to get a dict of center files

#     Args:
#         syn (synapseclient.Synapse): Synapse connection
#         synid (str): Synapse Id of a folder
#         project_id (str): GENIE Synapse project id
#         center (str): GENIE center name
#         process (str, optional): Process type include "main", "mutation".
#                                  Defaults to "main".
#         downloadFile (bool, optional): Downloads the file. Defaults to True.

#     Returns:
#         dict: {entity_name: {
#             entity: Synapse.File
#             filetypeformat_object: FileTypeFormat
#             }
#     """
#     logger.info("GETTING {center} INPUT FILES".format(center=center))
#     clinical_pair_name = [
#         "data_clinical_supp_sample_{center}.txt".format(center=center),
#         "data_clinical_supp_patient_{center}.txt".format(center=center),
#     ]
#     clinicalpair_entities = []

#     center_files = synapseutils.walk(syn, synid)
#     prepared_center_files = {}

#     for _, _, entities in center_files:
#         for name, ent_synid in entities:
#             if name in clinical_pair_name:
#                 clinicalpair_entities.append(ent)
#                 continue

#             if name.endswith(".vcf") and process != "mutation":
#                 continue

#             ent = syn.get(ent_synid, downloadFile=downloadFile)

#             validator = validate.GenieValidationHelper(
#                 syn=syn,
#                 project_id=project_id,
#                 center=center,
#                 entitylist=[ent],
#                 format_registry=format_registry,
#                 genie_config=genie_config,
#                 ancillary_files=None,
#             )
#             filetype = validator.file_type
#             if validator.file_type not in validator._format_registry:
#                 continue
#             validator_cls = validator._format_registry[validator.file_type]
#             fileformat_validator = validator_cls(
#                 syn=validator._synapse_client,
#                 center=validator.center,
#                 genie_config=validator.genie_config,
#                 ancillary_files=None,
#             )

#             prepared_center_files[name] = {}
#             prepared_center_files[name]["entity"] = ent
#             prepared_center_files[name]["filetypeformat_object"] = fileformat_validator

#     # if the clinical files exist
#     if clinicalpair_entities:
#         # handling for just the clinical pair, can remove once we have separate classes
#         cli_validator = validate.GenieValidationHelper(
#             syn=syn,
#             project_id=project_id,
#             center=center,
#             entitylist=clinicalpair_entities,
#             format_registry=format_registry,
#             genie_config=genie_config,
#             ancillary_files=None,
#         )
#         cli_filetype = cli_validator.file_type
#         cli_validator_cls = cli_validator._format_registry[cli_validator.file_type]
#         cli_fileformat_validator = validator_cls(
#             syn=cli_validator._synapse_client,
#             center=cli_validator.center,
#             genie_config=cli_validator.genie_config,
#             ancillary_files=None,
#         )

#         prepared_center_files[name] = {}
#         prepared_center_files[name]["entity"] = clinicalpair_entities
#         prepared_center_files[name]["filetypeformat_object"] = cli_fileformat_validator
#     return prepared_center_files


# TODO: Add to validation.py
def validatefile(
    syn: synapseclient.Synapse,
    project_id: str,
    entities: List[synapseclient.File],
    validation_status_table: synapseclient.table.CsvFileTable,
    error_tracker_table: synapseclient.table.CsvFileTable,
    center: str,
    format_registry: Optional[dict] = None,
    genie_config: Optional[dict] = None,
    ancillary_files: Optional[list] = None,
):
    """Validate a list of entities.

    If a file has not changed, then it doesn't need to be validated.

    Args:
        syn: Synapse object
        project_id (str): GENIE Synapse project id
        entities: A list of entities for a single file 'type' (usually a single file, but clinical can have two)
        validation_status_table: Validation status dataframe
        error_tracker_table: Invalid files error tracking dataframe
        center: Center of interest
        format_registry (list, optional): GENIE file format registry.
                                                 Defaults to None.
        genie_config (list, optional): See example of genie config at
                                              ./genie_config.json. Defaults to None.
        ancillary_files (list): all files downloaded for validation

    Returns:
        tuple: input_status_list - status of input files,
               invalid_errors_list - error list
               messages_to_send - list of tuples with (filenames, message, file_users)

    """
    # TODO: Look into if errors should be thrown if these are None
    # Aka. should these actually be optional params
    if genie_config is None:
        genie_config = {}
    if format_registry is None:
        format_registry = {}

    # filepaths = [entity.path for entity in entities]
    filenames = [entity.name for entity in entities]

    logger.info("VALIDATING {filenames}".format(filenames=", ".join(filenames)))

    file_users = [entities[0].modifiedBy, entities[0].createdBy]

    check_file_status = check_existing_file_status(
        validation_status_table, error_tracker_table, entities
    )

    status_list = check_file_status["status_list"]
    error_list = check_file_status["error_list"]

    messages_to_send = []
    # Need to figure out to how to remove this
    # This must pass in filenames, because filetype is determined by entity
    # name not by actual path of file
    validator = validate.GenieValidationHelper(
        syn=syn,
        project_id=project_id,
        center=center,
        entitylist=entities,
        format_registry=format_registry,
        genie_config=genie_config,
        ancillary_files=ancillary_files,
    )
    filetype = validator.file_type
    if check_file_status["to_validate"]:
        # HACK: Don't download again if only_validate is not True, but all
        # files need to be downloaded currently when validation + processing
        # isn't split up
        # if entities[0].get("path") is None:
        #    validator.entitylist = [syn.get(entity) for entity in entities]

        valid_cls, message = validator.validate_single_file(
            oncotree_link=genie_config["oncotreeLink"], nosymbol_check=False
        )

        logger.info("VALIDATION COMPLETE")
        input_status_list, invalid_errors_list = _get_status_and_error_list(
            valid_cls.is_valid(), message, entities
        )
        # Send email the first time the file is invalid
        if invalid_errors_list:
            messages_to_send.append((filenames, message, file_users))
    else:
        input_status_list = [
            {"entity": entity, "status": status}
            for entity, status in zip(entities, status_list)
        ]
        invalid_errors_list = [
            {"entity": entity, "errors": errors}
            for entity, errors in zip(entities, error_list)
        ]
    # add in static filetype and center information
    for input_status in input_status_list:
        input_status.update({"fileType": filetype, "center": center})
    # An empty list is returned if there are no errors,
    # so nothing will be appended
    for invalid_errors in invalid_errors_list:
        invalid_errors.update({"fileType": filetype, "center": center})
    return input_status_list, invalid_errors_list, messages_to_send


# TODO: Create ProcessHelper class
# TODO: Add to transform.py
def processfiles(
    syn,
    validfiles,
    center,
    path_to_genie,
    processing="main",
    format_registry=None,
    genie_config=None,
):
    """Processing validated files

    Args:
        syn: Synapse object
        validfiles: pandas dataframe containing validated files
                    has 'id', 'path', and 'fileType' column
        center: GENIE center name
        path_to_genie: Path to GENIE workdir
        center_mapping_df: Center mapping dataframe
        oncotree_link: Link to oncotree
        databaseToSynIdMappingDf: Database to synapse id mapping dataframe
        processing: Processing type. Defaults to main
    """
    logger.info(f"PROCESSING {center} FILES: {len(validfiles)}")
    center_staging_folder = os.path.join(path_to_genie, center)
    center_staging_synid = genie_config["center_config"][center]["stagingSynId"]

    if not os.path.exists(center_staging_folder):
        os.makedirs(center_staging_folder)

    if processing == "main":
        for _, row in validfiles.iterrows():
            filetype = row["fileType"]
            # filename = os.path.basename(filePath)
            newpath = os.path.join(center_staging_folder, row["name"])
            # store = True
            # Table id can be None
            tableid = genie_config.get(filetype)

            if filetype is not None and filetype != "other":
                # Example GENIE config can be found in tests/conftest.py
                processor = format_registry[filetype](
                    syn=syn, center=center, genie_config=genie_config
                )
                processor.process(
                    filePath=row["path"],
                    newPath=newpath,
                    parentId=center_staging_synid,
                    databaseSynId=tableid,
                    fileSynId=row["id"],
                )
    else:
        process_mutation.process_mutation_workflow(
            syn=syn,
            center=center,
            validfiles=validfiles,
            genie_config=genie_config,
            workdir=path_to_genie,
        )

    logger.info("ALL DATA STORED IN DATABASE")


# TODO: Add to validation.py
def append_duplication_errors(duplicated_filesdf, user_message_dict):
    """Duplicated files can occur because centers can upload files with the
    same filename in different folders.  This is to append duplication
    errors to the list of errors to email

    Args:
        duplicated_filesdf: Dataframe of duplciated files
        user_message_dict: Dictionary containing list of error messages to
                           send to each user.

    Returns:
        Dictionary containing list of error messages to send to each user.
    """
    duplication_error = (
        "Duplicated filename! Files should be uploaded as new versions "
        "and the entire dataset should be uploaded."
    )
    if not duplicated_filesdf.empty:
        filenames = []
        users = []
        for entity in duplicated_filesdf["entity"]:
            users.append(entity.modifiedBy)
            users.append(entity.createdBy)
            filenames.append(entity.name)
        file_messages = dict(filenames=filenames, messages=duplication_error)
        # Must get unique set of users or there
        # will be duplicated error messages sent in the email
        for user in set(users):
            user_message_dict[user].append(file_messages)
    return user_message_dict


# TODO: Add to validation.py
def get_duplicated_files(validation_statusdf):
    """
    Check for duplicated files.  There should be no duplication,
    files should be uploaded as new versions and the entire dataset
    should be uploaded everytime

    #TODO: This is a custom GENIE function

    Args:
        validation_statusdf: dataframe with 'name' and 'id' column
        duplicated_error_message: Error message for duplicated files

    Returns:
        dataframe with 'id', 'name', 'errors', 'center', 'fileType'
        and 'entity' of duplicated files
    """
    # This is special
    logger.info("CHECK FOR DUPLICATED FILES")
    duplicated_filesdf = validation_statusdf[
        validation_statusdf["name"].duplicated(keep=False)
    ]
    # Define filename str vector
    filename_str = validation_statusdf.name.str
    # cbs/seg files should not be duplicated.
    cbs_seg_index = filename_str.endswith(("cbs", "seg"))
    cbs_seg_files = validation_statusdf[cbs_seg_index]
    if len(cbs_seg_files) > 1:
        duplicated_filesdf = pd.concat([duplicated_filesdf, cbs_seg_files])
    # clinical files should not be duplicated.
    clinical_index = filename_str.startswith("data_clinical_supp")
    clinical_files = validation_statusdf[clinical_index]
    if len(clinical_files) > 2:
        duplicated_filesdf = pd.concat([duplicated_filesdf, clinical_files])
    duplicated_filesdf.drop_duplicates("id", inplace=True)
    logger.info("THERE ARE {} DUPLICATED FILES".format(len(duplicated_filesdf)))
    duplicated_filesdf["errors"] = DUPLICATED_FILE_ERROR
    return duplicated_filesdf


# TODO: Add to validation.py
def build_validation_status_table(input_valid_statuses: List[dict]):
    """Build validation status dataframe

    Args:
        input_valid_statuses: list of file validation status

    Returns:
        Validation status dataframe

    """
    status_table_columns = [
        "id",
        "path",
        "md5",
        "status",
        "name",
        "modifiedOn",
        "fileType",
        "center",
        "version",
        "entity",
    ]
    input_status_rows = []
    for input_status in input_valid_statuses:
        entity = input_status["entity"]
        row = {
            "id": entity.id,
            "path": entity.path,
            "md5": entity.md5,
            "status": input_status["status"],
            "name": entity.name,
            "modifiedOn": entity_date_to_timestamp(entity.properties.modifiedOn),
            "fileType": input_status["fileType"],
            "center": input_status["center"],
            "version": entity.versionNumber,
            "entity": entity,
        }
        input_status_rows.append(row)
    if input_status_rows:
        input_valid_statusdf = pd.DataFrame(input_status_rows)
    else:
        input_valid_statusdf = pd.DataFrame(
            input_status_rows, columns=status_table_columns
        )
    return input_valid_statusdf


# TODO: Add to validation.py
def build_error_tracking_table(invalid_errors: List[dict]):
    """Build error tracking dataframe

    Args:
        invalid_errors: list of file invalid errors

    Returns:
        Error tracking dataframe

    """
    error_table_columns = [
        "id",
        "errors",
        "name",
        "fileType",
        "center",
        "version",
        "entity",
    ]
    invalid_error_rows = []
    for invalid_error in invalid_errors:
        entity = invalid_error["entity"]
        row = {
            "id": entity.id,
            "errors": invalid_error["errors"],
            "name": entity.name,
            "fileType": invalid_error["fileType"],
            "center": invalid_error["center"],
            "version": entity.versionNumber,
            "entity": entity,
        }
        invalid_error_rows.append(row)
    if invalid_error_rows:
        invalid_errorsdf = pd.DataFrame(invalid_error_rows)
    else:
        invalid_errorsdf = pd.DataFrame(invalid_error_rows, columns=error_table_columns)
    return invalid_errorsdf


# TODO: Add to validation.py
def update_status_and_error_tables(
    syn,
    input_valid_statusdf,
    invalid_errorsdf,
    validation_status_table,
    error_tracker_table,
):
    """
    Update validation status and error tracking table

    Args:
        syn: Synapse object
        center: Center
        input_valid_status: list of lists of validation status
        invalid_errors: List of lists of invalid errors
        validation_status_table: Synapse table query of validation status
        error_tracker_table: Synapse table query of error tracker

    """
    logger.info("UPDATE VALIDATION STATUS DATABASE")

    load._update_table(
        syn,
        error_tracker_table.asDataFrame(),
        invalid_errorsdf,
        error_tracker_table.tableId,
        ["id"],
        to_delete=True,
    )

    load._update_table(
        syn,
        validation_status_table.asDataFrame(),
        input_valid_statusdf,
        validation_status_table.tableId,
        ["id"],
        to_delete=True,
    )


# TODO: Add to validation.py
def _update_tables_content(validation_statusdf, error_trackingdf):
    """Update validation status and error tracking dataframes with duplicated
    files.  Also update the error table to only contain errors - centers
    may have fixed their files so will want to remove old errors.

    Args:
        validation_statusdf: Validation status dataframe
        error_trackingdf: Error tracking dataframe

    Returns:
        dict: validation_statusdf: Updated validation status dataframe
              error_trackingdf: Updated error tracking dataframe
              duplicated_filesdf:  Duplicated files dataframe

    """
    # Get duplicated files
    duplicated_filesdf = get_duplicated_files(validation_statusdf)
    # index of all duplicated files
    duplicated_idx = validation_statusdf["id"].isin(duplicated_filesdf["id"])
    validation_statusdf["status"][duplicated_idx] = "INVALID"
    duplicated_idx = error_trackingdf["id"].isin(duplicated_filesdf["id"])
    error_trackingdf["errors"][duplicated_idx] = DUPLICATED_FILE_ERROR

    # Old errors are pulled down in validation, so obtain list of
    # files with duplicated file errors
    dup_ids = error_trackingdf["id"][
        error_trackingdf["errors"] == DUPLICATED_FILE_ERROR
    ]
    # Checks to see if the old duplicated files are still duplicated
    remove_ids = dup_ids[~dup_ids.isin(duplicated_filesdf["id"])]

    # Remove fixed duplicated files
    error_trackingdf = error_trackingdf[~error_trackingdf["id"].isin(remove_ids)]
    validation_statusdf = validation_statusdf[
        ~validation_statusdf["id"].isin(remove_ids)
    ]

    # Append duplicated file errors
    duplicated_filesdf["id"].isin(error_trackingdf["id"][duplicated_idx])
    error_trackingdf = pd.concat(
        [error_trackingdf, duplicated_filesdf[error_trackingdf.columns]]
    )
    # Remove duplicates if theres already an error that exists for the file
    error_trackingdf.drop_duplicates("id", inplace=True)

    # Since old errors are retained, make sure to only update
    # files that are actually invalid
    invalid_ids = validation_statusdf["id"][validation_statusdf["status"] == "INVALID"]
    error_trackingdf = error_trackingdf[error_trackingdf["id"].isin(invalid_ids)]
    # Fill blank file type values with 'other'
    error_trackingdf["fileType"].fillna("other", inplace=True)
    validation_statusdf["fileType"].fillna("other", inplace=True)

    return {
        "validation_statusdf": validation_statusdf,
        "error_trackingdf": error_trackingdf,
        "duplicated_filesdf": duplicated_filesdf,
    }


# TODO: Add to validation.py
def validation(
    syn,
    project_id,
    center,
    process,
    center_files,
    format_registry,
    genie_config,
    ancillary_files=None,
) -> pd.DataFrame:
    """
    Validation of all center files

    Args:
        syn: Synapse object
        center: Center name
        process: main, mutation
        center_files:
        format_registry:
        genie_config:
        ancillary_files:

    Returns:
        pd.DataFrame: Dataframe of valid GENIE files
    """
    logger.info(f"{center} has uploaded {len(center_files)} files.")
    validation_status_synid = genie_config["validationStatus"]
    error_tracker_synid = genie_config["errorTracker"]

    # Make sure the vcf validation statuses don't get wiped away
    # If process is not vcf, the vcf files are not downloaded
    # TODO: Add parameter to exclude types
    exclude_type = "vcf" if process != "mutation" else ""
    # id, md5, status, name, center, modifiedOn, fileType
    validation_status_table = syn.tableQuery(
        f"SELECT * FROM {validation_status_synid} where "
        f"center = '{center}' and fileType <> '{exclude_type}'"
    )
    # id, center, errors, name, fileType
    error_tracker_table = syn.tableQuery(
        f"SELECT * FROM {error_tracker_synid} where "
        f"center = '{center}' and fileType <> '{exclude_type}'"
    )

    input_valid_statuses = []
    invalid_errors = []

    # This default dict will capture all the error messages to send to
    # particular users
    user_message_dict = defaultdict(list)

    for ents in center_files:
        status, errors, messages_to_send = validatefile(
            syn=syn,
            project_id=project_id,
            entities=ents,
            validation_status_table=validation_status_table,
            error_tracker_table=error_tracker_table,
            center=center,
            format_registry=format_registry,
            genie_config=genie_config,
            ancillary_files=ancillary_files,
        )

        input_valid_statuses.extend(status)
        if errors is not None:
            invalid_errors.extend(errors)

        if messages_to_send:
            logger.debug("Collating messages to send to users.")
            for filenames, messages, users in messages_to_send:
                file_messages = dict(filenames=filenames, messages=messages)
                # Must get unique set of users or there
                # will be duplicated error messages sent in the email
                for user in set(users):
                    user_message_dict[user].append(file_messages)

    validation_statusdf = build_validation_status_table(input_valid_statuses)

    error_trackingdf = build_error_tracking_table(invalid_errors)

    new_tables = _update_tables_content(validation_statusdf, error_trackingdf)

    validation_statusdf = new_tables["validation_statusdf"]
    error_trackingdf = new_tables["error_trackingdf"]
    duplicated_filesdf = new_tables["duplicated_filesdf"]

    # In GENIE, we not only want to send out file format errors, but
    # also when there are duplicated errors.  The function below will
    # append duplication errors as an email to send to users (if applicable)
    user_message_dict = append_duplication_errors(duplicated_filesdf, user_message_dict)

    for user, message_objs in user_message_dict.items():
        logger.debug("Sending messages to user {user}.".format(user=user))
        # HACK: Avoid sending emails for staging project
        if project_id != "syn22033066":
            _send_validation_error_email(syn=syn, user=user, message_objs=message_objs)
    # \n write out new lines when they exist in the middle of a column
    # So the \n never gets uploaded into synapse table
    # change the delimiting to '|'.
    error_trackingdf["errors"] = [
        error.replace("\n", "|") for error in error_trackingdf["errors"]
    ]
    update_status_and_error_tables(
        syn=syn,
        input_valid_statusdf=validation_statusdf,
        invalid_errorsdf=error_trackingdf,
        validation_status_table=validation_status_table,
        error_tracker_table=error_tracker_table,
    )
    valid_filesdf = validation_statusdf.query('status == "VALIDATED"')
    return valid_filesdf[["id", "path", "fileType", "name"]]


# TODO: etl.py
def center_input_to_database(
    syn: Synapse,
    project_id: str,
    center: str,
    process: str,
    only_validate: bool,
    delete_old: bool = False,
    format_registry: Optional[dict] = None,
    genie_config: Optional[dict] = None,
):
    """Processing per center

    Args:
        syn (Synapse): Synapse connection
        project_id (str): GENIE Synapse project id
        center (str): GENIE center
        process (str): main or mutation processing
        only_validate (bool): Only validate or not
        delete_old (bool, optional): Delete old files. Defaults to False.
        format_registry (dict, optional): GENIE file format registry.
                                                 Defaults to None.
        genie_config (dict, optional): See example of genie config at
                                              ./genie_config.json. Defaults to None.
    """
    # TODO: Look into if errors should be thrown if these are None
    # Aka. should these actually be optional params
    if genie_config is None:
        genie_config = {}
    if format_registry is None:
        format_registry = {}

    if only_validate:
        log_path = os.path.join(
            process_functions.SCRIPT_DIR, f"{center}_validation_log.txt"
        )
    else:
        log_path = os.path.join(
            process_functions.SCRIPT_DIR, f"{center}_{process}_log.txt"
        )
    # Set up logger to write to a log file as well as streaming logs
    logFormatter = logging.Formatter(
        "%(asctime)s [%(name)s][%(levelname)s] %(message)s"
    )
    fileHandler = logging.FileHandler(log_path, mode="w")
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    # ----------------------------------------
    # Start processing
    # ----------------------------------------

    # path_to_genie = os.path.realpath(
    #   os.path.join(process_functions.SCRIPT_DIR, "../")
    # )
    # HACK:
    # Make the synapse cache dir the genie input folder for now
    # The main reason for this is because the .synaspeCache dir
    # is mounted by batch
    path_to_genie = os.path.expanduser("~/.synapseCache")
    # Create input and staging folders
    if not os.path.exists(os.path.join(path_to_genie, center, "input")):
        os.makedirs(os.path.join(path_to_genie, center, "input"))
    if not os.path.exists(os.path.join(path_to_genie, center, "staging")):
        os.makedirs(os.path.join(path_to_genie, center, "staging"))

    if delete_old:
        process_functions.rmFiles(os.path.join(path_to_genie, center))

    center_input_synid = genie_config["center_config"][center]["inputSynId"]
    logger.info("Center: " + center)
    center_files = extract.get_center_input_files(
        syn=syn,
        synid=center_input_synid,
        center=center,
        process=process,
        # HACK: Don't download all the files when only validate
        # downloadFile=(not only_validate),
    )

    # ancillary_files = get_ancillary_files(
    #    syn=syn,
    #    synid=center_input_synid,
    #    project_id=project_id,
    #    center=center,
    #    process=process,
    #    format_registry=format_registry,
    #    genie_config=genie_config,
    # )

    # only validate if there are center files
    if center_files:
        validFiles = validation(
            syn=syn,
            project_id=project_id,
            center=center,
            process=process,
            center_files=center_files,
            format_registry=format_registry,
            genie_config=genie_config,
            ancillary_files=center_files,
        )
    else:
        logger.info(f"{center} has not uploaded any files")
        return

    if len(validFiles) > 0 and not only_validate:
        # Reorganize so BED file are always validated and processed first
        bed_files = validFiles["fileType"] == "bed"
        beds = validFiles[bed_files]
        validFiles = pd.concat([beds, validFiles])
        validFiles.drop_duplicates(inplace=True)
        # merge clinical files into one row
        clinical_ind = validFiles["fileType"] == "clinical"
        if clinical_ind.any():
            clinical_files = validFiles[clinical_ind].to_dict(orient="list")
            # The [] implies the values in the dict as a list
            merged_clinical = pd.DataFrame([clinical_files])
            merged_clinical["fileType"] = "clinical"
            merged_clinical["name"] = f"data_clinical_supp_{center}.txt"
            validFiles = pd.concat([validFiles[~clinical_ind], merged_clinical])

        processTrackerSynId = genie_config["processTracker"]
        # Add process tracker for time start
        processTrackerDf = extract.get_syntabledf(
            syn=syn,
            query_string=(
                f"SELECT timeStartProcessing FROM {processTrackerSynId} "
                f"where center = '{center}' and "
                f"processingType = '{process}'"
            ),
        )
        if processTrackerDf.empty:
            new_rows = [
                [
                    center,
                    str(int(time.time() * 1000)),
                    str(int(time.time() * 1000)),
                    process,
                ]
            ]

            syn.store(synapseclient.Table(processTrackerSynId, new_rows))
        else:
            processTrackerDf["timeStartProcessing"].iloc[0] = str(
                int(time.time() * 1000)
            )
            syn.store(synapseclient.Table(processTrackerSynId, processTrackerDf))

        # Start transformations
        processfiles(
            syn=syn,
            validfiles=validFiles,
            center=center,
            path_to_genie=path_to_genie,
            processing=process,
            format_registry=format_registry,
            genie_config=genie_config,
        )

        # Should add in this process end tracking before the deletion of samples
        processTrackerDf = extract.get_syntabledf(
            syn=syn,
            query_string=(
                f"SELECT timeEndProcessing FROM {processTrackerSynId} "
                f"where center = '{center}' and "
                f"processingType = '{process}'"
            ),
        )
        processTrackerDf["timeEndProcessing"].iloc[0] = str(int(time.time() * 1000))
        syn.store(synapseclient.Table(processTrackerSynId, processTrackerDf))

        logger.info("SAMPLE/PATIENT RETRACTION")
        toRetract.retract(syn, project_id=project_id)

    else:
        messageOut = (
            f"{center} does not have any valid files"
            if not only_validate
            else f"ONLY VALIDATION OCCURED FOR {center}"
        )
        logger.info(messageOut)

    # Store and remove log file
    load.store_file(syn=syn, filepath=log_path, parentid=genie_config["logs"])
    os.remove(log_path)
    logger.info("ALL PROCESSES COMPLETE")
