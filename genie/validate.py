#!/usr/bin/env python

from genie import PROCESS_FILES
import synapseclient
import os
import getpass
import logging
logger = logging.getLogger('genie')


def synapse_login():
    """
    This function logs into synapse for you if credentials are saved.
    If not saved, then user is prompted username and password.

    :returns:     Synapseclient object
    """
    try:
        syn = synapseclient.login(silent=True)
    except Exception:
        logger.info(
            "Please provide your synapse username/email and password "
            "(You will only be prompted once)")
        Username = os.getenv("GENIE_USER")
        Password = os.getenv("GENIE_PASS")
        if Username is None or Password is None:
            Username = raw_input("Username: ")
            Password = getpass.getpass()
        syn = synapseclient.login(
            email=Username, password=Password, rememberMe=True, silent=True)
    return(syn)


def get_filetype(syn, path_list, center):
    '''
    Get the file type of the file by validating its filename

    Args:
        syn: Synapse object
        path_list: list of filepaths to center files
        center: Participating Center

    Returns:
        str: File type of input files
    '''
    filetype = None
    for file_format in PROCESS_FILES:
        try:
            filetype = PROCESS_FILES[file_format](
                syn, center).validateFilename(path_list)
        except AssertionError:
            continue
        # If valid filename, return file type.
        if filetype is not None:
            break
    return(filetype)


def validate(syn, filePath, center, threads, fileType=None,
             oncotree_url=None, uploadToSynapse=None,
             testing=False, noSymbolCheck=False):
    """
    This performs the validation of files

    :returns:   Text with the errors of the chosen file
    """
    # CHECK: Fail if filename is incorrect
    validator = PROCESS_FILES[fileType](syn, center, threads)
    if not offline:
        try:
            validator.validateFilename(filePath)
        except AssertionError as e:
            raise ValueError(
                "Your filename is incorrect!\n{}\n"
                "Please change your filename before you run "
                "the validator again.".format(e))
    total_error, warning = validator.validate(
        filePathList=filePath, oncotreeLink=oncotree_url,
        testing=testing, noSymbolCheck=noSymbolCheck)

    # Complete error message
    message = "----------------ERRORS----------------\n"
    if total_error == "":
        message = "YOUR FILE IS VALIDATED!\n"
        logger.info(message)
        valid = True
    else:
        for errors in total_error.split("\n"):
            if errors != '':
                logger.error(errors)
        message += total_error
        valid = False
    if warning != "":
        for warn in warning.split("\n"):
            if warn != '':
                logger.warning(warn)
        message += "-------------WARNINGS-------------\n" + warning
    if valid and uploadToSynapse is not None:
        logger.info("Uploading file to %s" % uploadToSynapse)
        for path in filePath:
            syn.store(synapseclient.File(path, parent=uploadToSynapse))
    return(message, valid)


def main(syn, filepath, center, thread, oncotreelink, filetype,
         upload_to_synapse, testing, nosymbol_check):
    if testing:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn11600968')
    else:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn10967259')

    databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()
    synId = databaseToSynIdMappingDf.Id[
        databaseToSynIdMappingDf['Database'] == "centerMapping"]
    center_mapping = syn.tableQuery('SELECT * FROM %s' % synId[0])
    center_mapping_df = center_mapping.asDataFrame()
    assert center in center_mapping_df.center.tolist(), \
        "Must specify one of these centers: {}".format(
            ", ".join(center_mapping_df.center))

    if oncotreelink is None:
        oncoLink = databaseToSynIdMappingDf['Id'][
            databaseToSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
        oncoLinkEnt = syn.get(oncoLink)
        oncotreelink = oncoLinkEnt.externalURL

    if upload_to_synapse is not None:
        if filetype is None:
            raise ValueError(
                "If you specify the uploadToSynapse option, your filename "
                "must be named correctly")
        else:
            try:
                syn.get(upload_to_synapse)
            except synapseclient.exceptions.SynapseHTTPError:
                raise ValueError(
                    "Provided Synapse id must be your input folder Synapse id "
                    "or a Synapse Id of a folder inside your input directory")

    validate(syn,
             filepath,
             center,
             thread,
             oncotreelink,
             filetype,
             upload_to_synapse,
             testing,
             nosymbol_check)
