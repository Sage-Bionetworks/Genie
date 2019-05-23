#!/usr/bin/env python

from genie import PROCESS_FILES
import synapseclient
import logging
logger = logging.getLogger('genie')


def determine_filetype(syn, filepathlist, center):
    '''
    Get the file type of the file by validating its filename

    Args:
        syn: Synapse object
        filepathlist: list of filepaths to center files
        center: Participating Center

    Returns:
        str: File type of input files
    '''
    filetype = None
    # Loop through file formats
    for file_format in PROCESS_FILES:
        try:
            validator = PROCESS_FILES[file_format](syn, center)
            filetype = validator.validateFilename(filepathlist)
        except AssertionError:
            continue
        # If valid filename, return file type.
        if filetype is not None:
            break
    if filetype is None:
        raise ValueError(
                "Your filename is incorrect! "
                "Please change your filename before you run "
                "the validator or specify the --filetype.")
    return(filetype)


def validate(syn, filepathlist, center, filetype=None, threads=1,
             oncotreelink=None, testing=False, nosymbol_check=False):
    """
    This performs the validation of files

    :returns:   Text with the errors of the chosen file
    """
    if filetype is None:
        filetype = determine_filetype(syn, filepathlist, center)

    validator = PROCESS_FILES[filetype](syn, center, threads)
    total_error, warning = validator.validate(
        filePathList=filepathlist, oncotreeLink=oncotreelink,
        testing=testing, noSymbolCheck=nosymbol_check)

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

    return(message, valid)


# def main(syn, filepathlist, center, filetype, thread, oncotreelink,
#          parentid, testing, nosymbol_check):
#     if testing:
#         databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn11600968')
#     else:
#         databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn10967259')

#     databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()
#     synId = databaseToSynIdMappingDf.Id[
#         databaseToSynIdMappingDf['Database'] == "centerMapping"]
#     center_mapping = syn.tableQuery('SELECT * FROM %s' % synId[0])
#     center_mapping_df = center_mapping.asDataFrame()
#     assert center in center_mapping_df.center.tolist(), \
#         "Must specify one of these centers: {}".format(
#             ", ".join(center_mapping_df.center))

#     if oncotreelink is None:
#         oncoLink = databaseToSynIdMappingDf['Id'][
#             databaseToSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
#         oncoLinkEnt = syn.get(oncoLink)
#         oncotreelink = oncoLinkEnt.externalURL

#     message, valid = validate(
#         syn,
#         filepathlist,
#         center,
#         filetype,
#         thread,
#         oncotreelink,
#         testing,
#         nosymbol_check)

#     if parentid is not None:
#         if filetype is None:
#             raise ValueError(
#                 "If you specify the uploadToSynapse option, your filename "
#                 "must be named correctly")
#         else:
#             try:
#                 syn.get(parentid)
#             except synapseclient.exceptions.SynapseHTTPError:
#                 raise ValueError(
#                     "Provided Synapse id must be your input folder Synapse id "
#                     "or a Synapse Id of a folder inside your input directory")

#     if valid and parentid is not None:
#         logger.info("Uploading file to {}".format(parentid))
#         for path in filepathlist:
#             syn.store(synapseclient.File(path, parent=parentid))