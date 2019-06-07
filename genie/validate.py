#!/usr/bin/env python3

from genie import PROCESS_FILES
import synapseclient
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def determine_filetype(syn, filepathlist, center):
    '''
    Get the file type of the file by validating its filename

    Args:
        syn: Synapse object
        filepathlist: list of filepaths to center files
        center: Participating Center

    Returns:
        str: File type of input files.  None if no filetype found
    '''
    filetype = None
    # Loop through file formats
    for file_format in PROCESS_FILES:
        validator = PROCESS_FILES[file_format](syn, center)
        try:
            filetype = validator.validateFilename(filepathlist)
        except AssertionError:
            continue
        # If valid filename, return file type.
        if filetype is not None:
            break
    return(filetype)


def determine_validity_and_log(total_error, warning):
    '''
    Determines the validity of the file based on the
    the error message

    Args:
        total_error: string of file errors
        warning: string of file warnings

    Returns:
        valid - Boolean value of validation status
        message - error + warning
    '''
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
    return(valid, message)


def validate_single_file(syn,
                         filepathlist,
                         center,
                         filetype=None,
                         oncotreelink=None,
                         testing=False,
                         nosymbol_check=False):
    """
    This function determines the filetype of a file
    if filetype is not specified and logs the validation errors and
    warnings of a file.

    Args:
        syn: Synapse object
        filepathlist: List of files (only a list because of clinical file)
        center: Center name
        filetype: By default, filetype is determined from the filename.
                  filetype can be specified to avoid filename validation
        oncotreelink: Oncotree URL.
        testing:  Specify to invoke testing environment
        nosymbol_check: Do not check hugo symbols of fusion and cna file.
                        Default is False.

    Returns:
        message - errors and warnings
        valid - Boolean value of validation status
    """
    if filetype is None:
        filetype = determine_filetype(syn, filepathlist, center)
    if filetype not in PROCESS_FILES:
        raise ValueError(
            "Your filename is incorrect! "
            "Please change your filename before you run "
            "the validator or specify --filetype if you are "
            "running the validator locally")

    validator = PROCESS_FILES[filetype](syn, center)
    total_error, warning = validator.validate(
        filePathList=filepathlist, oncotreeLink=oncotreelink,
        testing=testing, noSymbolCheck=nosymbol_check)

    # Complete error message
    valid, message = determine_validity_and_log(total_error, warning)

    return(valid, message, filetype)


def perform_validate(syn, args):
    if args.testing:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn11600968')
    else:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn10967259')

    databasetosynid_mappingdf = databaseToSynIdMapping.asDataFrame()
    synid = databasetosynid_mappingdf.Id[
        databasetosynid_mappingdf['Database'] == "centerMapping"]
    center_mapping = syn.tableQuery('SELECT * FROM {}'.format(synid[0]))
    center_mapping_df = center_mapping.asDataFrame()
    assert args.center in center_mapping_df.center.tolist(), \
        "Must specify one of these centers: {}".format(
            ", ".join(center_mapping_df.center))

    if args.oncotreelink is None:
        oncolink = databasetosynid_mappingdf['Id'][
            databasetosynid_mappingdf['Database'] == 'oncotreeLink'].values[0]
        oncolink_ent = syn.get(oncolink)
        args.oncotreelink = oncolink_ent.externalURL

    valid, message, filetype = validate_single_file(
        syn, args.filepath, args.center, args.filetype,
        args.oncotreelink, args.testing, args.nosymbol_check)

    if args.parentid is not None:
        if args.filetype is not None:
            raise ValueError(
                "If you used --parentid, you must not use "
                "--filetype")
        if valid:
            try:
                syn.get(args.parentid)
            except synapseclient.exceptions.SynapseHTTPError:
                raise ValueError(
                    "Provided Synapse id must be your input folder Synapse id "
                    "or a Synapse Id of a folder inside your input directory")
            logger.info("Uploading file to {}".format(args.parentid))
            for path in args.filepath:
                syn.store(synapseclient.File(path, parent=args.parentid))
