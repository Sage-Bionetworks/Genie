#!/usr/bin/env python3

from genie import PROCESS_FILES
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def determine_filetype(syn, filepathlist, center, raise_error=True):
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
        validator = PROCESS_FILES[file_format](syn, center)
        try:
            filetype = validator.validateFilename(filepathlist)
        except AssertionError:
            continue
        # If valid filename, return file type.
        if filetype is not None:
            break
    if filetype is None and raise_error:
        raise ValueError(
            "Your filename is incorrect! "
            "Please change your filename before you run "
            "the validator or specify --filetype if you are "
            "running the validator locally")
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


def validate_single_file_workflow(syn,
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
        filetype = \
            determine_filetype(syn, filepathlist, center, raise_error=True)
    else:
        if filetype not in PROCESS_FILES:
            raise ValueError("Must specify correct filetype")

    validator = PROCESS_FILES[filetype](syn, center)
    total_error, warning = validator.validate(
        filePathList=filepathlist, oncotreeLink=oncotreelink,
        testing=testing, noSymbolCheck=nosymbol_check)

    # Complete error message
    valid, message = determine_validity_and_log(total_error, warning)

    return(valid, message, filetype)
