#!/usr/bin/env python3

import logging

import synapseclient
from synapseclient.exceptions import SynapseHTTPError

from .config import PROCESS_FILES
from . import process_functions

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Validator(object):
    """A validator helper class for a center's files.
    """

    def __init__(self, syn, center, filepathlist, format_registry=PROCESS_FILES):
        """A validator helper class for a center's files.

        Args:
            center: Participating Center
        """

        self._synapse_client = syn
        self.filepathlist = filepathlist
        self.center = center
        self._format_registry = format_registry
        self.file_type = self.determine_filetype()

    def determine_filetype(self):
        '''
        Get the file type of the file by validating its filename

        Args:
            syn: Synapse object
            filepathlist: list of filepaths to center files

        Returns:
            str: File type of input files.  None if no filetype found
        '''
        filetype = None
        # Loop through file formats
        for file_format in self._format_registry:
            validator = self._format_registry[file_format](self._synapse_client, self.center)
            try:
                filetype = validator.validateFilename(self.filepathlist)
            except AssertionError:
                continue
            # If valid filename, return file type.
            if filetype is not None:
                break
        return(filetype)

    def validate_single_file(self, filetype=None,
                             oncotreelink=None, testing=False, 
                             nosymbol_check=False):
        """
        This function determines the filetype of a single submitted 'file'.
        The 'file' should be one of those defined in config.PROCESS_FILES and 
        may actually be composed of multiple files.
        if filetype is not specified and logs the validation errors and
        warnings of a file.

        Args:
            filepathlist: List of local paths to files.
            center: Center name
            filetype: If None, filetype is determined from the filename.
                    If specified, filetype determination is skipped.
            oncotreelink: Oncotree URL.
            testing: Specify to invoke testing environment
            nosymbol_check: Do not check hugo symbols of fusion and cna file.

        Returns:
            message: errors and warnings
            valid: Boolean value of validation status
            filetype: String of the type of the file
        """
        if filetype is None:
            filetype = self.file_type

        if filetype not in self._format_registry:
            valid = False
            errors = "Your filename is incorrect! Please change your filename before you run the validator or specify --filetype if you are running the validator locally"
            warnings = ""
        else:
            validator = self._format_registry[filetype](self._synapse_client, self.center)
            valid, errors, warnings = validator.validate(filePathList=self.filepathlist, 
                                                         oncotreeLink=oncotreelink,
                                                         testing=testing,
                                                         noSymbolCheck=nosymbol_check)

        # Complete error message
        message = collect_errors_and_warnings(errors, warnings)

        return(valid, message, filetype)


        

# Validates annotations on Synapse
# def validateAnnotations(fileList):
#     logger.info("VALIDATING ANNOTATIONS")
#     notcorrect = []
#     for i,ID in enumerate(fileList['entity.id']):
#         foo = syn.get(ID, downloadFile=False)
#         required_annot = ["center","dataType","fileType","disease","consortium",
#         "platform","tissueSource","organism","dataSubType"]
#         check = [annot for annot in required_annot if foo.annotations.has_key(annot)]
#         if len(check) != len(required_annot):
#             notcorrect.append(fileList.iloc[i]['entity.id'])
#     if len(notcorrect) >0:
#         return(False)
#     else:
#         return(True)


def collect_errors_and_warnings(errors, warnings):
    '''Aggregates error and warnings into a string.

    Args:
        errors: string of file errors, separated by new lines.
        warnings: string of file warnings, separated by new lines.

    Returns:
        message - errors + warnings
    '''
    # Complete error message
    message = "----------------ERRORS----------------\n"
    if errors == "":
        message = "YOUR FILE IS VALIDATED!\n"
        logger.info(message)
    else:
        for error in errors.split("\n"):
            if error != '':
                logger.error(error)
        message += errors
    if warnings != "":
        for warning in warnings.split("\n"):
            if warning != '':
                logger.warning(warning)
        message += "-------------WARNINGS-------------\n" + warnings
    return message


def get_config(syn, synid):
    '''
    Get Synapse database to Table mapping in dict
    '''
    config = syn.tableQuery('SELECT * FROM {}'.format(synid))
    configdf = config.asDataFrame()
    configdf.index = configdf['Database']
    config_dict = configdf.to_dict()
    return(config_dict['Id'])


def _check_parentid_permission_container(syn, parentid):
    '''
    Check permission / container
    Currently only checks if a user has READ permissions...
    '''
    if parentid is not None:
        try:
            syn_ent = syn.get(parentid, downloadFile=False)
            # If not container, throw an assertion
            assert synapseclient.entity.is_container(syn_ent)
        except (SynapseHTTPError, AssertionError):
            raise ValueError(
                "Provided Synapse id must be your input folder Synapse id "
                "or a Synapse Id of a folder inside your input directory")


def _check_center_input(center, center_list):
    '''
    Check center input

    Args:
        center: Center name
        center_list: List of allowed centers
    '''
    if center not in center_list:
        raise ValueError(
            "Must specify one of these centers: {}".format(
                ", ".join(center_list)))


def _get_oncotreelink(syn, databasetosynid_mappingdf, oncotreelink=None):
    '''
    Get oncotree link unless a link is specified by the user

    Args:
        syn: Synapse object
        databasetosynid_mappingdf: database to synid mapping
        oncotreelink: link to oncotree. Default is None
    '''
    if oncotreelink is None:
        oncolink = databasetosynid_mappingdf.query(
            'Database == "oncotreeLink"').Id
        oncolink_ent = syn.get(oncolink.iloc[0])
        oncotreelink = oncolink_ent.externalURL
    return(oncotreelink)


def _upload_to_synapse(syn, filepaths, valid, parentid=None):
    '''
    Upload to synapse if parentid is specified and valid

    Args:
        syn: Synapse object
        filepaths: List of file paths
        valid: Boolean value for validity of file
        parentid: Synapse id of container. Default is None
    '''
    if parentid is not None and valid:
        logger.info("Uploading file to {}".format(parentid))
        for path in filepaths:
            file_ent = synapseclient.File(path, parent=parentid)
            ent = syn.store(file_ent)
            logger.info("Stored to {}".format(ent.id))


def _perform_validate(syn, args):
    """This is the main entry point to the genie command line tool.
    """

    # Check parentid argparse
    _check_parentid_permission_container(syn, args.parentid)

    databasetosynid_mappingdf = process_functions.get_synid_database_mappingdf(
        syn, test=args.testing)

    synid = databasetosynid_mappingdf.query('Database == "centerMapping"').Id

    center_mapping = syn.tableQuery('select * from {}'.format(synid.iloc[0]))
    center_mapping_df = center_mapping.asDataFrame()

    # Check center argparse
    _check_center_input(args.center, center_mapping_df.center.tolist())

    args.oncotreelink = _get_oncotreelink(syn, databasetosynid_mappingdf,
                                          oncotreelink=args.oncotreelink)

    validator = Validator(syn=syn, center=args.center, filepathlist=args.filepath)
    valid, message, filetype = validator.validate_single_file(args.filetype,
        args.oncotreelink, args.testing, args.nosymbol_check)

    # Upload to synapse if parentid is specified and valid
    _upload_to_synapse(syn, args.filepath, valid, parentid=args.parentid)
