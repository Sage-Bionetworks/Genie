#!/usr/bin/env python3

import logging

import synapseclient
from synapseclient.exceptions import SynapseHTTPError

from .config import PROCESS_FILES
from . import process_functions

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class ValidationHelper(object):

    # Used for the kwargs in validate_single_file
    # Overload this per class
    _validate_kwargs = []

    def __init__(self, syn, center, filepathlist,
                 format_registry=PROCESS_FILES,
                 testing=False):
        """A validator helper class for a center's files.

        Args:
            syn: a synapseclient.Synapse object
            center: The participating center name.
            filepathlist: a list of file paths.
            format_registry: A dictionary mapping file format name to the format class.
            testing: Run in testing mode.
        """

        self._synapse_client = syn
        self.filepathlist = filepathlist
        self.center = center
        self._format_registry = format_registry
        self.file_type = self.determine_filetype()
        self.testing = testing

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
            validator = self._format_registry[file_format](self._synapse_client, self.center,
                                                           testing=self.testing)
            try:
                filetype = validator.validateFilename(self.filepathlist)
            except AssertionError:
                continue
            # If valid filename, return file type.
            if filetype is not None:
                break
        return(filetype)

    def validate_single_file(self, **kwargs):
        """Validate a submitted file unit.

        Returns:
            message: errors and warnings
            valid: Boolean value of validation status
            filetype: String of the type of the file
        """

        if self.file_type not in self._format_registry:
            valid = False
            errors = "Your filename is incorrect! Please change your filename before you run the validator or specify --filetype if you are running the validator locally"
            warnings = ""
        else:
            mykwargs = {}
            for required_parameter in self._validate_kwargs:
                assert required_parameter in kwargs.keys(), \
                    "%s not in parameter list" % required_parameter
                mykwargs[required_parameter] = kwargs[required_parameter]

            validator_cls = self._format_registry[self.file_type]
            validator = validator_cls(self._synapse_client, self.center,
                                      testing=self.testing)
            valid, errors, warnings = validator.validate(filePathList=self.filepathlist,
                                                         **mykwargs)

        # Complete error message
        message = collect_errors_and_warnings(errors, warnings)

        return(valid, message, self.file_type)


class GenieValidationHelper(ValidationHelper):
    """A validator helper class for AACR Project Genie.
    """

    _validate_kwargs = ['oncotreeLink', 'noSymbolCheck']

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

    validator = GenieValidationHelper(syn=syn, center=args.center, 
                                      filepathlist=args.filepath)
    mykwargs = dict(oncotreeLink=args.oncotreelink, noSymbolCheck=args.nosymbol_check)
    valid, message, filetype = validator.validate_single_file(**mykwargs)

    # Upload to synapse if parentid is specified and valid
    _upload_to_synapse(syn, args.filepath, valid, parentid=args.parentid)
