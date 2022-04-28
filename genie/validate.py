#!/usr/bin/env python3
import logging

import synapseclient
from synapseclient.core.exceptions import SynapseHTTPError

from . import config, example_filetype_format, process_functions

logger = logging.getLogger(__name__)


class ValidationHelper(object):

    # Used for the kwargs in validate_single_file
    # Overload this per class
    _validate_kwargs = []

    def __init__(
        self,
        syn,
        project_id,
        center,
        entitylist,
        format_registry=None,
        file_type=None,
        genie_config=None,
    ):
        """A validator helper class for a center's files.

        Args:
            syn: a synapseclient.Synapse object
            project_id: Synapse Project ID where files are stored and configured.
            center: The participating center name.
            filepathlist: a list of file paths.
            format_registry: A dictionary mapping file format name to the
                             format class.
            file_type: Specify file type to skip filename validation
        """
        self._synapse_client = syn
        self._project = syn.get(project_id)
        self.entitylist = entitylist
        self.center = center
        self._format_registry = format_registry
        self.file_type = self.determine_filetype() if file_type is None else file_type
        self.genie_config = genie_config

    def determine_filetype(self):
        """Gets the file type of the file by validating its filename

        Args:
            syn: Synapse object
            filepathlist: list of filepaths to center files

        Returns:
            str: File type of input files.  None if no filetype found

        """
        filetype = None
        # Loop through file formats
        for file_format in self._format_registry:
            validator = self._format_registry[file_format](
                self._synapse_client, self.center
            )
            try:
                filenames = [entity.name for entity in self.entitylist]
                filetype = validator.validateFilename(filenames)
            except AssertionError:
                continue
            # If valid filename, return file type.
            if filetype is not None:
                break
        return filetype

    def validate_single_file(self, **kwargs):
        """Validate a submitted file unit.

        Returns:
            message: errors and warnings
            valid: Boolean value of validation status
        """

        if self.file_type not in self._format_registry:
            valid_result_cls = example_filetype_format.ValidationResults(
                errors="Your filename is incorrect! Please change your filename before you run the validator or specify --filetype if you are running the validator locally",
                warnings="",
            )
        else:
            mykwargs = {}
            for required_parameter in self._validate_kwargs:
                assert required_parameter in kwargs.keys(), (
                    "%s not in parameter list" % required_parameter
                )
                mykwargs[required_parameter] = kwargs[required_parameter]
                mykwargs["project_id"] = self._project.id

            validator_cls = self._format_registry[self.file_type]
            validator = validator_cls(
                syn=self._synapse_client,
                center=self.center,
                genie_config=self.genie_config,
            )
            filepathlist = [entity.path for entity in self.entitylist]
            valid_result_cls = validator.validate(filePathList=filepathlist, **mykwargs)

        # Complete error message
        message = valid_result_cls.collect_errors_and_warnings()
        return (valid_result_cls, message)


# TODO: Remove this at some point
class GenieValidationHelper(ValidationHelper):
    """A validator helper class for AACR Project Genie."""

    _validate_kwargs = ["nosymbol_check"]


def get_config(syn, synid):
    """Gets Synapse database to Table mapping in dict

    Args:
        syn: Synapse connection
        synid: Synapse id of database mapping table

    Returns:
        dict: {'databasename': 'synid'}

    """
    config = syn.tableQuery("SELECT * FROM {}".format(synid))
    configdf = config.asDataFrame()
    configdf.index = configdf["Database"]
    config_dict = configdf.to_dict()
    return config_dict["Id"]


def _check_parentid_permission_container(syn, parentid):
    """Checks permission / container
    # TODO: Currently only checks if a user has READ permissions
    """
    if parentid is not None:
        try:
            syn_ent = syn.get(parentid, downloadFile=False)
            # If not container, throw an assertion
            assert synapseclient.entity.is_container(syn_ent)
        except (SynapseHTTPError, AssertionError):
            raise ValueError(
                "Provided Synapse id must be your input folder Synapse id "
                "or a Synapse Id of a folder inside your input directory"
            )


def _check_center_input(center, center_list):
    """Checks center input

    Args:
        center: Center name
        center_list: List of allowed centers

    Raises:
        ValueError: If specify a center not part of the center list

    """
    if center not in center_list:
        raise ValueError(
            "Must specify one of these " f"centers: {', '.join(center_list)}"
        )


def _upload_to_synapse(syn, filepaths, valid, parentid=None):
    """
    Upload to synapse if parentid is specified and valid

    Args:
        syn: Synapse object
        filepaths: List of file paths
        valid: Boolean value for validity of file
        parentid: Synapse id of container. Default is None

    """
    if parentid is not None and valid:
        logger.info("Uploading file to {}".format(parentid))
        for path in filepaths:
            file_ent = synapseclient.File(path, parent=parentid)
            ent = syn.store(file_ent)
            logger.info("Stored to {}".format(ent.id))


def _perform_validate(syn, args):
    """This is the main entry point to the genie command line tool."""

    # Check parentid argparse
    _check_parentid_permission_container(syn=syn, parentid=args.parentid)
    genie_config = process_functions.get_genie_config(
        syn=syn, project_id=args.project_id
    )
    # HACK: Modify oncotree link config
    genie_config["oncotreeLink"] = process_functions._get_oncotreelink(
        syn=syn, genie_config=genie_config, oncotree_link=args.oncotree_link
    )
    # Check center argparse
    _check_center_input(args.center, list(genie_config["center_config"].keys()))

    format_registry = config.collect_format_types(args.format_registry_packages)
    logger.debug(f"Using {format_registry} file formats.")
    entity_list = [
        synapseclient.File(name=filepath, path=filepath, parentId=None)
        for filepath in args.filepath
    ]
    validator = GenieValidationHelper(
        syn=syn,
        project_id=args.project_id,
        center=args.center,
        entitylist=entity_list,
        format_registry=format_registry,
        file_type=args.filetype,
        genie_config=genie_config,
    )
    mykwargs = dict(
        nosymbol_check=args.nosymbol_check,
        project_id=args.project_id,
    )
    valid, message = validator.validate_single_file(**mykwargs)

    # Upload to synapse if parentid is specified and valid
    _upload_to_synapse(syn, args.filepath, valid, parentid=args.parentid)
