#!/usr/bin/env python3
import logging

import pandas as pd
import synapseclient
from synapseclient.core.exceptions import SynapseHTTPError

from genie import (
    config,
    example_filetype_format,
    extract,
    load,
    process_functions,
    transform,
)

logger = logging.getLogger(__name__)


ACCEPTED_CHROMOSOMES = list(map(str, range(1, 23))) + ["X", "Y", "MT"]


# TODO: move to models.py
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


# TODO: Currently only checks if a user has READ permissions
def _check_parentid_permission_container(syn, parentid):
    """Checks permission / container"""
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


def _validate_chromosome(
    df: pd.DataFrame,
    col: str,
    fileformat: str,
    allow_chr: bool = True,
    allow_na: bool = False,
) -> tuple:
    """Validate chromosome values

    Args:
        df (pd.DataFrame): Dataframe
        col (str): Column header for column containing chromosome values
        fileformat (str): GENIE supported file format
        allow_chr (bool): whether the chr prefix is allowed in the values
        allow_na (bool): whether NA/blanks are allowed in the values

    Returns:
        tuple: errors and warnings
    """
    have_column = process_functions.checkColExist(df, col)
    errors = ""
    warnings = ""
    if have_column:
        nochr = ["chr" in i for i in df[col] if isinstance(i, str)]
        if sum(nochr) > 0:
            if allow_chr:
                warnings += f"{fileformat}: Should not have the chr prefix in front of chromosomes.\n"
            else:
                errors += f"{fileformat}: Should not have the chr prefix in front of chromosomes.\n"
        # correct_chromosomes = [
        #     str(chrom).replace("chr", "") in accepted_chromosomes
        #     for chrom in df[col]
        # ]
        # preserve NAs
        df[col] = transform._convert_float_col_with_nas_to_int(df, col)
        df[col] = transform._convert_col_with_nas_to_str(df, col)
        df[col] = [val.replace("chr", "") if pd.notna(val) else val for val in df[col]]
        warning, error = process_functions.check_col_and_values(
            df=df,
            col=col,
            possible_values=ACCEPTED_CHROMOSOMES,
            filename=fileformat,
            na_allowed=allow_na,
        )
        errors += error
        warnings += warning
    return (errors, warnings)


# TODO: specify all the arguments instead of using args.
# TODO: This is a cli function...
def _perform_validate(syn, args):
    """This is the main entry point to the genie command line tool."""

    # Check parentid argparse
    _check_parentid_permission_container(syn=syn, parentid=args.parentid)
    genie_config = extract.get_genie_config(syn=syn, project_id=args.project_id)
    # HACK: Modify oncotree link config
    # TODO: Remove oncotree_link parameter from this function
    genie_config["oncotreeLink"] = extract._get_oncotreelink(
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
    if valid and args.parentid is not None:
        logger.info(f"Uploading files to {args.parentid}")
        load.store_files(syn=syn, filepaths=args.filepath, parentid=args.parentid)
