"""
Validation helper functions
"""
import datetime
import logging
from typing import List

import pandas as pd
from py import process
import synapseclient
from synapseclient.core.exceptions import SynapseHTTPError

from . import config, process_functions

logger = logging.getLogger(__name__)


class ValidationHelper(object):

    # Used for the kwargs in validate_single_file
    # Overload this per class
    _validate_kwargs = []

    def __init__(
        self, syn, project_id, center, entitylist, format_registry=None, file_type=None
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
            valid = False
            errors = "Your filename is incorrect! Please change your filename before you run the validator or specify --filetype if you are running the validator locally"
            warnings = ""
        else:
            mykwargs = {}
            for required_parameter in self._validate_kwargs:
                assert required_parameter in kwargs.keys(), (
                    "%s not in parameter list" % required_parameter
                )
                mykwargs[required_parameter] = kwargs[required_parameter]
                mykwargs["project_id"] = self._project.id

            validator_cls = self._format_registry[self.file_type]
            validator = validator_cls(self._synapse_client, self.center)
            filepathlist = [entity.path for entity in self.entitylist]
            valid, errors, warnings = validator.validate(
                filePathList=filepathlist, **mykwargs
            )

        # Complete error message
        message = collect_errors_and_warnings(errors, warnings)

        return (valid, message)


# TODO: Remove this at some point
class GenieValidationHelper(ValidationHelper):
    """A validator helper class for AACR Project Genie."""

    _validate_kwargs = ["oncotree_link", "nosymbol_check"]


def collect_errors_and_warnings(errors, warnings):
    """Aggregates error and warnings into a string.

    Args:
        errors: string of file errors, separated by new lines.
        warnings: string of file warnings, separated by new lines.

    Returns:
        message - errors + warnings
    """
    # Complete error message
    message = "----------------ERRORS----------------\n"
    if errors == "":
        message = "YOUR FILE IS VALIDATED!\n"
        logger.info(message)
    else:
        for error in errors.split("\n"):
            if error != "":
                logger.error(error)
        message += errors
    if warnings != "":
        for warning in warnings.split("\n"):
            if warning != "":
                logger.warning(warning)
        message += "-------------WARNINGS-------------\n" + warnings
    return message


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


def _get_oncotreelink(syn, databasetosynid_mappingdf, oncotree_link=None):
    """
    Gets oncotree link unless a link is specified by the user

    Args:
        syn: Synapse object
        databasetosynid_mappingdf: database to synid mapping
        oncotree_link: link to oncotree. Default is None

    Returns:
        oncotree link

    """
    if oncotree_link is None:
        oncolink = databasetosynid_mappingdf.query('Database == "oncotreeLink"').Id
        oncolink_ent = syn.get(oncolink.iloc[0])
        oncotree_link = oncolink_ent.externalURL
    return oncotree_link


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
    _check_parentid_permission_container(syn, args.parentid)

    databasetosynid_mappingdf = process_functions.get_synid_database_mappingdf(
        syn, project_id=args.project_id
    )

    synid = databasetosynid_mappingdf.query('Database == "centerMapping"').Id

    center_mapping = syn.tableQuery("select * from {}".format(synid.iloc[0]))
    center_mapping_df = center_mapping.asDataFrame()

    # Check center argparse
    _check_center_input(args.center, center_mapping_df.center.tolist())

    args.oncotree_link = _get_oncotreelink(
        syn, databasetosynid_mappingdf, oncotree_link=args.oncotree_link
    )

    format_registry = config.collect_format_types(args.format_registry_packages)
    logger.debug("Using {} file formats.".format(format_registry))
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
    )
    mykwargs = dict(
        oncotree_link=args.oncotree_link,
        nosymbol_check=args.nosymbol_check,
        project_id=args.project_id,
    )
    valid, message = validator.validate_single_file(**mykwargs)

    # Upload to synapse if parentid is specified and valid
    _upload_to_synapse(syn, args.filepath, valid, parentid=args.parentid)


# def _check_year(year: int) -> bool:
#     """_summary_

#     Args:
#         year (int): _description_

#     Returns:
#         bool: _description_
#     """
#     try:
#         year = datetime.datetime.strptime(str(int(year)), "%Y").year
#         return year
#     except Exception:
#         return False


def check_year(
    df: pd.DataFrame,
    col: int,
    center: str = None,
    allowed_string_values: list = None,
) -> List[dict]:
    """Check year values

    Args:
        df (pd.DataFrame): Clinical dataframe
        col (int): Column name
        center (str, optional): Center name. Defaults to None.
        allowed_string_values (list, optional): list of other allowed string values.
                                                Defaults to None.

    Returns:
        List[dict]: (row index, summary, detailed, error or warning)
    """
    year_now = datetime.datetime.utcnow().year
    # Generate summary error
    summary = (
        f"Please double check your {col} "
        "column, it must be an integer in YYYY format "
        f"<= {year_now}"
    )
    # Tack on allowed string values
    if allowed_string_values:
        summary += " or '{}'.\n".format("', '".join(allowed_string_values))
    else:
        summary += ".\n"

    row_errors = []
    if allowed_string_values is None:
        allowed_string_values = []
    if process_functions.checkColExist(df, col):
        # Deal with pre-redacted values and other allowed strings
        # first because can't int(text) because there are
        # instances that have <YYYY
        year_series = df[col][~df[col].isin(allowed_string_values)]
        for index, year in year_series.iteritems():
            try:
                # If year is greater than current year, it is invalid
                invalid_year = datetime.datetime.strptime(
                    str(int(year)), "%Y"
                ).year > year_now
                # Make sure that none of the years are greater than the current
                # year.  It can be the same, but can't future years.
            except Exception:
                invalid_year = True
            if invalid_year:
                detailed = f"{year} is not a valid year."
                row_errors.append(
                    {
                        "index": index,
                        "summary": summary,
                        "detailed": detailed,
                        "check_level": "error"
                    }
                )
    return row_errors


def check_required_columns(
    df: pd.DataFrame,
    cols: List[int],
    center: str = None,
    allowed_string_values: list = None
) -> List[dict]:
    """_summary_

    Args:
        df (pd.DataFrame): _description_
        cols List[int]: List of required columns
        center (str, optional): _description_. Defaults to None.
        allowed_string_values (list, optional): _description_. Defaults to None.

    Returns:
        List[dict]: _description_
    """
    row_errors = []
    for col in cols:
        if not process_functions.checkColExist(df, col):
            row_errors.append(
                {
                    "index": None,
                    "summary": f"Must have {col} column.",
                    "detailed": None,
                    "check_level": "error"
                }
            )
    return row_errors


def check_empty_rows(
    df: pd.DataFrame,
    cols: List[int],
    center: str = None,
    allowed_string_values: list = None
) -> List[tuple]:
    row_errors = []
    empty_rows = df.isnull().sum(axis=1) == len(df.columns)
    for index, value in empty_rows[empty_rows].iteritems():
        # Remove completely empty rows to speed up processing
        row_errors.append(
            {
                "index": index,
                "summary": "No empty rows allowed.\n",
                "detailed": "No empty rows allowed.\n",
                "check_level": "error"
            }
        )
    return row_errors


def check_duplicated_values(
    df: pd.DataFrame,
    cols: List[int],
    center: str = None,
    allowed_string_values: list = None
) -> List[dict]:
    """_summary_

    Args:
        df (pd.DataFrame): _description_
        cols (List[int]): _description_
        center (str, optional): _description_. Defaults to None.
        allowed_string_values (list, optional): _description_. Defaults to None.

    Returns:
        List[dict]: _description_
    """
    row_errors = []
    for col in cols:
        is_duplicated = df[col][df[col].duplicated()]
        for index, value in is_duplicated.iteritems():
            row_errors.append(
                {
                    "index": index,
                    "summary": (
                        "No duplicated SAMPLE_ID "
                        "allowed.\nIf there are no duplicated "
                        "SAMPLE_IDs, and both sample and patient files are "
                        "uploaded, then please check to make sure no duplicated "
                        "PATIENT_IDs exist in the patient clinical file.\n"
                    ),
                    "detailed": f"{value} is duplicated!",
                    "check_level": "error"
                }
            )
    return row_errors
