#!/usr/bin/env python3
# noqa pylint: disable=line-too-long
"""genie cli"""
import argparse
import logging

from . import create_case_lists, process_functions, validate, write_invalid_reasons
from . import __version__


logger = logging.getLogger(__name__)


def perform_create_case_list(syn, args):
    """CLI to create case lists given a clinical and assay information file"""
    create_case_lists.main(
        args.clinical_file_name,
        args.assay_info_file_name,
        args.output_dir,
        args.study_id,
    )


def perform_get_file_errors(syn, args):
    """CLI to get invalid reasons"""
    project = syn.get(args.project_id)
    db_mapping = syn.tableQuery(f"select * from {project.dbMapping[0]}")
    db_mappingdf = db_mapping.asDataFrame()
    error_tracker_synid = db_mappingdf["Id"][
        db_mappingdf["Database"] == "errorTracker"
    ][0]
    center_errors = write_invalid_reasons.get_center_invalid_errors(
        syn, error_tracker_synid
    )
    print(center_errors[args.center])


def build_parser():
    parser = argparse.ArgumentParser(description="GENIE processing")

    parser.add_argument(
        "-v", "--version", action="version", version=f"genie {__version__}"
    )

    subparsers = parser.add_subparsers(
        title="commands",
        description="The following commands are available: ",
        help='For additional help use: "genie <COMMAND> -h"',
    )

    parser_validate = subparsers.add_parser(
        "validate", help="Validates GENIE file formats. "
    )

    parser_validate.add_argument(
        "filepath",
        type=str,
        nargs="+",
        help="File(s) that you are validating. "
        "If you have separate clinical sample and patient files, "
        "you must provide both files when validating.",
    )

    parser_validate.add_argument("center", type=str, help="Contributing Centers")

    validate_group = parser_validate.add_mutually_exclusive_group()

    validate_group.add_argument(
        "--filetype",
        type=str,
        help="Use the --filetype {FILETYPE} parameter to ignore filename validation. "
        "By default, the validator uses the filename to match "
        "the file format.  If your filename is incorrectly named, "
        "it will be invalid. "
        "Options: [maf, vcf, clinical, assayinfo, bed, cna, sv, seg, mutationsInCis]",
    )

    validate_group.add_argument(
        "--parentid",
        type=str,
        default=None,
        help="Synapse id of center input folder. "
        "If specified, your valid files will be uploaded "
        "to this directory.",
    )

    parser_validate.add_argument(
        "--oncotree_link",
        type=str,
        help="Specify an oncotree url when validating your clinical "
        "file "
        "(e.g: https://oncotree.info/api/tumorTypes/tree?version=oncotree_2021_11_02). "
        "By default the oncotree version used will be specified in this entity: "
        "syn13890902",
    )

    parser_validate.add_argument(
        "--nosymbol-check",
        action="store_true",
        help="Ignores specific post-processing validation criteria related to HUGO symbols "
        "in the structural variant and cna files.",
    )

    # TODO: remove this default when private genie project is ready
    parser_validate.add_argument(
        "--project_id",
        type=str,
        default="syn3380222",
        help="FOR DEVELOPER USE ONLY: Synapse Project ID where data is stored. "
        "(default: %(default)s).",
    )

    parser_validate.add_argument(
        "--format_registry_packages",
        type=str,
        nargs="+",
        default=["genie_registry"],
        help="FOR DEVELOPER USE ONLY: Python package name(s) to get valid file formats "
        "from (default: %(default)s).",
    )

    parser_validate.set_defaults(func=validate._perform_validate)

    parser_create_case = subparsers.add_parser(
        "create-case-lists", help="Creates cBioPortal case list files"
    )
    parser_create_case.add_argument(
        "clinical_file_name", type=str, help="Clinical file path"
    )
    parser_create_case.add_argument(
        "assay_info_file_name", type=str, help="gene matrix file path"
    )
    parser_create_case.add_argument("output_dir", type=str, help="Output directory")
    parser_create_case.add_argument("study_id", type=str, help="Output directory")
    parser_create_case.set_defaults(func=perform_create_case_list)

    parser_get_invalid = subparsers.add_parser(
        "get-file-errors", help="Get the file invalid reasons for a specific center"
    )
    parser_get_invalid.add_argument("center", type=str, help="Contributing Centers")
    parser_get_invalid.add_argument(
        "--project_id",
        type=str,
        default="syn3380222",
        help="Synapse Project ID where data is stored. (default: %(default)s).",
    )
    parser_get_invalid.set_defaults(func=perform_get_file_errors)

    return parser


def main():
    """Invoke"""
    args = build_parser().parse_args()
    syn = process_functions.synapse_login()
    # func has to match the set_defaults
    args.func(syn, args)


if __name__ == "__main__":
    main()
