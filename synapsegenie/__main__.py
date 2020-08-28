#!/usr/bin/env python3
# noqa pylint: disable=line-too-long
"""genie cli"""
import argparse
import logging

import synapseclient

from synapsegenie import (bootstrap, config, input_to_database,
                          process_functions, validate, write_invalid_reasons)

from .__version__ import __version__


logger = logging.getLogger('genie')


def synapse_login(username=None, password=None):
    """
    This function logs into synapse for you if credentials are saved.
    If not saved, then user is prompted username and password.

    :returns:     Synapseclient object
    """
    try:
        syn = synapseclient.login(silent=True)
    except Exception:
        if username is None and password is None:
            raise ValueError("Please specify --syn_user, --syn_pass to specify your Synapse "
                             "login. Please view https://docs.synapse.org/articles/client_configuration.html"
                             "to learn about logging into Synapse via the Python client.")
        syn = synapseclient.login(email=username,
                                  password=password,
                                  silent=True)
    return syn


def bootstrap_infra(syn, args):
    """Create GENIE-like infrastructure"""
    bootstrap.main(syn)


def process_cli_wrapper(syn, args):
    """Process CLI wrapper"""
    process(syn, args.process, args.project_id, center=args.center,
            pemfile=args.pemfile, delete_old=args.delete_old,
            only_validate=args.only_validate, debug=args.debug,
            format_registry_packages=args.format_registry_packages)


def process(syn, process, project_id, center=None, pemfile=None,
            delete_old=False, only_validate=False, debug=False,
            format_registry_packages=None):
    """Process files"""
    # Get the Synapse Project where data is stored
    # Should have annotations to find the table lookup
    project = syn.get(project_id)
    database_to_synid_mapping_synid = project.annotations.get("dbMapping", "")

    databaseToSynIdMapping = syn.tableQuery(
        'SELECT * FROM {}'.format(database_to_synid_mapping_synid[0]))
    databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()

    center_mapping_id = process_functions.getDatabaseSynId(
        syn, "centerMapping",
        databaseToSynIdMappingDf=databaseToSynIdMappingDf
    )

    center_mapping = syn.tableQuery(f'SELECT * FROM {center_mapping_id}')
    center_mapping_df = center_mapping.asDataFrame()

    if center is not None:
        assert center in center_mapping_df.center.tolist(), (
            "Must specify one of these centers: {}".format(
                ", ".join(center_mapping_df.center)))
        centers = [center]
    else:
        center_mapping_df = center_mapping_df[
            ~center_mapping_df['inputSynId'].isnull()
        ]
        # release is a bool column
        center_mapping_df = center_mapping_df[center_mapping_df['release']]
        centers = center_mapping_df.center

    format_registry = config.collect_format_types(format_registry_packages)

    for process_center in centers:
        input_to_database.center_input_to_database(
            syn, project_id, process_center, process,
            only_validate, databaseToSynIdMappingDf,
            center_mapping_df,
            delete_old=delete_old,
            format_registry=format_registry
        )

    error_tracker_synid = process_functions.getDatabaseSynId(
        syn, "errorTracker", databaseToSynIdMappingDf=databaseToSynIdMappingDf
    )
    # Only write out invalid reasons if the center
    # isnt specified and if only validate
    if center is None and only_validate:
        logger.info("WRITING INVALID REASONS TO CENTER STAGING DIRS")
        write_invalid_reasons.write_invalid_reasons(
            syn, center_mapping_df, error_tracker_synid
        )

def build_parser():
    """Build CLI parsers"""
    parser = argparse.ArgumentParser(description='GENIE processing')

    parser.add_argument("--syn_user", type=str, help='Synapse username')

    parser.add_argument("--syn_pass", type=str, help='Synapse password')

    parser.add_argument('-v', '--version', action='version',
                        version='genie {}'.format(__version__))

    subparsers = parser.add_subparsers(title='commands',
                                       description='The following commands are available:',
                                       help='For additional help: "genie <COMMAND> -h"')

    parser_validate = subparsers.add_parser('validate', help='Validates GENIE file formats')

    parser_validate.add_argument("filepath", type=str, nargs="+",
                                 help='File(s) that you are validating.'
                                      'If you validation your clinical files and you have both sample and '
                                      'patient files, you must provide both')

    parser_validate.add_argument("center", type=str, help='Contributing Centers')

    parser_validate.add_argument("--format_registry_packages", type=str, nargs="+",
                                 default=["genie"],
                                 help="Python package name(s) to get valid file formats from (default: %(default)s).")

    parser_validate.add_argument("--oncotree_link", type=str, help="Link to oncotree code")

    validate_group = parser_validate.add_mutually_exclusive_group()

    validate_group.add_argument("--filetype", type=str,
                                help='By default, the validator uses the filename to match '
                                     'the file format.  If your filename is incorrectly named, '
                                     'it will be invalid.  If you know the file format you are '
                                     'validating, you can ignore the filename validation and skip '
                                     'to file content validation. '
                                     'Note, the filetypes with SP at '
                                     'the end are for special sponsored projects.')

    validate_group.add_argument("--parentid", type=str, default=None,
                                help='Synapse id of center input folder. '
                                     'If specified, your valid files will be uploaded '
                                     'to this directory.')

    # TODO: remove this default when private genie project is ready
    parser_validate.add_argument("--project_id", type=str,
                                 default="syn3380222",
                                 help='Synapse Project ID where data is stored.')

    parser_validate.add_argument("--nosymbol-check", action='store_true',
                                 help='Do not check hugo symbols of fusion and cna file')

    parser_validate.set_defaults(func=validate._perform_validate)

    parser_bootstrap = subparsers.add_parser('bootstrap-infra',
                                            help='Create GENIE-like infra')
    parser_bootstrap.add_argument(
        "--format_registry_packages", type=str, nargs="+",
        default=["example_registry"],
        help="Python package name(s) to get valid file formats from "
             "(default: %(default)s)."
    )

    parser_bootstrap.set_defaults(func=bootstrap_infra)

    parser_process = subparsers.add_parser(
        'process', help='Process files'
    )
    parser_process.add_argument(
        "process", choices=['main']
    )
    parser_process.add_argument(
        "--project_id",
        help="Synapse Project ID where data is stored.",
        required=True
    )
    parser_process.add_argument(
        '--center', help='The centers'
    )
    parser_process.add_argument(
        "--pemfile", type=str,
        help="Path to PEM file (genie.pem)"
    )
    parser_process.add_argument(
        "--delete_old", action='store_true',
        help="Delete all old processed and temp files"
    )
    parser_process.add_argument(
        "--only_validate", action='store_true',
        help="Only validate the files, don't process"
    )
    parser_process.add_argument(
        "--debug", action='store_true',
        help="Add debug mode to synapse"
    )
    # DEFAULT PARAMS
    parser_process.add_argument(
        "--format_registry_packages", type=str, nargs="+",
        default=["example_registry"],
        help="Python package name(s) to get valid file formats from "
             "(default: %(default)s)."
    )
    parser_process.set_defaults(func=process_cli_wrapper)

    return parser


def main():
    """Invoke"""
    args = build_parser().parse_args()
    syn = synapse_login(args.syn_user, args.syn_pass)
    # func has to match the set_defaults
    args.func(syn, args)



if __name__ == "__main__":
    main()
