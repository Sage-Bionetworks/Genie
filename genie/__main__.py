#!/usr/bin/env python3
# noqa pylint: disable=line-too-long
"""genie cli"""
import argparse
import logging

import synapseclient

# import genie.config
from . import create_case_lists, validate
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


def perform_create_case_list(syn, args):
    """CLI to create case lists given a clinical and assay information file"""
    create_case_lists.main(args.clinical_file_name,
                           args.assay_info_file_name,
                           args.output_dir,
                           args.study_id)



def build_parser():
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

    parser_create_case = subparsers.add_parser(
        'create-case-lists', help='Creates cBioPortal case list files'
    )
    parser_create_case.add_argument("clinical_file_name",
                                    type=str,
                                    help="Clinical file path")
    parser_create_case.add_argument("assay_info_file_name",
                                    type=str,
                                    help="gene matrix file path")
    parser_create_case.add_argument("output_dir",
                                    type=str,
                                    help="Output directory")
    parser_create_case.add_argument("study_id",
                                    type=str,
                                    help="Output directory")
    parser_create_case.set_defaults(func=perform_create_case_list)

    return parser


def main():
    """Invoke"""
    args = build_parser().parse_args()
    syn = synapse_login(args.syn_user, args.syn_pass)
    # func has to match the set_defaults
    args.func(syn, args)


if __name__ == "__main__":
    main()
