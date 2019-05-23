#!/usr/bin/env python
import argparse
import genie


def perform_validate(syn, args):
    

def perform_main(syn, args):
    if 'func' in args:
        try:
            args.func(syn, args)
        except Exception:
            raise


def build_parser():
    parser = argparse.ArgumentParser(description='GENIE processing')

    subparsers = parser.add_subparsers(
        title='commands',
        description='The following commands are available:',
        help='For additional help: "genie <COMMAND> -h"')

    parser_validate = subparsers.add_parser(
        'validate',
        help='Validates a GENIE file')

    parser_validate.add_argument(
        "fileType",
        type=str,
        choices=genie.PROCESS_FILES.keys(),
        help='Filetypes that you are validating. Note, the filetypes with SP at \
            the end are for special sponsored projects')

    parser_validate.add_argument(
        "file",
        type=str,
        nargs="+",
        help='File(s) that you are validating.  \
        If you validation your clinical files and you have both sample and \
        patient files, you must provide both')

    parser_validate.add_argument(
        "center",
        type=str,
        help='Contributing Centers')

    parser_validate.add_argument(
        "--thread",
        type=int,
        required=False,
        default=1,
        help='Number of threads used in validation symbols')

    parser_validate.add_argument(
        "--offline",
        action='store_true',
        help='No validation of filenames')

    parser_validate.add_argument(
        "--uploadToSynapse",
        type=str,
        default=None,
        help='Will upload the file to the synapse directory of users choice')

    parser_validate.add_argument(
        "--oncotreeLink",
        type=str,
        help="Link to oncotree code")

    parser_validate.add_argument(
        "--noSymbolCheck",
        action='store_true',
        help='Do not check hugo symbols of fusion and cna file')

    parser_validate.add_argument(
        "--testing",
        action='store_true',
        help='Put in testing mode')

    parser_validate.set_defaults(func=genie.validate.perform_validate)
    return(parser)


def main():
    args = build_parser().parse_args()
    syn = genie.validate.synapse_login()
    perform_main(syn, args)


if __name__ == "__main__":
    main()
