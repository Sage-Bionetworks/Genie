"""Script to invoke creating maf database"""
#! /usr/bin/env python3
import argparse
import logging
import time

from genie import process_functions


def main(filetype: str, projectid: str, archive_projectid: str = None,
         pemfile: str = None, debug: str = False):
    """Creates maf database

    Args:
        filetype: Type of table to create
        projectid: Project Id
        achive_projectid: Archive Project id
        pemfile: Path to PEM file (genie.pem)
        debug: Add debug mode to synapse
    """
    syn = process_functions.synLogin(pemfile, debug=debug)
    table_name = 'Narrow MAF {} Database'.format(time.time())
    filetype = "vcf2maf"
    if archive_projectid is None:
        archive_projectid = projectid

    new_tables = process_functions.create_new_fileformat_table(
        syn, filetype, table_name, projectid, archive_projectid
    )
    # Remove can download permissions from project GENIE team
    # this wouldn't be necessary if the database is just moved to a project
    # that isn't viewable by GENIE peeps.
    syn.setPermissions(new_tables['newdb_ent'].id, 3326313, [])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('project_id', type=str,
                        help='Synapse Project ID where data is stored.',
                        required=True)
    parser.add_argument(
        '--archive_project_id', type=str,
        help='Synapse Project ID where you want old table to be stored'
             'Defaults to project_id'
    )
    parser.add_argument("--pemFile", type=str,
                        help="Path to PEM file (genie.pem)")
    parser.add_argument("--debug", action='store_true',
                        help="Add debug mode to synapse")
    args = parser.parse_args()

    main(filetype=args.filetype, projectid=args.project_id,
         archive_projectid=args.archive_project_id,
         pemfile=args.pemFile, debug=args.debug)
