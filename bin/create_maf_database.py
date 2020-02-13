"""Script to invoke creating maf database"""
#! /usr/bin/env python3
import argparse
import logging
import time

from genie import process_functions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# This function can further be abstracted once things are further fixed
def main(db_mapping_id, pemfile, debug):
    """Creates maf database"""
    syn = process_functions.synLogin(pemfile, debug=debug)

    query = 'SELECT * FROM {}'.format(db_mapping_id)
    database_mapping = syn.tableQuery(query)
    database_mappingdf = database_mapping.asDataFrame()
    # New table name
    table_name = 'Narrow MAF {} Database'.format(time.time())

    #Get main project id
    projectid = process_functions.getDatabaseSynId(
        syn, "main", databaseToSynIdMappingDf=database_mappingdf)

    new_tables = process_functions.create_new_fileformat_table(syn, database_mapping,
                                                               "vcf2maf",
                                                               table_name,
                                                               projectid,
                                                               "syn7208886")

    # Remove can download permissions from project GENIE team
    # this wouldn't be necessary if the database is just moved to a project
    # that isn't viewable by GENIE peeps.
    syn.setPermissions(new_tables['newdb_ent'].id, 3326313, [])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--db_mapping_id', type=str,
                        help='Database to ID mapping table Synapse ID.')
    parser.add_argument("--pemFile", type=str,
                        help="Path to PEM file (genie.pem)")
    parser.add_argument("--debug", action='store_true',
                        help="Add debug mode to synapse")

    args = parser.parse_args()

    main(db_mapping_id=args.db_mapping_id, pemfile=args.pemFile,
         debug=args.debug)
