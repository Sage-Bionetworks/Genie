#! /usr/bin/env python3
import argparse
import logging

from genie import input_to_database
from genie import process_functions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main(db_mapping_id, pemfile):

    syn = process_functions.synLogin(pemfile, debug=debug)

    query = 'SELECT * FROM {}'.format(db_mapping_id)
    databaseToSynIdMapping = syn.tableQuery(query)
    databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()

    databaseToSynIdMappingDf = \
            input_to_database.create_and_archive_maf_database(syn, databaseToSynIdMappingDf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--db_mapping_id', type=str,
                        help='Database to ID mapping table Synapse ID.')
    parser.add_argument("--pemFile", type=str,
                        help="Path to PEM file (genie.pem)")
                     
    args = parser.parse_args()

    main(db_mapping_id=args.db_mapping_id, pemfile=args.pemFile)
