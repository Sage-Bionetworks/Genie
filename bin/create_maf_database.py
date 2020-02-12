"""Script to invoke creating maf database"""
#! /usr/bin/env python3
import argparse
import logging
import time

import synapseclient

from genie import process_functions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# def _create_maf_db(syn, foo):
#     maf_database_ent = syn.get(maf_database_synid)
#     print(maf_database_ent)
#     maf_columns = list(syn.getTableColumns(maf_database_synid))
#     schema = synapseclient.Schema(
#         name='Narrow MAF {current_time} Database'.format(
#             current_time=time.time()),
#         columns=maf_columns,
#         parent=process_functions.getDatabaseSynId(
#             syn, "main",
#             databaseToSynIdMappingDf=database_synid_mappingdf))
#     schema.primaryKey = maf_database_ent.primaryKey
#     new_maf_database = syn.store(schema)

# TODO: Should split this into 3 funcitons
# so that unit tests are easier to write


def create_and_archive_maf_database(syn, database_synid_mappingdf):
    '''
    Creates new MAF database and archives the old database in the staging site

    Args:
        syn: Synapse object
        databaseToSynIdMappingDf: Database to synapse id mapping dataframe

    Return:
        Editted database to synapse id mapping dataframe
    '''
    maf_database_synid = process_functions.getDatabaseSynId(
        syn, "vcf2maf", databaseToSynIdMappingDf=database_synid_mappingdf)
    maf_database_ent = syn.get(maf_database_synid)
    maf_columns = list(syn.getTableColumns(maf_database_synid))
    schema = synapseclient.Schema(
        name='Narrow MAF {current_time} Database'.format(
            current_time=time.time()),
        columns=maf_columns,
        parent=process_functions.getDatabaseSynId(
            syn, "main", databaseToSynIdMappingDf=database_synid_mappingdf))
    schema.primaryKey = maf_database_ent.primaryKey
    new_maf_database = syn.store(schema)
    # Store in the new database synid
    database_synid_mappingdf['Id'][
        database_synid_mappingdf[
            'Database'] == 'vcf2maf'] = new_maf_database.id
    # Only get the vcf2maf row so that only this row is updated in the
    # mapping table
    vcf2maf_mappingdf = database_synid_mappingdf[
        database_synid_mappingdf['Database'] == 'vcf2maf']

    # Update this synid later (This synid needs to not be hardcoded)
    syn.store(synapseclient.Table("syn10967259", vcf2maf_mappingdf))
    # Move and archive old mafdatabase (This is the staging synid)
    maf_database_ent.parentId = "syn7208886"
    maf_database_ent.name = "ARCHIVED " + maf_database_ent.name
    syn.store(maf_database_ent)
    # maf_database_synid = new_maf_database.id
    # Remove can download permissions from project GENIE team
    syn.setPermissions(new_maf_database.id, 3326313, [])
    return database_synid_mappingdf


def main(db_mapping_id, pemfile, debug):
    """Creates maf database"""
    syn = process_functions.synLogin(pemfile, debug=debug)

    query = 'SELECT * FROM {}'.format(db_mapping_id)
    database_mapping = syn.tableQuery(query)
    database_mappingdf = database_mapping.asDataFrame()
    database_mappingdf = create_and_archive_maf_database(syn,
                                                         database_mappingdf)


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
