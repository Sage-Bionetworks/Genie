#! /usr/bin/env python3
import argparse
import csv
import json
import logging
import multiprocessing.dummy

import synapseclient

from genie import input_to_database
from genie import process_functions

logger = logging.getLogger(__name__)

def get_rows(view_result):
    rows = csv.DictReader(open(view_result.filepath))

    for ent_row in rows:
        # pop these
        for pop_this in ['ROW_ID', 'ROW_VERSION', 'ROW_ETAG']:
            ent_row.pop(pop_this)

        yield dict(id=ent_row['id'], name=ent_row['name'])

def get_entity(ent_row, syn=None, downloadFile=False):
    """Utility to get an entity remotely or create a local one.

    """

    if syn:
        ent_synid = ent_row['id']
        ent = syn.get(ent_synid, downloadFile=downloadFile)
    else:
        ent = synapseclient.File(**ent_row)
    return ent

def get_center_input_files(syn, synid, center, db_mapping_id,
                           synapse_entity=True, downloadFile=True):
    '''
    This function uses a file view to get a list of tuples of center files.

    Args:
        syn: Synapse object
        viewid: Synapse File View ID that lists all files
        synid: Center input folder synid
        center: Center name
    Returns:
        List of entities with the correct format to pass into validation
    '''
    logger.info("GETTING {center} INPUT FILES".format(center=center))
    clinical_pair_name = [
        "data_clinical_supp_sample_{center}.txt".format(center=center),
        "data_clinical_supp_patient_{center}.txt".format(center=center)]

    view_res = get_tracking_table(syn, db_mapping_id, "centerFileView",
                                  "WHERE parentId = '{synid}'".format(synid=synid))

    logger.debug(f"Got view to {view_res.filepath}")

    center_files = get_rows(view_res)
    logger.debug("Opened view")

    clinicalpair_entities = []
    prepared_center_file_list = []

    if synapse_entity:
        syn_conn = syn
    else:
        syn_conn = None

    pool = multiprocessing.dummy.Pool(8)

    ents = pool.map(lambda x: get_entity(x, syn=syn_conn, 
                                         downloadFile=downloadFile),
                    center_files)

    for ent in ents:
        # Clinical file can come as two files.
        # The two files need to be merged together which is
        # why there is this format

        if ent.name in clinical_pair_name:
            clinicalpair_entities.append([ent])
            continue

        prepared_center_file_list.append([ent])

    if clinicalpair_entities:
        # clinicalpair_entities = [x for x in clinicalpair]
        prepared_center_file_list.append(clinicalpair_entities)

    return prepared_center_file_list

def check_existing_file_status(validation_status_table, error_tracker_table, entities):
    '''
    This function checks input files against the existing validation and error
    tracking dataframe

    Args:
        validation_status_table: Validation status Synapse Table query result
        error_tracker_table: Error tracking Synapse Table query result
        entities: list of center input entites

    Returns:
        dict: Input file status
            status_list: file validation status
            error_list: Errors of the files if they exist,
            to_validate: Boolean value for whether of not an input
                         file needs to be validated
    '''
    if len(entities) > 2:
        raise ValueError(
            "There should never be more than 2 files being validated.")

    statuses = []
    errors = []

    validation_statusdf = validation_status_table.asDataFrame()
    error_trackerdf = error_tracker_table.asDataFrame()

    for ent in entities:
        to_validate = False
        # Get the current status and errors from the tables.
        current_status = validation_statusdf[validation_statusdf['id'] == ent.id]
        current_error = error_trackerdf[error_trackerdf['id'] == ent.id]

        if current_status.empty:
            to_validate = True
        else:
            # This to_validate is here, because the following is a
            # sequential check of whether files need to be validated
            statuses.append(current_status['status'].values[0])
            if current_error.empty:
                to_validate = \
                    current_status['status'].values[0] == "INVALID"
            else:
                errors.append(current_error['errors'].values[0])
            # Add Name check here (must add name of the entity as a column)
            if current_status['md5'].values[0] != ent.md5 or \
               current_status['name'].values[0] != ent.name:
                to_validate = True
            else:
                status_str = "{filename} ({id}) FILE STATUS IS: {filestatus}"
                logger.debug(status_str.format(filename=ent.name, id=ent.id,
                                              filestatus=current_status['status'].values[0]))

    return({
        'status_list': statuses,
        'error_list': errors,
        'to_validate': to_validate})



def flatten_status(status_obj):
    ents = status_obj['entity']
    statuses = status_obj['status']['status_list']
    to_validate = status_obj['status']['to_validate']
    error_list = status_obj['status']['error_list']

    assert len(ents) == len(statuses), "Number of entities and statuses not equal."
    
    for (ent, status, error) in zip(ents, statuses, error_list):
        yield(dict(entity=ent, status=status, error=error, to_validate=to_validate))

def get_tracking_table(syn, db_mapping_id, table_name, where_clause=None):
    """Get the contents of one of the status tables in the database mapping.
    """

    logger.debug(f"Getting {table_name} from db mapping table {db_mapping_id}")

    db_mapping_tbl = syn.tableQuery(
        'SELECT * FROM {}'.format(db_mapping_id))
    db_mapping_df = db_mapping_tbl.asDataFrame()

    logger.debug(f"Got db mapping table {db_mapping_id}")

    synid = process_functions.getDatabaseSynId(
        syn, table_name,
        databaseToSynIdMappingDf=db_mapping_df)

    logger.debug(f"Got id of {table_name} table: {synid}")

    where = where_clause if where_clause is not None else ""

    tbl = syn.tableQuery(
        "SELECT * FROM {synid} {where}".format(synid=synid, where=where))

    logger.debug(f"Ran table query to get tracking table {synid}")

    return tbl

def main(project, center=None):

    syn = synapseclient.login(silent=True)
    syn_project = syn.get(project)
    db_mapping_synid = syn_project.annotations.dbMapping[0]

    center_mapping_tbl = get_tracking_table(syn,
                                            db_mapping_id=db_mapping_synid,
                                            table_name="centerMapping")
    center_mapping_df = center_mapping_tbl.asDataFrame()
    logger.info("Got center mapping table")
    if center is not None:
        assert center in center_mapping_df.center.tolist(), (
            "Must specify one of these centers: {}".format(
                ", ".join(center_mapping_df.center)))
        center_mapping_df = center_mapping_df[center_mapping_df.center == center]
    else:
        center_mapping_df = \
            center_mapping_df[~center_mapping_df['inputSynId'].isnull()]

    logger.debug(center_mapping_df)

    syn.table_query_timeout = 50000

    status_table = get_tracking_table(syn,
                                      db_mapping_id=db_mapping_synid,
                                      table_name="validationStatus")

    error_table = get_tracking_table(syn,
                                     db_mapping_id=db_mapping_synid,
                                     table_name="errorTracker")

    entity_status_list = []

    pool = multiprocessing.dummy.Pool(8)

    for index, row in center_mapping_df.iterrows():
        center = row['center']
        synid = row['inputSynId']
        logger.debug(f"{index} {center} {synid}")
        input_files = get_center_input_files(syn,
                                             synid=synid,
                                             db_mapping_id=db_mapping_synid,
                                             center=center,
                                             downloadFile=False)


        statuses = pool.map(
            lambda x: input_to_database.check_existing_file_status(entities=x,
                                                                   validation_status_table=status_table,
                                                                   error_tracker_table=error_table),
                        input_files)

        for entity, status in zip(input_files, statuses):
            entity_status_list.append(dict(entity=entity, status=status))

    for obj in entity_status_list:
        flat_obj = flatten_status(obj)
        for row in flat_obj:
            res = dict(entity_id=row['entity'].id,
                       status=row['status'], 
                       to_validate=row['to_validate'])
            print(json.dumps(res))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get input files.')
    parser.add_argument('--projectId', type=str,
                        help='The Synapse Project for data submission and validation.')
    parser.add_argument('--center', help='The centers', type=str, default=None)
    args = parser.parse_args()

    main(project=args.projectId, center=args.center)
