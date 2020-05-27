"""Write invalid reasons"""
import logging
import os

import synapseclient

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def write_file_invalid_reasons(x, syn, error_file):
    ent = syn.get(x['id'], downloadFile=False)
    errors = x['errors']
    showname = "\t" + ent.name + " (%s):\n\n" % ent.id
    error_file.write(showname)
    error_file.write(errors + "\n\n\n")


def write(syn, center_mapping_df, error_tracker_synid):
    for center in center_mapping_df['center']:
        logger.info(center)
        staging_synid = center_mapping_df['stagingSynId'][
            center_mapping_df['center'] == center][0]
        with open(center + "_errors.txt", 'w') as errorfile:
            error_tracker = syn.tableQuery(
                f"SELECT * FROM {error_tracker_synid} where "
                f"center = '{center}'"
            )
            error_trackerdf = error_tracker.asDataFrame()
            error_trackerdf.apply(
                lambda row: write_file_invalid_reasons(row, syn, errorfile),
                axis=1
            )
            if error_trackerdf.empty:
                errorfile.write("No errors!")
        ent = synapseclient.File(center + "_errors.txt",
                                 parentId=staging_synid)
        syn.store(ent)
        os.remove(center + "_errors.txt")
