import argparse
import logging
import os

import synapseclient

from . import process_functions

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def write_file_invalid_reasons(x, syn, error_file):
    ent = syn.get(x['id'],downloadFile=False)
    errors = x['errors']
    showname = "\t" + ent.name + " (%s):\n\n" % ent.id
    error_file.write(showname)
    error_file.write(errors + "\n\n\n")

def write_invalid_reasons(syn, center_mapping_df, error_tracker_synid):
    for center in center_mapping_df['center']:
        logger.info(center)
        staging_synid = center_mapping_df['stagingSynId'][center_mapping_df['center'] == center][0]
        with open(center + "_errors.txt", 'w') as errorfile:
            error_tracker = syn.tableQuery("SELECT * FROM %s where center = '%s'"  % (error_tracker_synid, center))
            error_trackerdf = error_tracker.asDataFrame()
            error_trackerdf.apply(lambda x: write_file_invalid_reasons(x, syn, errorfile),axis=1)
            if error_trackerdf.empty:
                errorfile.write("No errors!")
        ent = synapseclient.File(center + "_errors.txt", parentId = staging_synid)
        syn.store(ent)
        os.remove(center + "_errors.txt")   


def main():
    parser = argparse.ArgumentParser(description='Write invalid reasons')
    parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")
    parser.add_argument("--debug", action='store_true',help="Synapse Debug Feature")
    args = parser.parse_args()
    
    syn = process_functions.synLogin(args.pemFile, debug=args.debug)
    center_mapping = syn.tableQuery('SELECT * FROM syn10061452 where inputSynId is not null and release is true')
    center_mappingdf = center_mapping.asDataFrame()
    error_tracker_synid = "syn10153306"
    write_invalid_reasons(syn, center_mappingdf, error_tracker_synid)

if __name__ == "__main__":
    main()