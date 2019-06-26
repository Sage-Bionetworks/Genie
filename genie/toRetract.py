#! /usr/bin/env python
import argparse

import synapseclient

from . import process_functions

def retract_samples(syn, database_synid, col, remove_values):
    '''
    Helper retraction function that help remove values from a column in
    a database

    params:
        syn: synapse object
        databse_synid:  synapse id of database
        col: column in database
        remove_values: list of values to remove from the column 
    '''

    schema = syn.get(database_synid)
    remove_values_query = "','".join(remove_values)
    remove_rows = syn.tableQuery("select %s from %s where %s in ('%s')" % (col, database_synid, col, remove_values_query))
    if len(remove_rows.asRowSet()['rows']) > 0:
        syn.delete(remove_rows.asRowSet())
    else:
        print("Nothing to retract")


def retract(syn, test=False):
    '''
    Main retraction function
    
    params:
        syn: synapse object
        test: use test files or main files. Default is False
    '''

    patientRetract = syn.tableQuery('select * from %s' % process_functions.getDatabaseSynId(syn, "patientRetraction", test=test))
    patientRetractIds = patientRetract.asDataFrame()
    #grab all clinical samples that belong to patients in the patient clinical file and append to sample list
    sampleClinical = syn.tableQuery('select * from %s' % process_functions.getDatabaseSynId(syn, "sample", test=test))
    sampleClinicalDf = sampleClinical.asDataFrame()
    appendSamples = sampleClinicalDf['SAMPLE_ID'][sampleClinicalDf['PATIENT_ID'].isin(patientRetractIds.geniePatientId)]
    
    sampleRetract = syn.tableQuery('select * from %s' % process_functions.getDatabaseSynId(syn, "sampleRetraction", test=test))
    sampleRetractIds = sampleRetract.asDataFrame()

    allRetractedSamples = sampleRetractIds['genieSampleId'].append(appendSamples)

    #Only need to retract clinical data, because the rest of the data is filtered by clinical data
    #Sample Clinical Data
    retract_samples(syn,process_functions.getDatabaseSynId(syn,"sample", test=test),"SAMPLE_ID",allRetractedSamples)
    #Patient Clinical Data
    retract_samples(syn, process_functions.getDatabaseSynId(syn, "patient", test=test),"PATIENT_ID",patientRetractIds['geniePatientId'])

def main():
    '''
    Main block with argparse and calls the main retract function
    '''
    parser = argparse.ArgumentParser(description='Sample retraction')
    parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")
    parser.add_argument("--test", action='store_true',help="Run test")
    parser.add_argument("--debug", action='store_true',help="Synapse Debug Feature")
    args = parser.parse_args()
    syn = process_functions.synLogin(args.pemFile, debug=args.debug)
    retract(syn, args.test)

if __name__ == "__main__":
    main()


