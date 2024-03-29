#! /usr/bin/env python
import argparse

import pandas as pd

from . import extract, process_functions


def retract_samples(syn, database_synid, col, remove_values):
    """
    Helper retraction function that help remove values from a column in
    a database

    params:
        syn: synapse object
        databse_synid:  synapse id of database
        col: column in database
        remove_values: list of values to remove from the column
    """
    # schema = syn.get(database_synid)
    remove_values_query = "','".join(remove_values)
    remove_rows = syn.tableQuery(
        "select %s from %s where %s in ('%s')"
        % (col, database_synid, col, remove_values_query)
    )
    if len(remove_rows.asRowSet()["rows"]) > 0:
        syn.delete(remove_rows.asRowSet())
    else:
        print("Nothing to retract")


def retract(syn, project_id):
    """
    Main retraction function

    params:
        syn: synapse object
        project_id: Synapse Project ID with a database mapping table
    """

    pat_retraction_synid = extract.getDatabaseSynId(
        syn=syn, tableName="patientRetraction", project_id=project_id
    )
    patientRetractIds = extract.get_syntabledf(
        syn=syn, query_string=f"select * from {pat_retraction_synid}"
    )
    # grab all clinical samples that belong to patients in the patient
    # clinical file and append to sample list
    sample_synid = extract.getDatabaseSynId(
        syn=syn, tableName="sample", project_id=project_id
    )
    sampleClinicalDf = extract.get_syntabledf(
        syn=syn, query_string=f"select * from {sample_synid}"
    )

    appendSamples = sampleClinicalDf["SAMPLE_ID"][
        sampleClinicalDf["PATIENT_ID"].isin(patientRetractIds.geniePatientId)
    ]

    sample_retract_synid = extract.getDatabaseSynId(
        syn=syn, tableName="sampleRetraction", project_id=project_id
    )
    sampleRetractIds = extract.get_syntabledf(
        syn=syn, query_string=f"select * from {sample_retract_synid}"
    )
    allRetractedSamples = pd.concat([sampleRetractIds["genieSampleId"], appendSamples])

    # Only need to retract clinical data, because the rest of the data is filtered by clinical data
    # Sample Clinical Data
    retract_samples(
        syn=syn,
        database_synid=sample_synid,
        col="SAMPLE_ID",
        remove_values=allRetractedSamples,
    )
    # Patient Clinical Data
    retract_samples(
        syn=syn,
        database_synid=extract.getDatabaseSynId(
            syn, tableName="patient", project_id=project_id
        ),
        col="PATIENT_ID",
        remove_values=patientRetractIds["geniePatientId"],
    )


# TODO: move to cli
def main():
    """
    Main block with argparse and calls the main retract function
    """
    parser = argparse.ArgumentParser(description="Sample retraction")
    parser.add_argument("--project_id", type=str, help="Synapse Project ID to use.")
    parser.add_argument("--debug", action="store_true", help="Synapse Debug Feature")
    args = parser.parse_args()
    syn = process_functions.synapse_login(debug=args.debug)
    retract(syn, args.project_id)


# TODO: Remove and move to cli
if __name__ == "__main__":
    main()
