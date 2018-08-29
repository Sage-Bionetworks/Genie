import synapseclient
import pandas as pd
import argparse
import os
from genie import process_functions

def retraction(syn, databaseSynId, sampleKey, deleteSamples):
	schema = syn.get(databaseSynId)
	deleteSamplesQuery = "','".join(deleteSamples)
	database = syn.tableQuery("select %s from %s where %s in ('%s')" % (sampleKey, databaseSynId, sampleKey, deleteSamplesQuery))
	if len(database.asRowSet()['rows']) > 0:
		syn.delete(database.asRowSet())
	else:
		print("Nothing to retract")

def getDatabaseSynId(syn, tableName, test=False):
	if test:
		databaseToSynIdMappingSynId = 'syn11600968'
	else:
		databaseToSynIdMappingSynId = 'syn10967259'
	databaseToSynIdMapping = syn.tableQuery('SELECT * FROM %s' % databaseToSynIdMappingSynId)
	databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()
	synId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == tableName]
	return(synId.values[0])

def main():
	parser = argparse.ArgumentParser(description='Sample retraction')
	parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")
	parser.add_argument("--test", action='store_true',help="Run test")
	args = parser.parse_args()
	syn = process_functions.synLogin(args)

	patientRetract = syn.tableQuery('select * from %s' % getDatabaseSynId(syn, "patientRetraction", test=args.test))
	patientRetractIds = patientRetract.asDataFrame()
	#grab all clinical samples that belong to patients in the patient clinical file and append to sample list
	sampleClinical = syn.tableQuery('select * from %s' % getDatabaseSynId(syn, "sample", test=args.test))
	sampleClinicalDf = sampleClinical.asDataFrame()
	appendSamples = sampleClinicalDf['SAMPLE_ID'][sampleClinicalDf['PATIENT_ID'].isin(patientRetractIds.geniePatientId)]
	
	sampleRetract = syn.tableQuery('select * from %s' % getDatabaseSynId(syn, "sampleRetraction", test=args.test))
	sampleRetractIds = sampleRetract.asDataFrame()

	allRetractedSamples = sampleRetractIds['genieSampleId'].append(appendSamples)

	#Sample Clinical Data
	retraction(syn,getDatabaseSynId(syn,"sample", test=args.test),"SAMPLE_ID",allRetractedSamples)
	#Patient Clinical Data
	retraction(syn, getDatabaseSynId(syn, "patient", test=args.test),"PATIENT_ID",patientRetractIds['geniePatientId'])
	#Fusions
	#retraction(syn,getDatabaseSynId(syn, "fusion", test=args.test),"TUMOR_SAMPLE_BARCODE",allRetractedSamples)
	#CNA
	#retraction(syn,getDatabaseSynId(syn, "cna", test=args.test),"TUMOR_SAMPLE_BARCODE",allRetractedSamples)
	#VCF2MAF
	#retraction(syn,getDatabaseSynId(syn, "vcf2maf", test=args.test),"Tumor_Sample_Barcode",allRetractedSamples)
	#Vital status
	#retraction(syn, getDatabaseSynId(syn, "vitalStatus", test=args.test),"PATIENT_ID",patientRetractIds['geniePatientId'])
	#SEG
	#retraction(syn,getDatabaseSynId(syn, "seg", test=args.test),"ID",allRetractedSamples)

if __name__ == "__main__":
	main()


