import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys
from genie.clinical import clinical


def test_processing():
	def createMockTable(dataframe):
		table = mock.create_autospec(synapseclient.table.CsvFileTable)
		table.asDataFrame.return_value= dataframe
		return(table)

	def table_query_results(*args):
		return(table_query_results_map[args])

	withNan = pd.DataFrame(dict(CODE=[1,2,3,4,float('nan')],
					  	        CBIO_LABEL=['Test','Why','foo','Me','Unknown'],
					  	        DESCRIPTION=['non','asdf','asdf','asdff','asdfasdf']))

	noNan = pd.DataFrame(dict(CODE=[1,2,3,4],
					  	      CBIO_LABEL=['Test','Why','foo','Me'],
					  	      DESCRIPTION=['non','asdf','asdf','asdff']))

	table_query_results_map = {
	("SELECT * FROM syn7434222",) : createMockTable(withNan),
	("SELECT * FROM syn7434236",) : createMockTable(withNan),
	("SELECT * FROM syn7434242",) : createMockTable(withNan),
	("SELECT * FROM syn7434273",) : createMockTable(noNan)
	}	

	syn = mock.create_autospec(synapseclient.Synapse) 
	syn.tableQuery.side_effect=table_query_results

	clin = clinical(syn, "SAGE")

	patientCols = ["PATIENT_ID","SEX","PRIMARY_RACE","SECONDARY_RACE",
				   "TERTIARY_RACE","ETHNICITY","BIRTH_YEAR","CENTER"]
	sampleCols = ["SAMPLE_ID","PATIENT_ID","AGE_AT_SEQ_REPORT","ONCOTREE_CODE","SAMPLE_TYPE",
				  "SEQ_ASSAY_ID",'SEQ_DATE','SAMPLE_TYPE_DETAILED','SEQ_YEAR']

	expectedPatientDf = pd.DataFrame(dict(PATIENT_ID=["GENIE-SAGE-ID1","GENIE-SAGE-ID2","GENIE-SAGE-ID3","GENIE-SAGE-ID4","GENIE-SAGE-ID5"],
										 SEX=['Test','Why','foo','Me','Unknown'],
										 PRIMARY_RACE=['Test','Why','foo','Me','Unknown'],
										 SECONDARY_RACE=['Test','Why','foo','Me','Unknown'],
										 TERTIARY_RACE=['Test','Why','foo','Me','Unknown'],
										 ETHNICITY=['Test','Why','foo','Me','Unknown'],
										 BIRTH_YEAR=[1990,1990,1990,1990,1990],
										 CENTER=["SAGE","SAGE","SAGE","SAGE","SAGE"]))
	#TEST patient processing
	patientDf = pd.DataFrame(dict(PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								 SEX=[1,2,3,4,float('nan')],
								 PRIMARY_RACE=[1,2,3,4,float('nan')],
								 SECONDARY_RACE=[1,2,3,4,float('nan')],
								 TERTIARY_RACE=[1,2,3,4,float('nan')],
								 ETHNICITY=[1,2,3,4,float('nan')],
								 BIRTH_YEAR=[1990,1990,1990,1990,1990],
								 CENTER=["FOO","FOO","FOO","FOO","FOO"]))
	clinicalTemplate = pd.DataFrame(columns=patientCols)
	patient=True
	newPatientDf = clin._process(patientDf, clinicalTemplate)

	assert expectedPatientDf.equals(newPatientDf)

	#TEST for sample processing
	expectedSampleDf = pd.DataFrame(dict(SAMPLE_ID=["GENIE-SAGE-ID1-1","GENIE-SAGE-ID2-1","GENIE-SAGE-ID3-1","GENIE-SAGE-ID4-1","GENIE-SAGE-ID5-1"],
										 PATIENT_ID=["GENIE-SAGE-ID1","GENIE-SAGE-ID2","GENIE-SAGE-ID3","GENIE-SAGE-ID4","GENIE-SAGE-ID5"],
										 AGE_AT_SEQ_REPORT=[100000,100000,100000,100000,100000],
										 ONCOTREE_CODE=['AMPCA','AMPCA','AMPCA','AMPCA','AMPCA'],
										 SAMPLE_TYPE=['Test','Why','foo','Me','Me'],
										 CENTER=["SAGE","SAGE","SAGE","SAGE","SAGE"],
										 SAMPLE_TYPE_DETAILED=['non','asdf','asdf','asdff','asdff'],
										 SEQ_ASSAY_ID=['SAGE-1','SAGE-1','SAGE-1','SAGE-1','SAGE-1'],
										 SEQ_DATE=['Jan-2012','Apr-2013','Jul-2014','Oct-2015','Release'],
										 SEQ_YEAR=[2012,2013,2014,2015,pd.np.nan]))
	clinicalTemplate = pd.DataFrame(columns=sampleCols)
	patient=False
	sampleDf = pd.DataFrame(dict(SAMPLE_ID=["ID1-1","ID2-1","ID3-1","ID4-1","ID5-1"],
								 PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								 AGE_AT_SEQ_REPORT=[100000,100000,100000,100000,100000],
								 ONCOTREE_CODE=['AMPCA','AMPCA','AMPCA','AMPCA','AMPCA'],
								 SAMPLE_TYPE=[1,2,3,4,4],
								 SEQ_ASSAY_ID=['SAGE-1','SAGE-1','SAGE-1','SAGE-1','SAGE-1'],
								 SEQ_DATE=['Jan-2012','Apr-2013','JUL-2014','Oct-2015','release']))

	newSampleDf = clin._process(sampleDf, clinicalTemplate)
	assert expectedSampleDf.equals(newSampleDf[expectedSampleDf.columns])

def test_validation():
	def createMockTable(dataframe):
		table = mock.create_autospec(synapseclient.table.CsvFileTable)
		table.asDataFrame.return_value= dataframe
		return(table)

	def table_query_results(*args):
		return(table_query_results_map[args])

	withNan = pd.DataFrame(dict(CODE=[1,2,3,4,float('nan')],
					  	        CBIO_LABEL=['Test','Why','foo','Me','Unknown'],
					  	        DESCRIPTION=['non','asdf','asdf','asdff','asdfasdf']))

	noNan = pd.DataFrame(dict(CODE=[1,2,3,4],
					  	      CBIO_LABEL=['Test','Why','foo','Me'],
					  	      DESCRIPTION=['non','asdf','asdf','asdff']))
	
	sexDf = pd.DataFrame(dict(CODE=[1,2,float('nan')],
				  	        CBIO_LABEL=['Male','Female','Unknown'],
				  	        DESCRIPTION=['Male','Female','Not coded']))


	table_query_results_map = {
	("SELECT * FROM syn7434222",) : createMockTable(sexDf),
	("SELECT * FROM syn7434236",) : createMockTable(withNan),
	("SELECT * FROM syn7434242",) : createMockTable(withNan),
	("SELECT * FROM syn7434273",) : createMockTable(noNan)
	}	

	syn = mock.create_autospec(synapseclient.Synapse) 
	syn.tableQuery.side_effect=table_query_results

	oncotree_url = 'http://oncotree.mskcc.org/api/tumor_types.txt?version=oncotree_latest_stable'
	json_oncotreeurl = "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"
	clin = clinical(syn, "SAGE")
	assert_raises(AssertionError, clin.validateFilename, ["foo"])
	assert_raises(AssertionError, clin.validateFilename, ["foo","data_clinical_supp_sample_SAGE.txt"])
	assert clin.validateFilename(["data_clinical_supp_SAGE.txt"]) == "clinical"
	assert clin.validateFilename(["data_clinical_supp_sample_SAGE.txt","data_clinical_supp_patient_SAGE.txt"]) == "clinical"


	patientDf = pd.DataFrame(dict(PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								 SEX=[1,2,1,2,float('nan')],
								 PRIMARY_RACE=[1,2,3,4,float('nan')],
								 SECONDARY_RACE=[1,2,3,4,float('nan')],
								 TERTIARY_RACE=[1,2,3,4,float('nan')],
								 ETHNICITY=[1,2,3,4,float('nan')],
								 BIRTH_YEAR=[float('nan'),1990,1990,1990,1990],
								 CENTER=["FOO","FOO","FOO","FOO","FOO"]))

	sampleDf = pd.DataFrame(dict(SAMPLE_ID=["ID1-1","ID2-1","ID3-1","ID4-1","ID5-1"],
								 PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								 AGE_AT_SEQ_REPORT=[100000,100000,100000,100000,100000],
								 ONCOTREE_CODE=['AMPCA','AMPCA','AMPCA','AMPCA','AMPCA'],
								 SAMPLE_TYPE=[1,2,3,4,4],
								 SEQ_ASSAY_ID=['SAGE-1','SAGE-1','SAGE-1','SAGE-1','SAGE-1'],
								 SEQ_DATE=['Jan-2013','ApR-2013','Jul-2013','Oct-2013','release']))

	error, warning = clin.validate_helper(patientDf, sampleDf, oncotree_url)
	assert error == ""
	assert warning == ""
	error, warning = clin.validate_helper(patientDf, sampleDf, json_oncotreeurl)
	assert error == ""
	assert warning == ""
	


	#TEST MISSING COLUMNS
	sampleDf = pd.DataFrame()
	patientDf = pd.DataFrame()
	error, warning = clin.validate_helper(patientDf, sampleDf, json_oncotreeurl)
	expectedErrors = ("Sample: clinical file must have SAMPLE_ID column.\n"
					  "Sample: clinical file must have AGE_AT_SEQ_REPORT column.\n"
					  "Sample: clinical file must have ONCOTREE_CODE column.\n"
					  "Sample: clinical file must have SAMPLE_TYPE column.\n"
					  "Sample: clinical file must have SEQ_ASSAY_ID column.\n"
					  "Sample: clinical file must SEQ_DATE column\n"
					  "Patient: clinical file must have BIRTH_YEAR column.\n"
					  "Patient: clinical file must have PATIENT_ID column.\n"
					  "Sample: clinical file must have PATIENT_ID column.\n"
					  "Patient: clinical file must have SEX column.\n")

	expectedWarnings = ("Patient: clinical file doesn't have PRIMARY_RACE column. A blank column will be added\n"
						"Patient: clinical file doesn't have SECONDARY_RACE column. A blank column will be added\n"
						"Patient: clinical file doesn't have TERTIARY_RACE column. A blank column will be added\n"
						"Patient: clinical file doesn't have ETHNICITY column. A blank column will be added\n")
	assert error == expectedErrors
	assert warning == expectedWarnings

	#TEST 
	sampleDf = pd.DataFrame(dict(SAMPLE_ID=[float('nan'),"ID2-1","ID3-1","ID4-1","ID5-1"],
								 PATIENT_ID=[float('nan'),"ID2","ID3","ID7","ID5"],
								 AGE_AT_SEQ_REPORT=[10,100000,100000,100000,100000],
								 ONCOTREE_CODE=['AMPCAD','TESTIS','AMPCA','AMPCA','UCEC'],
								 SAMPLE_TYPE=[1,2,3,4,float('nan')],
								 SEQ_ASSAY_ID=[float('nan'),'Sage-1','SAGE-1','SAGE-1','SAGE-1'],
								 SEQ_DATE=['Jane-2013','Jan-2013','Jan-2013','Jan-2013','Jan-2013']))

	patientDf = pd.DataFrame(dict(PATIENT_ID=["ID6","ID2","ID3",float("nan"),"ID5"],
								 SEX=[1,2,1,5,float('nan')],
								 PRIMARY_RACE=[1,2,3,6,float('nan')],
								 SECONDARY_RACE=[1,2,3,6,float('nan')],
								 TERTIARY_RACE=[1,2,3,6,float('nan')],
								 ETHNICITY=[1,2,3,6,float('nan')],
								 BIRTH_YEAR=[1990,1990,1990,1990,1990],
								 CENTER=["FOO","FOO","FOO","FOO","FOO"]))

	error, warning = clin.validate_helper(patientDf, sampleDf, json_oncotreeurl)
	expectedErrors = ("Sample: There can't be any blank values for SAMPLE_ID\n"
					  "Sample: Please double check that all your ONCOTREE CODES exist in the mapping. You have 1 samples that don't map. These are the codes that don't map: AMPCAD\n"
					  "Sample: Please double check your SAMPLE_TYPE column. No null values allowed.\n"
					  "Sample: Please normalize your SEQ_ASSAY_ID names.  You have these SEQ_ASSAY_IDs: Sage-1, SAGE-1.\n"
					  "Sample: SEQ_DATE must be one of five values- For Jan-March: use Jan-YEAR. For Apr-June: use Apr-YEAR. For July-Sep: use Jul-YEAR. For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) For values that don't have SEQ_DATES that you want released use 'release'.\n"
					  "Patient: There can't be any blank values for PATIENT_ID\n"
					  "Sample: There can't be any blank values for PATIENT_ID\n"
					  "Sample: PATIENT_ID's much be contained in the SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n"
					  "Sample: All samples must have associated patient information. These samples are missing patient data: ID4-1\n"
					  "Patient: Please double check your PRIMARY_RACE column.  This column must be these values 1, 2, 3, 4, or blank.\n"
					  "Patient: Please double check your SECONDARY_RACE column.  This column must be these values 1, 2, 3, 4, or blank.\n"
					  "Patient: Please double check your TERTIARY_RACE column.  This column must be these values 1, 2, 3, 4, or blank.\n"
					  "Patient: Please double check your SEX column.  This column must be these values 1, 2, or blank.\n"
					  "Patient: Please double check your ETHNICITY column.  This column must be these values 1, 2, 3, 4, or blank.\n")
	expectedWarnings = ("Sample: Some SAMPLE_IDs have conflicting SEX and ONCOTREE_CODES: ID2-1,ID5-1\n"
						"Sample: Please double check your SEQ_ASSAY_ID columns, there are empty rows.\n"
						"Sample: All patients must have associated sample information. These patients are missing sample data: ID6\n")
	assert error == expectedErrors
	assert warning == expectedWarnings