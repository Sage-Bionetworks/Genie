import synapseclient
import pandas as pd
import mock
from genie.clinical import clinical
import pytest

def createMockTable(dataframe):
	table = mock.create_autospec(synapseclient.table.CsvFileTable)
	table.asDataFrame.return_value= dataframe
	return(table)

def table_query_results(*args):
	return(table_query_results_map[args])

with_nan = pd.DataFrame(dict(CODE=[1,2,3,4,float('nan')],
							CBIO_LABEL=['Test','Why','foo','Me','Unknown'],
							DESCRIPTION=['non','asdf','asdf','asdff','asdfasdf']))

no_nan = pd.DataFrame(dict(CODE=[1,2,3,4],
						  CBIO_LABEL=['Test','Why','foo','Me'],
						  DESCRIPTION=['non','asdf','asdf','asdff']))

sexdf = pd.DataFrame(dict(CODE=[1,2,float('nan')],
						  CBIO_LABEL=['Male','Female','Unknown'],
						  DESCRIPTION=['Male','Female','Not coded']))

table_query_results_map = {
("SELECT * FROM syn7434222",) : createMockTable(sexdf),
("SELECT * FROM syn7434236",) : createMockTable(with_nan),
("SELECT * FROM syn7434242",) : createMockTable(with_nan),
("SELECT * FROM syn7434273",) : createMockTable(no_nan)
}   

syn = mock.create_autospec(synapseclient.Synapse) 
syn.tableQuery.side_effect=table_query_results
clin_class = clinical(syn, "SAGE")

oncotree_url = 'http://oncotree.mskcc.org/api/tumor_types.txt?version=oncotree_latest_stable'
json_oncotreeurl = "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"

def test_filetype():
	assert clin_class._fileType == "clinical"

def test_incorrect_validatefilename():
	with pytest.raises(AssertionError):
		clin_class.validateFilename(['foo'])
		clin_class.validateFilename(["foo","data_clinical_supp_sample_SAGE.txt"])

def test_correct_validatefilename():
	assert clin_class.validateFilename(["data_clinical_supp_SAGE.txt"]) == "clinical"
	assert clin_class.validateFilename(["data_clinical_supp_sample_SAGE.txt","data_clinical_supp_patient_SAGE.txt"]) == "clinical"

def test_patient__process():
	expected_patientdf = pd.DataFrame(dict(PATIENT_ID=["GENIE-SAGE-ID1","GENIE-SAGE-ID2","GENIE-SAGE-ID3","GENIE-SAGE-ID4","GENIE-SAGE-ID5"],
									   SEX=['Male','Female','Male','Female','Unknown'],
									   PRIMARY_RACE=['Test','Why','foo','Me','Unknown'],
									   SECONDARY_RACE=['Test','Why','foo','Me','Unknown'],
									   TERTIARY_RACE=['Test','Why','foo','Me','Unknown'],
									   ETHNICITY=['Test','Why','foo','Me','Unknown'],
									   BIRTH_YEAR=[1990,1990,1990,1990,1990],
									   CENTER=["SAGE","SAGE","SAGE","SAGE","SAGE"]))
	#TEST patient processing
	#Clinical file headers are capitalized prior to processing
	patientdf = pd.DataFrame(dict(PATIENT_Id=["ID1","ID2","ID3","ID4","ID5"],
							   sex=[1,2,1,2,float('nan')],
							   PRIMARY_RACE=[1,2,3,4,float('nan')],
							   Secondary_RACE=[1,2,3,4,float('nan')],
							   TERTIARY_RACE=[1,2,3,4,float('nan')],
							   ETHNICITY=[1,2,3,4,float('nan')],
							   BIRTH_YEAR=[1990,1990,1990,1990,1990],
							   CENTER=["FOO","FOO","FOO","FOO","FOO"]))
	patient_cols = ["PATIENT_ID","SEX","PRIMARY_RACE","SECONDARY_RACE",
				"TERTIARY_RACE","ETHNICITY","BIRTH_YEAR","CENTER"]
	clinical_template = pd.DataFrame(columns=patient_cols)
	patient=True
	new_patientdf = clin_class._process(patientdf, clinical_template)
	assert expected_patientdf.equals(new_patientdf[expected_patientdf.columns])

def test_sample__process():
	#TEST for sample processing
	expected_sampledf = pd.DataFrame(dict(SAMPLE_ID=["GENIE-SAGE-ID1-1","GENIE-SAGE-ID2-1","GENIE-SAGE-ID3-1","GENIE-SAGE-ID4-1","GENIE-SAGE-ID5-1"],
									   PATIENT_ID=["GENIE-SAGE-ID1","GENIE-SAGE-ID2","GENIE-SAGE-ID3","GENIE-SAGE-ID4","GENIE-SAGE-ID5"],
									   AGE_AT_SEQ_REPORT=[100000,100000,100000,100000,100000],
									   ONCOTREE_CODE=['AMPCA','AMPCA','AMPCA','AMPCA','AMPCA'],
									   SAMPLE_TYPE=['Test','Why','foo','Me','Me'],
									   CENTER=["SAGE","SAGE","SAGE","SAGE","SAGE"],
									   SAMPLE_TYPE_DETAILED=['non','asdf','asdf','asdff','asdff'],
									   SEQ_ASSAY_ID=['SAGE-1','SAGE-1','SAGE-1','SAGE-1','SAGE-1'],
									   SEQ_DATE=['Jan-2012','Apr-2013','Jul-2014','Oct-2015','Release'],
									   SEQ_YEAR=[2012,2013,2014,2015,pd.np.nan]))
	sample_cols = ["SAMPLE_ID","PATIENT_ID","AGE_AT_SEQ_REPORT","ONCOTREE_CODE","SAMPLE_TYPE",
			   "SEQ_ASSAY_ID",'SEQ_DATE','SAMPLE_TYPE_DETAILED','SEQ_YEAR']

	clinical_template = pd.DataFrame(columns=sample_cols)
	patient=False
	sampledf = pd.DataFrame(dict(SAMPLE_ID=["ID1-1","ID2-1","ID3-1","ID4-1","ID5-1"],
							   PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
							   Age_AT_SEQ_REPORT=[100000,100000,100000,100000,100000],
							   ONCOTree_CODE=['AMPCA','AMPCA','AMPCA','AMPCA','AMPCA'],
							   SAMPLE_TYPE=[1,2,3,4,4],
							   SEQ_ASSAY_ID=['SAGE-1','SAGE-1','SAGE-1','SAGE-1','SAGE-1'],
							   SEQ_DATE=['Jan-2012','Apr-2013','JUL-2014','Oct-2015','release']))

	new_sampledf = clin_class._process(sampledf, clinical_template)
	assert expected_sampledf.equals(new_sampledf[expected_sampledf.columns])

def test_perfect__validate():
	patientdf = pd.DataFrame(dict(PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								  SEX=[1,2,1,2,float('nan')],
								  PRIMARY_RACE=[1,2,3,4,float('nan')],
								  SECONDARY_RACE=[1,2,3,4,float('nan')],
								  TERTIARY_RACE=[1,2,3,4,float('nan')],
								  ETHNICITY=[1,2,3,4,float('nan')],
								  BIRTH_YEAR=[float('nan'),1990,1990,1990,1990],
								  CENTER=["FOO","FOO","FOO","FOO","FOO"]))

	sampledf = pd.DataFrame(dict(SAMPLE_ID=["ID1-1","ID2-1","ID3-1","ID4-1","ID5-1"],
								 PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								 AGE_AT_SEQ_REPORT=[100000,100000,100000,100000,100000],
								 ONCOTREE_CODE=['AMPCA','AMPCA','AMPCA','AMPCA','AMPCA'],
								 SAMPLE_TYPE=[1,2,3,4,4],
								 SEQ_ASSAY_ID=['SAGE-1-1','SAGE-SAGE-1','SAGE-1','SAGE-1','SAGE-1'],
								 SEQ_DATE=['Jan-2013','ApR-2013','Jul-2013','Oct-2013','release']))

	clinicaldf = patientdf.merge(sampledf, on="PATIENT_ID")
	error, warning = clin_class._validate(clinicaldf, oncotree_url)
	assert error == ""
	assert warning == ""
	error, warning = clin_class._validate(clinicaldf, json_oncotreeurl)
	assert error == ""
	assert warning == ""

def test_missingcols__validate():
	clincaldf = pd.DataFrame()
	error, warning = clin_class._validate(clincaldf, json_oncotreeurl)
	expected_errors = ("Sample: clinical file must have SAMPLE_ID column.\n"
					"Sample/Clinical: clinical files must have PATIENT_ID column.\n"
					"Sample: clinical file must have AGE_AT_SEQ_REPORT column.\n"
					"Sample: clinical file must have ONCOTREE_CODE column.\n"
					"Sample: clinical file must have SAMPLE_TYPE column.\n"
					"Sample: clinical file must have SEQ_ASSAY_ID column.\n"
					"Sample: clinical file must SEQ_DATE column\n"
					"Patient: clinical file must have BIRTH_YEAR column.\n"
					"Patient: clinical file must have SEX column.\n")

	expected_warnings = ("Patient: clinical file doesn't have PRIMARY_RACE column. A blank column will be added\n"
					  "Patient: clinical file doesn't have SECONDARY_RACE column. A blank column will be added\n"
					  "Patient: clinical file doesn't have TERTIARY_RACE column. A blank column will be added\n"
					  "Patient: clinical file doesn't have ETHNICITY column. A blank column will be added\n")
	assert error == expected_errors
	assert warning == expected_warnings

def test_errors__validate():
	sampledf = pd.DataFrame(dict(SAMPLE_ID=[float('nan'),"ID2-1","ID3-1","ID4-1","ID5-1"],
	                           PATIENT_ID=["ID6","ID2","ID3",float('nan'),"ID5"],
	                           AGE_AT_SEQ_REPORT=[10,100000,100000,100000,100000],
	                           ONCOTREE_CODE=['AMPCAD','TESTIS','AMPCA','AMPCA','UCEC'],
	                           SAMPLE_TYPE=[1,2,3,4,float('nan')],
	                           SEQ_ASSAY_ID=[float('nan'),'Sage-1','SAGE-1','S-SAGE-1','SAGE-1'],
	                           SEQ_DATE=['Jane-2013','Jan-2013','Jan-2013','Jan-2013','Jan-2013']))
	patientdf = pd.DataFrame(dict(PATIENT_ID=["ID6","ID2","ID3",float("nan"),"ID5"],
	                           SEX=[1,2,1,5,float('nan')],
	                           PRIMARY_RACE=[1,2,3,6,float('nan')],
	                           SECONDARY_RACE=[1,2,3,6,float('nan')],
	                           TERTIARY_RACE=[1,2,3,6,float('nan')],
	                           ETHNICITY=[1,2,3,6,float('nan')],
	                           BIRTH_YEAR=[1990,1990,1990,1990,1990],
	                           CENTER=["FOO","FOO","FOO","FOO","FOO"]))
	clinicaldf = patientdf.merge(sampledf, on="PATIENT_ID")

	error, warning = clin_class._validate(clinicaldf, json_oncotreeurl)
	expected_errors = ("Sample: PATIENT_ID's much be contained in the SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n"
	                "Patient: All samples must have associated patient information and no null patient ids allowed. These samples are missing patient data: ID4-1\n"
	                "Sample: Please double check that all your ONCOTREE CODES exist in the mapping. You have 1 samples that don't map. These are the codes that don't map: AMPCAD\n"
	                "Sample: Please double check your SAMPLE_TYPE column. No null values allowed.\n"
	                "Sample: Please make sure your SEQ_ASSAY_IDs start with your center abbreviation: S-SAGE-1.\n"
	                "Sample: SEQ_DATE must be one of five values- For Jan-March: use Jan-YEAR. For Apr-June: use Apr-YEAR. For July-Sep: use Jul-YEAR. For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) For values that don't have SEQ_DATES that you want released use 'release'.\n"
	                "Patient: Please double check your PRIMARY_RACE column.  This column must be these values 1, 2, 3, 4, or blank.\n"
	                "Patient: Please double check your SECONDARY_RACE column.  This column must be these values 1, 2, 3, 4, or blank.\n"
	                "Patient: Please double check your TERTIARY_RACE column.  This column must be these values 1, 2, 3, 4, or blank.\n"
	                "Patient: Please double check your SEX column.  This column must be these values 1, 2, or blank.\n"
	                "Patient: Please double check your ETHNICITY column.  This column must be these values 1, 2, 3, 4, or blank.\n")
	expected_warnings = ("Sample: All patients must have associated sample information. These patients are missing sample data: ID6\n"
	                  "Sample: Some SAMPLE_IDs have conflicting SEX and ONCOTREE_CODES: ID2-1,ID5-1\n"
	                  "Sample: Please double check your SEQ_ASSAY_ID columns, there are empty rows.\n")
	assert error == expected_errors
	assert warning == expected_warnings
