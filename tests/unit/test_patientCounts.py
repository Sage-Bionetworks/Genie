import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys

from genie.patientCounts import patientCounts


def test_processing():

    #def _process(self, patientCountsDf):
	oncotree_url = 'http://oncotree.mskcc.org/api/tumor_types.txt?version=oncotree_latest_stable'
	json_oncotreeurl = "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"

	syn = mock.create_autospec(synapseclient.Synapse) 

	pc = patientCounts(syn, "SAGE")

	expectedpcDf = pd.DataFrame(dict(CENTER=['SAGE','SAGE','SAGE','SAGE','SAGE'],
									 ONCOTREE_CODE=['LUAD','AMPCA','CHOL','IHCH','UCA'],
									 NUM_PATIENTS_PD1_PDL1=[1,2,3,4,3],
									 PRIMARY_CODE=['LUNG','AMPULLA_OF_VATER','BILIARY_TRACT','BILIARY_TRACT','BLADDER']))

	pcDf = pd.DataFrame(dict(CENTER=['foo','foo','foo','foo','foo'],
							 ONCOTREE_CODE=['LUAD','AMPCA','CHOL','IHCH','UCA'],
							 NUM_PATIENTS_PD1_PDL1=[1,2,3,4,3]))
	
	newpcDf = pc._process(pcDf, oncotree_url)
	assert expectedpcDf.equals(newpcDf[expectedpcDf.columns])



def test_validation():
	syn = mock.create_autospec(synapseclient.Synapse) 

	pc = patientCounts(syn, "SAGE")

	assert_raises(AssertionError, pc.validateFilename, ["foo"])
	assert pc.validateFilename(["patient_counts.txt"]) == "patientCounts"

	pcDf = pd.DataFrame(dict(CENTER=['foo','foo','foo','foo','foo'],
							 ONCOTREE_CODE=['LUAD','AMPCA','CHOL','IHCH','UCA'],
							 NUM_PATIENTS_PD1_PDL1=[1,2,3,4,3]))
	
	oncotree_url = 'http://oncotree.mskcc.org/api/tumor_types.txt?version=oncotree_latest_stable'
	json_oncotreeurl = "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"

	error, warning = pc._validate(pcDf, oncotree_url)
	assert error == ""
	assert warning == ""
	error, warning = pc._validate(pcDf, json_oncotreeurl)
	assert error == ""
	assert warning == ""

	pcDf = pd.DataFrame(dict(ONCOTREE_CODE=['LUAD','LUAD','CHOLD','IHCH','UCA'],
							 NUM_PATIENTS_PD1_PDL1=[1,2,3,4,float('nan')]))

	error, warning = pc._validate(pcDf, json_oncotreeurl)
	expectedErrors = ("Patient Counts: Must not have any duplicated ONCOTREE CODES.\n"
					  "Patient Counts: Please double check that all your ONCOTREE CODES exist in the mapping. You have 1 codes that don't map. These are the codes that don't map: CHOLD\n"
					  "Patient Counts: Must not have any null values, and must be all integers.\n")
	assert error == expectedErrors
	assert warning == ""

	pcDf = pd.DataFrame(dict(ONCOTREE_CODE=['LUAD','AMPCA','CHOL','IHCH','UCA']))
	error, warning = pc._validate(pcDf, json_oncotreeurl)

	expectedErrors = ("Patient Counts: File must have NUM_PATIENTS_PD1_PDL1 column.\n")

	assert error == expectedErrors
	assert warning == ""

	pcDf = pd.DataFrame(dict(NUM_PATIENTS_PD1_PDL1=[1,2,3,4,5]))
	error, warning = pc._validate(pcDf, json_oncotreeurl)

	expectedErrors = ("Patient Counts: File must have ONCOTREE_CODE column.\n")
	assert error == expectedErrors
	assert warning == ""


#test_validation()
#test_processing()
