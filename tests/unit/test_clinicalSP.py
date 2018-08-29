import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys
from genie.clinicalSP import clinicalSP


def test_processing():

	syn = mock.create_autospec(synapseclient.Synapse) 


	clin = clinicalSP(syn, "SAGE")

	expectedclinDf = pd.DataFrame(dict(PATIENT_ID=["GENIE-SAGE-ID1","GENIE-SAGE-ID2","GENIE-SAGE-ID3","GENIE-SAGE-ID4","GENIE-SAGE-ID5"],
								 SAMPLE_ID=["GENIE-SAGE-ID1-1","GENIE-SAGE-ID2-1","GENIE-SAGE-ID3-1","GENIE-SAGE-ID4-1","GENIE-SAGE-ID5-1"],
								 CENTER = ["SAGE","SAGE","SAGE","SAGE","SAGE"],
								 SEQ_ASSAY_ID=["SAGE-TEST","SAGE-TEST","SAGE-TEST","SAGE-TEST","SAGE-TEST"]))
	#TEST patient processing
	clinDf = pd.DataFrame(dict(PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								 SAMPLE_ID=["ID1-1","ID2-1","ID3-1","ID4-1","ID5-1"],
								 SEQ_ASSAY_ID=["SAGE-TEST","SAGE-TEST","SAGE-TEST","SAGE-TEST","SAGE-TEST"]))

	clinDf = clin._process(clinDf)
	assert expectedclinDf.equals(clinDf[expectedclinDf.columns])

def test_validation():

	syn = mock.create_autospec(synapseclient.Synapse) 

	clin = clinicalSP(syn, "SAGE")
	assert_raises(AssertionError, clin.validateFilename, ["foo"])
	assert clin.validateFilename(["nonGENIE_data_clinical.txt"]) == "clinicalSP"


	clinDf = pd.DataFrame(dict(PATIENT_ID=["ID1","ID2","ID3","ID4","ID5"],
								 SAMPLE_ID=["ID1-1","ID2-1","ID3-1","ID4-1","ID5-1"],
								 SEQ_ASSAY_ID=["SAGE-TEST","SAGE-TEST","SAGE-TEST","SAGE-TEST","SAGE-TEST"]))

	error, warning = clin._validate(clinDf)
	assert error == ""
	assert warning == ""

	clinDf = pd.DataFrame()

	error, warning = clin._validate(clinDf)
	expectedErrors = ("nonGENIE_data_clinical.txt: File must have SAMPLE_ID column.\n"
					  "nonGENIE_data_clinical.txt: File must have SEQ_ASSAY_ID column.\n"
					  "nonGENIE_data_clinical.txt: File must have PATIENT_ID column.\n")
	assert error == expectedErrors
	assert warning == ""

	clinDf = pd.DataFrame(dict(PATIENT_ID=[float('nan'),"ID2","ID3","ID4","ID5"],
								 SAMPLE_ID=["ID1-1",float('nan'),"ID3-1","ID4-1","ID5-1"],
								 SEQ_ASSAY_ID=[float('nan'),"SAGE-TEST","SAGE-TEST","SAGE-TEST","SAGE-TEST"]))

	error, warning = clin._validate(clinDf)

	expectedErrors = ("nonGENIE_data_clinical.txt: There can't be any blank values for SAMPLE_ID\n"
					  "nonGENIE_data_clinical.txt: There can't be any blank values for PATIENT_ID\n")
	expectedWarnings = ("nonGENIE_data_clinical.txt: Please double check your SEQ_ASSAY_ID columns, there are empty rows.\n")
	
	assert error == expectedErrors
	assert warning == expectedWarnings
