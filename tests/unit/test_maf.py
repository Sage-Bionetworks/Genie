import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR,"../../processing_sage"))

from maf import maf
from mafSP import mafSP

def test_processing():

	syn = mock.create_autospec(synapseclient.Synapse) 

	mafClass = maf(syn, "SAGE")
	mafDf = pd.DataFrame(dict(Center=["foo","dsdf","sdf"],
							 Tumor_Sample_Barcode=["GENIE-SAGE-1-3","1-2","3-2"],
							 Sequence_Source=["3","e","sd"],
							 Sequencer=["dsf","sdf","d"],
							 Validation_Status=["Unknown","unknown","f"]))
	formattedMafDf = mafClass.formatMAF(mafDf)

	expectedMafDf = pd.DataFrame(dict(Center=["SAGE","SAGE","SAGE"],
							 Tumor_Sample_Barcode=["GENIE-SAGE-1-3","GENIE-SAGE-1-2","GENIE-SAGE-3-2"],
							 Sequence_Source=[pd.np.nan,pd.np.nan,pd.np.nan],
							 Sequencer=[pd.np.nan,pd.np.nan,pd.np.nan],
							 Validation_Status=['','',"f"]))
	assert expectedMafDf.equals(formattedMafDf[expectedMafDf.columns])

def test_validation():

	syn = mock.create_autospec(synapseclient.Synapse) 

	mafClass = maf(syn, "SAGE")
	mafSPClass = mafSP(syn, "SAGE")

	assert_raises(AssertionError, mafClass.validateFilename, ["foo"])
	assert mafClass.validateFilename(["data_mutations_extended_SAGE.txt"]) == "maf"


	correct_column_headers = ['CHROMOSOME','START_POSITION','REFERENCE_ALLELE','TUMOR_SAMPLE_BARCODE','T_ALT_COUNT'] #T_REF_COUNT + T_ALT_COUNT = T_DEPTH
	optional_headers = ['T_REF_COUNT','N_DEPTH','N_REF_COUNT','N_ALT_COUNT']
	tumors = ['TUMOR_SEQ_ALLELE2','TUMOR_SEQ_ALLELE1']

	mafDf = pd.DataFrame(dict(CHROMOSOME=[1,2,3,4,5],
							 START_POSITION=[1,2,3,4,2],
							 REFERENCE_ALLELE=["A","A","A","A","A"],
							 TUMOR_SAMPLE_BARCODE=["ID1-1","ID1-1","ID1-1","ID1-1","ID1-1"],
							 T_ALT_COUNT=[1,2,3,4,3],
							 T_DEPTH=[1,2,3,4,3],
							 T_REF_COUNT=[1,2,3,4,3],
							 N_DEPTH=[1,2,3,4,3],
							 N_REF_COUNT=[1,2,3,4,3],
							 N_ALT_COUNT=[1,2,3,4,3],
							 TUMOR_SEQ_ALLELE1=["A","A","A","A","A"],
							 TUMOR_SEQ_ALLELE2=[float('nan'),float('nan'),float('nan'),float('nan'),float('nan')]))

	error, warning = mafClass.validate_helper(mafDf)
	assert error == ""
	assert warning == ""

	mafDf = pd.DataFrame(dict(CHROMOSOME=[1,2,3,4,5],
							 START_POSITION=[1,2,3,4,2],
							 REFERENCE_ALLELE=["NA",float('nan'),"A","A","A"],
							 TUMOR_SAMPLE_BARCODE=["ID1-1","ID1-1","ID1-1","ID1-1","ID1-1"],
							 T_ALT_COUNT=[1,2,3,4,3],
							 T_DEPTH=[1,2,3,4,3],
							 N_REF_COUNT=[1,2,3,4,3],
							 N_ALT_COUNT=[1,2,3,4,3]))	

	error, warning = mafClass.validate_helper(mafDf)
	expectedErrors = ("Your mutation file must also have at least one of these headers: TUMOR_SEQ_ALLELE2 or TUMOR_SEQ_ALLELE1.\n"
					  "Your mutation file cannot have any empty REFERENCE_ALLELE values.\n")
	expectedWarnings = ("Your mutation file does not have the column headers that can give extra information to the processed mutation file: T_REF_COUNT, N_DEPTH.\n"
						"Your REFERENCE_ALLELE column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n")
	assert error == expectedErrors
	assert warning == expectedWarnings


	mafDf = pd.DataFrame(dict(START_POSITION=[1,2,3,4,2],
							 REFERENCE_ALLELE=["A","A","A","A","A"],
							 TUMOR_SAMPLE_BARCODE=["ID1-1","ID1-1","ID1-1","ID1-1","ID1-1"],
							 T_ALT_COUNT=[1,2,3,4,3],
							 N_DEPTH=[1,2,3,4,3],
							 N_REF_COUNT=[1,2,3,4,3],
							 N_ALT_COUNT=[1,2,3,4,3],
							 TUMOR_SEQ_ALLELE1=["NA",float('nan'),"A","A","A"]))

	error, warning = mafClass.validate_helper(mafDf)
	expectedErrors = ("First column header must be one of these: CHROMOSOME, HUGO_SYMBOL, TUMOR_SAMPLE_BARCODE.\n"
					  "If you are missing T_DEPTH, you must have T_REF_COUNT!\n"
					  "Your mutation file must at least have these headers: CHROMOSOME.\n"
					  "Your mutation file must have at least one of these: TUMOR_SEQ_ALLELE2, TUMOR_SEQ_ALLELE1 and at least one of those columns cannot have any empty values.\n")
	expectedWarnings = ("Your TUMOR_SEQ_ALLELE1 column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n"
						"Your mutation file does not have the column headers that can give extra information to the processed mutation file: T_REF_COUNT.\n")
	assert error == expectedErrors
	assert warning == expectedWarnings

	mafDf = pd.DataFrame(dict(START_POSITION=[1,2,3,4,2],
							 REFERENCE_ALLELE=["A","A","A","A","A"],
							 TUMOR_SAMPLE_BARCODE=["ID1-1","ID1-1","ID1-1","ID1-1","ID1-1"],
							 N_DEPTH=[1,2,3,4,3],
							 N_REF_COUNT=[1,2,3,4,3],
							 N_ALT_COUNT=[1,2,3,4,3],
							 TUMOR_SEQ_ALLELE2=["NA",float('nan'),"A","A","A"]))
	error, warning = mafSPClass.validate_helper(mafDf, SP=True)
	expectedErrors = ("First column header must be one of these: CHROMOSOME, HUGO_SYMBOL, TUMOR_SAMPLE_BARCODE.\n"
					  "Your mutation file must at least have these headers: CHROMOSOME.\n"
					  "Your mutation file must have at least one of these: TUMOR_SEQ_ALLELE2, TUMOR_SEQ_ALLELE1 and at least one of those columns cannot have any empty values.\n")
	expectedWarnings = ("Your TUMOR_SEQ_ALLELE2 column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n")
	assert error == expectedErrors
	assert warning == expectedWarnings