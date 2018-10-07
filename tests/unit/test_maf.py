import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys

from genie.maf import maf
from genie.mafSP import mafSP

# def test_processing():

# 	syn = mock.create_autospec(synapseclient.Synapse) 

# 	mafClass = maf(syn, "SAGE")
# 	mafDf = pd.DataFrame(dict(Center=["foo","dsdf","sdf"],
# 							 Tumor_Sample_Barcode=["GENIE-SAGE-1-3","1-2","3-2"],
# 							 Sequence_Source=["3","e","sd"],
# 							 Sequencer=["dsf","sdf","d"],
# 							 Validation_Status=["Unknown","unknown","f"]))
# 	formattedMafDf = mafClass.formatMAF(mafDf)

# 	expectedMafDf = pd.DataFrame(dict(Center=["SAGE","SAGE","SAGE"],
# 							 Tumor_Sample_Barcode=["GENIE-SAGE-1-3","GENIE-SAGE-1-2","GENIE-SAGE-3-2"],
# 							 Sequence_Source=[pd.np.nan,pd.np.nan,pd.np.nan],
# 							 Sequencer=[pd.np.nan,pd.np.nan,pd.np.nan],
# 							 Validation_Status=['','',"f"]))
# 	assert expectedMafDf.equals(formattedMafDf[expectedMafDf.columns])

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
							 TUMOR_SEQ_ALLELE2=["A","A","A","A","A"]))

	error, warning = mafClass.validate_helper(mafDf)
	assert error == ""
	assert warning == ""

	mafDf = pd.DataFrame(dict(CHROMOSOME=[1,"chr2","WT",4,5],
							 START_POSITION=[1,2,3,4,2],
							 REFERENCE_ALLELE=["NA",float('nan'),"A","A","A"],
							 TUMOR_SAMPLE_BARCODE=["ID1-1","ID1-1","ID1-1","ID1-1","ID1-1"],
							 T_ALT_COUNT=[1,2,3,4,3],
							 T_DEPTH=[1,2,3,4,3],
							 N_REF_COUNT=[1,2,3,4,3],
							 N_ALT_COUNT=[1,2,3,4,3]))	

	error, warning = mafClass.validate_helper(mafDf)
	expectedErrors = ("Mutation File: Must at least have these headers: TUMOR_SEQ_ALLELE2.\n"
					  "Mutation File: Cannot have any empty REFERENCE_ALLELE values.\n"
					  "Mutation File: CHROMOSOME column cannot have any values that start with 'chr' or any 'WT' values.\n")
	expectedWarnings = ("Mutation File: Does not have the column headers that can give extra information to the processed mutation file: T_REF_COUNT, N_DEPTH.\n"
						"Mutation File: Your REFERENCE_ALLELE column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n")
	assert error == expectedErrors
	assert warning == expectedWarnings


	mafDf = pd.DataFrame(dict(START_POSITION=[1,2,3,4,2],
							 REFERENCE_ALLELE=["A","A","A","A","A"],
							 TUMOR_SAMPLE_BARCODE=["ID1-1","ID1-1","ID1-1","ID1-1","ID1-1"],
							 T_ALT_COUNT=[1,2,3,4,3],
							 N_DEPTH=[1,2,3,4,3],
							 N_REF_COUNT=[1,2,3,4,3],
							 N_ALT_COUNT=[1,2,3,4,3],
							 TUMOR_SEQ_ALLELE2=["NA",float('nan'),"A","A","A"]))

	error, warning = mafClass.validate_helper(mafDf)
	expectedErrors = ("Mutation File: First column header must be one of these: CHROMOSOME, HUGO_SYMBOL, TUMOR_SAMPLE_BARCODE.\n"
					  "Mutation File: If you are missing T_DEPTH, you must have T_REF_COUNT!\n"
					  "Mutation File: Must at least have these headers: CHROMOSOME.\n"
					  "Mutation File: TUMOR_SEQ_ALLELE2 can't have any null values.\n")
	expectedWarnings = ("Mutation File: TUMOR_SEQ_ALLELE2 column contains 'NA' values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n"
						"Mutation File: Does not have the column headers that can give extra information to the processed mutation file: T_REF_COUNT.\n")
	assert error == expectedErrors
	assert warning == expectedWarnings