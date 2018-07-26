import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR,"../../processing"))

from sampleRetraction import sampleRetraction
from patientRetraction import patientRetraction

def test_processing():
	syn = mock.create_autospec(synapseclient.Synapse) 

	sr = sampleRetraction(syn, "SAGE")

	expectedsrDf = pd.DataFrame(dict(genieSampleId=["GENIE-SAGE-ID1-1","GENIE-SAGE-ID2-1","GENIE-SAGE-ID3-1","GENIE-SAGE-ID4-1","GENIE-SAGE-ID5-1"],
									 retractionDate=[1523039400000,1523039400000,1523039400000,1523039400000,1523039400000],
									 center=["SAGE","SAGE","SAGE","SAGE","SAGE"]))

	srDf = pd.DataFrame({0:["ID1-1","ID2-1","ID3-1","ID4-1","ID5-1"]})
	
	newsrDf = sr._process(srDf,"2018-04-06T18:30:00")
	assert expectedsrDf.equals(newsrDf[expectedsrDf.columns])

	pr = patientRetraction(syn, "SAGE")

	expectedprDf = pd.DataFrame(dict(geniePatientId=["GENIE-SAGE-ID1","GENIE-SAGE-ID2","GENIE-SAGE-ID3","GENIE-SAGE-ID4","GENIE-SAGE-ID5"],
									 retractionDate=[1523125800000,1523125800000,1523125800000,1523125800000,1523125800000],
									 center=["SAGE","SAGE","SAGE","SAGE","SAGE"]))

	prDf = pd.DataFrame({0:["ID1","ID2","ID3","ID4","ID5"]})

	newprDf = pr._process(prDf,"2018-04-07T18:30:00")
	assert expectedprDf.equals(newprDf[expectedprDf.columns])


def test_validation():
	syn = mock.create_autospec(synapseclient.Synapse) 

	sr = sampleRetraction(syn, "SAGE")
	assert_raises(AssertionError, sr.validateFilename, ["foo"])
	assert sr.validateFilename(["sampleRetraction.csv"]) == "sampleRetraction"

	pr = patientRetraction(syn, "SAGE")
	assert_raises(AssertionError, pr.validateFilename, ["foo"])
	assert pr.validateFilename(["patientRetraction.csv"]) == "patientRetraction"
