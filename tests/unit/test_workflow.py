import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR,"../../processing"))

from workflow import workflow

def test_processing():

	syn = mock.create_autospec(synapseclient.Synapse) 

	workflowClass = workflow(syn, "SAGE")
	pass

def test_validation():

	syn = mock.create_autospec(synapseclient.Synapse) 

	workflowClass = workflow(syn, "SAGE")

	assert_raises(AssertionError, workflowClass.validateFilename, ["foo"])
	assert_raises(AssertionError, workflowClass.validateFilename, ["SAGE-test.txt"])
	assert workflowClass.validateFilename(["SAGE-test.md"]) == "md"
