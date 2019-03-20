import synapseclient
import mock
import pytest
from genie.workflow import workflow

syn = mock.create_autospec(synapseclient.Synapse)
workflowClass = workflow(syn, "SAGE")


def test_processing():
    pass


def test_validation():
    with pytest.raises(AssertionError):
        workflowClass.validateFilename(["foo"])
        workflowClass.validateFilename(["SAGE-test.txt"])
    assert workflowClass.validateFilename(["SAGE-test.md"]) == "md"
