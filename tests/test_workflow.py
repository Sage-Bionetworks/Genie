from unittest import mock

import pytest
import synapseclient

from genie_registry.workflow import workflow

syn = mock.create_autospec(synapseclient.Synapse)
workflow_class = workflow(syn, "SAGE")


def test_processing():
    pass


@pytest.fixture(params=[(["foo"]), (["SAGE-test.txt"])])
def filename_fileformat_map(request):
    return request.param


def test_incorrect_validatefilename(filename_fileformat_map):
    filepath_list = filename_fileformat_map
    with pytest.raises(AssertionError):
        workflow_class.validateFilename(filepath_list)


def test_correct_validatefilename():
    assert workflow_class.validateFilename(["SAGE-test.md"]) == "md"
