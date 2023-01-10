import pytest

from genie_registry.workflow import workflow


@pytest.fixture
def workflow_class(syn):
    return workflow(syn, "SAGE")


def test_processing():
    pass


@pytest.fixture(params=[(["foo"]), (["SAGE-test.txt"])])
def filename_fileformat_map(request):
    return request.param


def test_incorrect_validatefilename(workflow_class, filename_fileformat_map):
    with pytest.raises(AssertionError):
        workflow_class.validateFilename(filename_fileformat_map)


def test_correct_validatefilename(workflow_class):
    assert workflow_class.validateFilename(["SAGE-test.md"]) == "md"
