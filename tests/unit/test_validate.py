import pytest
import mock
import synapseclient
from genie import validate
center = "SAGE"
syn = mock.create_autospec(synapseclient.Synapse)


@pytest.fixture(params=[
    # tuple with (input, expectedOutput)
    (["data_CNA_SAGE.txt"], "cna"),
    (["data_clinical_supp_SAGE.txt"], "clinical"),
    (["data_clinical_supp_sample_SAGE.txt",
      "data_clinical_supp_patient_SAGE.txt"], "clinical")])
def filename_fileformat_map(request):
    return request.param


def test_perfect_get_filetype(filename_fileformat_map):
    (filepath_list, fileformat) = filename_fileformat_map
    assert validate.determine_filetype(
        syn, filepath_list, center) == fileformat


# def test_wrongfilename_get_filetype():
#     assert input_to_database.get_filetype(syn, ['wrong.txt'], center) is None
