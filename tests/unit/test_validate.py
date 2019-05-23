import pytest
import mock
import synapseclient
import pytest
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


def test_perfect_determine_filetype(filename_fileformat_map):
    '''
    Tests determining of file type through filenames
    Parameters are passed in from filename_fileformat_map
    '''
    (filepath_list, fileformat) = filename_fileformat_map
    assert validate.determine_filetype(
        syn, filepath_list, center) == fileformat


def test_wrongfilename_determine_filetype():
    '''
    Tests ValueError is raised when wrong filename is passed
    '''
    with pytest.raises(
        ValueError,
        match="Your filename is incorrect! "
              "Please change your filename before you run "
              "the validator or specify --filetype if you are "
              "running the validator locally"):
        validate.determine_filetype(syn, ['wrong.txt'], center)


def test_valid_determine_validity_and_log():
    '''
    Tests if no error and warning strings are passed that
    returned valid and message is correct
    '''
    valid, message = \
        validate.determine_validity_and_log('', '')
    assert valid
    assert message == "YOUR FILE IS VALIDATED!\n"


def test_invalid_determine_validity_and_log():
    '''
    Tests if error and warnings strings are passed that
    returned valid and message is correct
    '''
    valid, message = \
        validate.determine_validity_and_log("error\nnow", 'warning\nnow')
    assert not valid
    assert message == (
        "----------------ERRORS----------------\n"
        "error\nnow"
        "-------------WARNINGS-------------\n"
        'warning\nnow')


def test_warning_determine_validity_and_log():
    '''
    Tests if no error but warnings strings are passed that
    returned valid and message is correct
    '''
    valid, message = \
        validate.determine_validity_and_log('', 'warning\nnow')
    assert valid
    assert message == (
        "YOUR FILE IS VALIDATED!\n"
        "-------------WARNINGS-------------\n"
        'warning\nnow')


def test_valid_validate_single_file_workflow():
    '''
    Tests that all the functions are run in validate single
    file workflow and all the right things are returned
    '''
    filepathlist = ['clinical.txt']
    error_string = ''
    warning_string = ''
    center = 'SAGE'
    expected_valid = True
    expected_message = "valid message here!"
    expected_filetype = "clinical"
    with mock.patch(
            "genie.validate.determine_filetype",
            return_value=expected_filetype) as mock_determine_filetype,\
        mock.patch(
            "genie.clinical.validate",
            return_value=(error_string, warning_string)) as mock_genie_class,\
        mock.patch(
            "genie.validate.determine_validity_and_log",
            return_value=(expected_valid, expected_message)) as mock_determine:

        valid, message, filetype = validate.validate_single_file_workflow(
            syn,
            filepathlist,
            center)

        assert valid == expected_valid
        assert message == expected_message
        assert filetype == expected_filetype

        mock_determine_filetype.assert_called_once_with(
            syn, filepathlist, center)

        mock_genie_class.assert_called_once_with(
            filePathList=filepathlist,
            oncotreeLink=None,
            testing=False,
            noSymbolCheck=False)

        mock_determine.assert_called_once_with(error_string, warning_string)
