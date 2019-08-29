import mock
import pytest

import pandas as pd
import synapseclient
from synapseclient.exceptions import SynapseHTTPError

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
    validator = validate.Validator(syn, center, filepath_list)

    assert validator.determine_filetype() == fileformat


def test_wrongfilename_noerror_determine_filetype():
    '''
    Tests None is passed back when wrong filename is passed
    when raise_error flag is False
    '''
    filepathlist = ['wrong.txt']
    validator = validate.Validator(syn, center=center, filepathlist=filepathlist)

    assert validator.file_type is None


def test_valid_collect_errors_and_warnings():
    '''
    Tests if no error and warning strings are passed that
    returned valid and message is correct
    '''
    message = validate.collect_errors_and_warnings('', '')
    assert message == "YOUR FILE IS VALIDATED!\n"


def test_invalid_collect_errors_and_warnings():
    '''
    Tests if error and warnings strings are passed that
    returned valid and message is correct
    '''
    message = validate.collect_errors_and_warnings("error\nnow", 'warning\nnow')
    assert message == (
        "----------------ERRORS----------------\n"
        "error\nnow"
        "-------------WARNINGS-------------\n"
        'warning\nnow')


def test_warning_collect_errors_and_warnings():
    '''
    Tests if no error but warnings strings are passed that
    returned valid and message is correct
    '''
    message = \
        validate.collect_errors_and_warnings('', 'warning\nnow')
    assert message == (
        "YOUR FILE IS VALIDATED!\n"
        "-------------WARNINGS-------------\n"
        'warning\nnow')


def test_valid_validate_single_file():
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
            "genie.validate.Validator.determine_filetype",
            return_value=expected_filetype) as mock_determine_filetype,\
        mock.patch(
            "genie.clinical.clinical.validate",
            return_value=(expected_valid, error_string, warning_string)) as mock_genie_class,\
        mock.patch(
            "genie.validate.collect_errors_and_warnings",
            return_value=expected_message) as mock_determine:

        validator = validate.Validator(syn, center=center, filepathlist=filepathlist)

        valid, message, filetype = validator.validate_single_file()

        assert valid == expected_valid
        assert message == expected_message
        assert filetype == expected_filetype

        mock_determine_filetype.assert_called_once_with()

        mock_genie_class.assert_called_once_with(
            filePathList=filepathlist,
            oncotreeLink=None,
            testing=False,
            noSymbolCheck=False)

        mock_determine.assert_called_once_with(error_string, warning_string)


def test_filetype_validate_single_file():
    '''
    Tests that if filetype is passed in that an error is thrown
    if it is an incorrect filetype
    '''
    filepathlist = ['clinical.txt']
    center = "SAGE"
    expected_error = "----------------ERRORS----------------\nYour filename is incorrect! Please change your filename before you run the validator or specify --filetype if you are running the validator locally"
    validator = validate.Validator(syn, center, filepathlist)

    valid, message, filetype = validator.validate_single_file(filetype="foobar")
    assert message == expected_error


def test_wrongfiletype_validate_single_file():
    '''
    Tests that if there is no filetype for the filename passed
    in, an error is thrown
    '''
    filepathlist = ['clinical.txt']
    center = "SAGE"
    expected_error = '----------------ERRORS----------------\nYour filename is incorrect! Please change your filename before you run the validator or specify --filetype if you are running the validator locally'

    with mock.patch(
            "genie.validate.Validator.determine_filetype",
            return_value=None) as mock_determine_filetype:
        validator = validate.Validator(syn=syn, center=center, 
                                       filepathlist=filepathlist)
        valid, message, filetype = validator.validate_single_file()
        
        assert message == expected_error
        mock_determine_filetype.assert_called_once_with()


def test_nopermission__check_parentid_permission_container():
    '''
    Error thrown if no permissions to access
    '''
    parentid = "syn123"
    with mock.patch.object(syn, "get", side_effect=SynapseHTTPError),\
        pytest.raises(
            ValueError,
            match="Provided Synapse id must be your input folder Synapse id "
                  "or a Synapse Id of a folder inside your input directory"):
        validate._check_parentid_permission_container(syn, parentid)


def test_notcontainer__check_parentid_permission_container():
    '''
    If input if synid of file, throw error
    '''
    parentid = "syn123"
    file_ent = synapseclient.File("foo", parentId=parentid)
    with mock.patch.object(syn, "get", return_value=file_ent),\
        pytest.raises(
            ValueError,
            match="Provided Synapse id must be your input folder Synapse id "
                  "or a Synapse Id of a folder inside your input directory"):
        validate._check_parentid_permission_container(syn, parentid)


def test_valid__check_parentid_permission_container():
    '''
    Test that parentid specified is a container and have permissions to access
    '''
    parentid = "syn123"
    folder_ent = synapseclient.Folder("foo", parentId=parentid)
    with mock.patch.object(syn, "get", return_value=folder_ent):
        validate._check_parentid_permission_container(syn, parentid)


def test_valid__check_center_input():
    center = "FOO"
    center_list = ["FOO", "WOW"]
    validate._check_center_input(center, center_list)


def test_invalid__check_center_input():
    center = "BARFOO"
    center_list = ["FOO", "WOW"]
    with pytest.raises(
            ValueError,
            match="Must specify one of these centers: {}".format(
                  ", ".join(center_list))):
        validate._check_center_input(center, center_list)


ONCOTREE_ENT = 'syn222'


class argparser:
    oncotreelink = "link"
    parentid = None
    filetype = None
    testing = False
    center = "try"
    filepath = "path.csv"
    nosymbol_check = False

    def asDataFrame(self):
        database_dict = {"Database": ["centerMapping", 'oncotreeLink'],
                         "Id": ["syn123", ONCOTREE_ENT],
                         "center": ["try", 'foo']}
        databasetosynid_mappingdf = pd.DataFrame(database_dict)
        return(databasetosynid_mappingdf)


def test_notnone_get_oncotreelink():
    '''
    Test link passed in by user is used
    '''
    arg = argparser()
    url = "https://www.synapse.org"
    link = validate._get_oncotreelink(syn, arg.asDataFrame(), oncotreelink=url)
    assert link == url


def test_none__getoncotreelink():
    '''
    Test oncotree link is gotten
    '''
    arg = argparser()
    url = "https://www.synapse.org"
    link = synapseclient.File("foo", parentId="foo", externalURL=url)
    with mock.patch.object(syn, "get", return_value=link) as patch_synget:
        oncolink = validate._get_oncotreelink(syn, arg.asDataFrame())
        patch_synget.assert_called_once_with(ONCOTREE_ENT)
        assert oncolink == url


def test_valid__upload_to_synapse():
    '''
    Test upload of file to synapse under right conditions
    '''
    ent = synapseclient.File(id="syn123", parentId="syn222")
    with mock.patch.object(syn, "store", return_value=ent) as patch_synstore:
        validate._upload_to_synapse(syn, ['foo'], True, parentid="syn123")
        patch_synstore.assert_called_once_with(
            synapseclient.File('foo', parent="syn123"))


def test_perform_validate():
    '''
    Make sure all functions are called
    '''
    arg = argparser()
    check_perm_call = "genie.validate._check_parentid_permission_container"
    check_get_db_call = "genie.process_functions.get_synid_database_mappingdf"
    check_center_call = "genie.validate._check_center_input"
    validate_file_call = "genie.validate.Validator.validate_single_file"
    get_oncotree_call = "genie.validate._get_oncotreelink"
    upload_to_syn_call = "genie.validate._upload_to_synapse"
    valid = True
    with mock.patch(check_perm_call) as patch_check_parentid,\
        mock.patch(
            check_get_db_call,
            return_value=arg.asDataFrame()) as patch_getdb,\
        mock.patch.object(
            syn,
            "tableQuery",
            return_value=arg) as patch_syn_tablequery,\
        mock.patch(check_center_call) as patch_check_center,\
        mock.patch(get_oncotree_call) as patch_get_onco,\
        mock.patch(
            validate_file_call,
            return_value=(valid, 'foo', 'foo')) as patch_validate,\
        mock.patch(
            upload_to_syn_call) as patch_syn_upload:
        validate._perform_validate(syn, arg)
        patch_check_parentid.assert_called_once_with(syn, arg.parentid)
        patch_getdb.assert_called_once_with(syn, test=arg.testing)
        patch_syn_tablequery.assert_called_once_with('select * from syn123')
        patch_check_center.assert_called_once_with(arg.center, ["try", "foo"])
        patch_get_onco.assert_called_once()
        patch_validate.assert_called_once_with(arg.filetype,
                                               arg.oncotreelink, arg.testing,
                                               arg.nosymbol_check)
        patch_syn_upload.assert_called_once_with(
            syn, arg.filepath, valid, parentid=arg.parentid)
