"""Tests validate.py"""
from unittest import mock
from unittest.mock import Mock, patch

import pandas as pd
import pytest
import synapseclient
from synapseclient.core.exceptions import SynapseHTTPError

from genie import validate, process_functions, example_filetype_format

CENTER = "SAGE"
syn = mock.create_autospec(synapseclient.Synapse)

CNA_ENT = synapseclient.File(
    name="data_CNA_SAGE.txt", path="data_CNA_SAGE.txt", parentId="syn12345"
)
CLIN_ENT = synapseclient.File(
    name="data_clinical_supp_SAGE.txt",
    path="data_clinical_supp_SAGE.txt",
    parentId="syn12345",
)
SAMPLE_ENT = synapseclient.File(
    name="data_clinical_supp_sample_SAGE.txt",
    path="data_clinical_supp_sample_SAGE.txt",
    parentId="syn12345",
)
PATIENT_ENT = synapseclient.File(
    name="data_clinical_supp_patient_SAGE.txt",
    path="data_clinical_supp_patient_SAGE.txt",
    parentId="syn12345",
)
WRONG_NAME_ENT = synapseclient.File(
    name="wrong.txt", path="data_clinical_supp_SAGE.txt", parentId="syn12345"
)


class FileFormat(example_filetype_format.FileTypeFormat):
    _fileType = "clinical"


def test_perfect_determine_filetype():
    """
    Tests determining of file type through filenames
    Parameters are passed in from filename_fileformat_map
    """
    filetype = "clincial"
    ent_list = [SAMPLE_ENT]
    with patch.object(FileFormat, "validateFilename", return_value=filetype):
        validator = validate.GenieValidationHelper(
            syn, None, CENTER, ent_list, format_registry={filetype: FileFormat}
        )
        assert validator.determine_filetype() == filetype


def test_wrongfilename_noerror_determine_filetype():
    """
    Tests None is passed back when wrong filename is passed
    when raise_error flag is False
    """
    ent_list = [WRONG_NAME_ENT]
    with patch.object(FileFormat, "validateFilename", side_effect=AssertionError):
        validator = validate.GenieValidationHelper(
            syn,
            project_id=None,
            center=CENTER,
            entitylist=ent_list,
            format_registry={"wrong": FileFormat},
        )
        assert validator.file_type is None


def test_valid_collect_errors_and_warnings():
    """
    Tests if no error and warning strings are passed that
    returned valid and message is correct
    """
    results = example_filetype_format.ValidationResults(errors="", warnings="")
    message = results.collect_errors_and_warnings()
    assert message == "YOUR FILE IS VALIDATED!\n"


def test_invalid_collect_errors_and_warnings():
    """
    Tests if error and warnings strings are passed that
    returned valid and message is correct
    """
    results = example_filetype_format.ValidationResults(
        errors="error\nnow", warnings="warning\nnow"
    )
    message = results.collect_errors_and_warnings()
    assert message == (
        "----------------ERRORS----------------\n"
        "error\nnow"
        "-------------WARNINGS-------------\n"
        "warning\nnow"
    )


def test_warning_collect_errors_and_warnings():
    """
    Tests if no error but warnings strings are passed that
    returned valid and message is correct
    """
    results = example_filetype_format.ValidationResults(
        errors="", warnings="warning\nnow"
    )
    message = results.collect_errors_and_warnings()
    assert message == (
        "YOUR FILE IS VALIDATED!\n"
        "-------------WARNINGS-------------\n"
        "warning\nnow"
    )


def test_valid_validate_single_file(genie_config):
    """
    Tests that all the functions are run in validate single
    file workflow and all the right things are returned
    """
    entitylist = [CLIN_ENT]
    results = example_filetype_format.ValidationResults(errors="", warnings="")
    expected_message = "valid message here!"
    expected_filetype = "clinical"
    project_ent = Mock(id="syn1234")
    with patch.object(syn, "get", return_value=project_ent), patch.object(
        validate.GenieValidationHelper,
        "determine_filetype",
        return_value=expected_filetype,
    ) as mock_determine_ftype, patch.object(
        FileFormat, "validate", return_value=results
    ) as mock_genie_class, patch.object(
        results, "collect_errors_and_warnings", return_value=expected_message
    ) as mock_determine:
        validator = validate.GenieValidationHelper(
            syn,
            project_id="syn1234",
            center=CENTER,
            entitylist=entitylist,
            format_registry={"clinical": FileFormat},
            genie_config=genie_config,
        )
        valid_cls, message = validator.validate_single_file(nosymbol_check=False)

        assert valid_cls == results
        assert message == expected_message
        assert validator.file_type == expected_filetype

        mock_determine_ftype.assert_called_once_with()

        mock_genie_class.assert_called_once_with(
            filePathList=[CLIN_ENT.path], nosymbol_check=False, project_id="syn1234"
        )

        mock_determine.assert_called_once_with()


def test_filetype_validate_single_file():
    """
    Tests that if filetype is passed in that an error is thrown
    if it is an incorrect filetype
    """
    entitylist = [WRONG_NAME_ENT]
    expected_error = (
        "----------------ERRORS----------------\n"
        "Your filename is incorrect! Please change your "
        "filename before you run the validator or specify "
        "--filetype if you are running the validator locally"
    )

    with patch.object(FileFormat, "validateFilename", side_effect=AssertionError):
        validator = validate.GenieValidationHelper(
            syn, None, CENTER, entitylist, format_registry={"wrong": FileFormat}
        )

        valid_cls, message = validator.validate_single_file()
        assert message == expected_error
        assert not valid_cls.is_valid()


def test_wrongfiletype_validate_single_file():
    """
    Tests that if there is no filetype for the filename passed
    in, an error is thrown
    """
    entitylist = [WRONG_NAME_ENT]
    expected_error = (
        "----------------ERRORS----------------\n"
        "Your filename is incorrect! Please change your "
        "filename before you run the validator or specify "
        "--filetype if you are running the validator locally"
    )

    with patch.object(
        validate.GenieValidationHelper, "determine_filetype", return_value=None
    ) as mock_determine_filetype:
        validator = validate.GenieValidationHelper(
            syn=syn,
            project_id=None,
            center=CENTER,
            entitylist=entitylist,
            format_registry={"wrong": Mock()},
        )
        valid, message = validator.validate_single_file()

        assert message == expected_error
        assert not valid.is_valid()
        mock_determine_filetype.assert_called_once_with()


def test_nopermission__check_parentid_permission_container():
    """Throws error if no permissions to access"""
    parentid = "syn123"
    with patch.object(syn, "get", side_effect=SynapseHTTPError), pytest.raises(
        ValueError,
        match="Provided Synapse id must be your input folder "
        "Synapse id or a Synapse Id of a folder inside "
        "your input directory",
    ):
        validate._check_parentid_permission_container(syn, parentid)


def test_notcontainer__check_parentid_permission_container():
    """Throws error if input if synid of file"""
    parentid = "syn123"
    file_ent = synapseclient.File("foo", parentId=parentid)
    with patch.object(syn, "get", return_value=file_ent), pytest.raises(
        ValueError,
        match="Provided Synapse id must be your input folder "
        "Synapse id or a Synapse Id of a folder inside "
        "your input directory",
    ):
        validate._check_parentid_permission_container(syn, parentid)


def test_valid__check_parentid_permission_container():
    """
    Test that parentid specified is a container and have permissions to access
    """
    parentid = "syn123"
    folder_ent = synapseclient.Folder("foo", parentId=parentid)
    with patch.object(syn, "get", return_value=folder_ent):
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
        match="Must specify one of these " "centers: {}".format(", ".join(center_list)),
    ):
        validate._check_center_input(center, center_list)


ONCOTREE_ENT = "syn222"


class argparser:
    oncotree_link = "link"
    parentid = None
    filetype = None
    project_id = None
    center = "try"
    filepath = "path.csv"
    nosymbol_check = False
    format_registry_packages = ["genie"]
    project_id = "syn1234"

    def asDataFrame(self):
        database_dict = {
            "Database": ["centerMapping", "oncotreeLink"],
            "Id": ["syn123", ONCOTREE_ENT],
            "center": ["try", "foo"],
        }
        databasetosynid_mappingdf = pd.DataFrame(database_dict)
        return databasetosynid_mappingdf


def test_valid__upload_to_synapse():
    """
    Test upload of file to synapse under right conditions
    """
    ent = synapseclient.File(id="syn123", parentId="syn222")
    with patch.object(syn, "store", return_value=ent) as patch_synstore:
        validate._upload_to_synapse(syn, ["foo"], True, parentid="syn123")
        patch_synstore.assert_called_once_with(
            synapseclient.File("foo", parent="syn123")
        )


def test_perform_validate(genie_config):
    """Make sure all functions are called"""
    arg = argparser()
    valid = True
    with patch.object(
        validate, "_check_parentid_permission_container", return_value=None
    ) as patch_check_parentid, patch.object(
        process_functions, "get_genie_config", return_value=genie_config
    ) as patch_get_config, patch.object(
        validate, "_check_center_input"
    ) as patch_check_center, patch.object(
        process_functions, "_get_oncotreelink"
    ) as patch_get_onco, patch.object(
        validate.GenieValidationHelper,
        "validate_single_file",
        return_value=(valid, "foo"),
    ) as patch_validate, patch.object(
        validate, "_upload_to_synapse"
    ) as patch_syn_upload:
        validate._perform_validate(syn, arg)
        patch_check_parentid.assert_called_once_with(syn=syn, parentid=arg.parentid)
        patch_get_config.assert_called_once_with(syn=syn, project_id=arg.project_id)
        patch_check_center.assert_called_once_with(arg.center, ["SAGE", "TEST"])
        patch_get_onco.assert_called_once()
        patch_validate.assert_called_once_with(
            nosymbol_check=arg.nosymbol_check, project_id=arg.project_id
        )
        patch_syn_upload.assert_called_once_with(
            syn, arg.filepath, valid, parentid=arg.parentid
        )
