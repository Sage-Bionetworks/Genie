from datetime import datetime
from unittest import mock
from unittest.mock import Mock, patch

import pandas as pd
import pytest
import synapseclient
import synapseutils

from genie import (
    input_to_database,
    example_filetype_format,
    process_functions,
    process_mutation,
)
from genie_registry.clinical import Clinical
from genie.validate import GenieValidationHelper


syn = mock.create_autospec(synapseclient.Synapse)
sample_clinical_synid = "syn2222"

sample_clinical_entity = synapseclient.File(
    path="data_clinical_supp_sample_SAGE.txt",
    id=sample_clinical_synid,
    parentId="syn45678",
    name="data_clinical_supp_sample_SAGE.txt",
    modifiedOn="2019-03-24T12:00:00.Z",
    md5="44444",
)

patient_clinical_synid = "syn11111"
patient_clinical_entity = synapseclient.File(
    path="data_clinical_supp_patient_SAGE.txt",
    id=patient_clinical_synid,
    parentId="syn45678",
    name="data_clinical_supp_patient_SAGE.txt",
)

vcf1synid = "syn6666"
vcf1_entity = synapseclient.File(
    path="GENIE-SAGE-1-1.vcf",
    id=vcf1synid,
    parentId="syn45678",
    name="GENIE-SAGE-1-1.vcf",
)
vcf2synid = "syn8888"
vcf2_entity = synapseclient.File(
    path="GENIE-SAGE-2-1.vcf",
    id=vcf2synid,
    parentId="syn45678",
    name="GENIE-SAGE-2-1.vcf",
)
first = (
    [("inputs", "syn12345")],
    [("vcfs", "syn33333")],
    [
        ("data_clinical_supp_sample_SAGE.txt", sample_clinical_synid),
        ("data_clinical_supp_patient_SAGE.txt", patient_clinical_synid),
    ],
)
second = (
    [("vcfs", "syn33333")],
    [],
    [("GENIE-SAGE-000-1111.vcf", vcf1synid), ("GENIE-SAGE-111-2222.vcf", vcf2synid)],
)
center = "SAGE"
oncotree_link = (
    "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"
)
center_input_synid = "syn9999"
center_staging_synid = "syn9999"
center_mapping = {
    "inputSynId": [center_input_synid],
    "stagingSynId": [center_input_synid],
    "center": [center],
}
center_mapping_df = pd.DataFrame(center_mapping)
validation_statusdf = pd.DataFrame(
    {
        "id": ["syn1234", "syn2345"],
        "status": ["VALID", "INVALID"],
        "md5": ["3333", "44444"],
        "name": ["first.txt", "second.txt"],
        "fileType": ["filetype1", "filetype2"],
    }
)
error_trackerdf = pd.DataFrame(
    {"id": ["syn2345"], "errors": ["Invalid file format"], "fileType": ["filetype1"]}
)
emptydf = pd.DataFrame(columns=["id"], dtype=str)


class mock_csv_query_result(object):
    def __init__(self, df):
        self.df = df

    def asDataFrame(self):
        return self.df


# def test_samename_rename_file():
#     '''Test that the file path is not renamed.
#     '''
#     filename = synapseclient.utils.make_bogus_data_file()
#     entity = synapseclient.File(path=filename,
#                                 id='syn012345',
#                                 parentId='syn45678',
#                                 name=os.path.basename(filename))
#     expectedpath = filename
#     new_entity = input_to_database.rename_file(entity)
#     assert new_entity.annotations.expectedPath == expectedpath
#     os.remove(filename)


# def test_diffname_rename_file():
#     '''Test that the file path is renamed.
#     '''
#     filename = synapseclient.utils.make_bogus_data_file()
#     entity = synapseclient.File(path=filename,
#                                 id='syn012345',
#                                 parentId='syn45678',
#                                 name='testname')

#     expectedpath = os.path.join(os.path.dirname(filename), "testname")
#     new_entity = input_to_database.rename_file(entity)
#     assert new_entity.annotations.expectedPath == expectedpath
#     os.remove(filename)


def walk_return():
    """
    Generator returned by synapseutils.walk
    """
    yield first
    yield second


def walk_return_empty():
    """
    Generator returned by synapseutils.walk
    """
    yield ([], [], [])


def test_main_get_center_input_files():
    """
    Test to make sure center input files are gotten
    excluding the vcf files since process main is specified
    """
    syn_get_effects = [sample_clinical_entity, patient_clinical_entity]
    expected_center_file_list = [syn_get_effects]

    calls = [
        mock.call(sample_clinical_synid, downloadFile=True),
        mock.call(patient_clinical_synid, downloadFile=True),
    ]

    with patch.object(
        synapseutils, "walk", return_value=walk_return()
    ) as patch_synapseutils_walk, patch.object(
        syn, "get", side_effect=syn_get_effects
    ) as patch_syn_get:
        center_file_list = input_to_database.get_center_input_files(
            syn, "syn12345", center
        )

        assert len(center_file_list) == len(expected_center_file_list)
        assert len(center_file_list[0]) == 2
        assert center_file_list == expected_center_file_list
        patch_synapseutils_walk.assert_called_once_with(syn, "syn12345")
        patch_syn_get.assert_has_calls(calls)


def test_mutation_get_center_input_files():
    """
    Test to make sure center input files are gotten
    including the vcf files since process vcf is specified
    """
    syn_get_effects = [
        sample_clinical_entity,
        patient_clinical_entity,
        vcf1_entity,
        vcf2_entity,
    ]
    expected_center_file_list = [
        [vcf1_entity],
        [vcf2_entity],
        [sample_clinical_entity, patient_clinical_entity],
    ]
    calls = [
        mock.call(sample_clinical_synid, downloadFile=True),
        mock.call(patient_clinical_synid, downloadFile=True),
        mock.call(vcf1synid, downloadFile=True),
        mock.call(vcf2synid, downloadFile=True),
    ]

    with patch.object(
        synapseutils, "walk", return_value=walk_return()
    ) as patch_synapseutils_walk, patch.object(
        syn, "get", side_effect=syn_get_effects
    ) as patch_syn_get:
        center_file_list = input_to_database.get_center_input_files(
            syn, "syn12345", center, process="mutation"
        )
        assert len(center_file_list) == len(expected_center_file_list)
        assert len(center_file_list[2]) == 2
        assert center_file_list == expected_center_file_list
        patch_synapseutils_walk.assert_called_once_with(syn, "syn12345")
        patch_syn_get.assert_has_calls(calls)


def test_empty_get_center_input_files():
    """
    Test that center input files is empty if directory
    pass in is empty
    """
    with patch.object(
        synapseutils, "walk", return_value=walk_return_empty()
    ) as patch_synapseutils_walk:
        center_file_list = input_to_database.get_center_input_files(
            syn, "syn12345", center, process="vcf"
        )
        assert center_file_list == []
        patch_synapseutils_walk.assert_called_once_with(syn, "syn12345")


# @pytest.fixture(params=[
#     # tuple with (input, expectedOutput)
#     (["data_CNA_SAGE.txt"], "cna"),
#     (["data_clinical_supp_SAGE.txt"], "clinical"),
#     (["data_clinical_supp_sample_SAGE.txt",
#       "data_clinical_supp_patient_SAGE.txt"], "clinical")])
# def filename_fileformat_map(request):
#     return request.param


# def test_perfect_get_filetype(filename_fileformat_map):
#     (filepath_list, fileformat) = filename_fileformat_map
#     assert input_to_database.get_filetype(
#         syn, filepath_list, center) == fileformat


# def test_wrongfilename_get_filetype():
#     assert input_to_database.get_filetype(syn, ['wrong.txt'], center) is None


def test_unvalidatedinput_check_existing_file_status():
    """
    Test the values returned by input that hasn't be validated
    """
    entity = synapseclient.Entity(id="syn1234")
    entities = [entity]

    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(emptydf), mock_csv_query_result(emptydf), entities
    )
    assert file_status["to_validate"]
    assert file_status["status_list"] == []
    assert file_status["error_list"] == []


def test_valid_check_existing_file_status():
    """
    Test the values returned by input that is already valid
    """
    entity = synapseclient.Entity(name="first.txt", id="syn1234", md5="3333")
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf),
        mock_csv_query_result(error_trackerdf),
        entities,
    )
    assert not file_status["to_validate"]
    assert file_status["status_list"] == ["VALID"]
    assert file_status["error_list"] == []


def test_invalid_check_existing_file_status():
    """
    Test the values returned by input that is invalid
    """
    entity = synapseclient.Entity(name="second.txt", id="syn2345", md5="44444")
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf),
        mock_csv_query_result(error_trackerdf),
        entities,
    )
    assert not file_status["to_validate"]
    assert file_status["status_list"] == ["INVALID"]
    assert file_status["error_list"] == ["Invalid file format"]


def test_nostorederrors_check_existing_file_status():
    """
    If there is no error uploaded, must re-validate file
    """
    entity = synapseclient.Entity(name="second.txt", id="syn2345", md5="44444")
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf),
        mock_csv_query_result(emptydf),
        entities,
    )
    assert file_status["to_validate"]
    assert file_status["status_list"] == ["INVALID"]
    assert file_status["error_list"] == []


def test_diffmd5validate_check_existing_file_status():
    """
    If md5 is different from stored md5, must re-validate file
    """
    entity = synapseclient.Entity(name="first.txt", id="syn1234", md5="44444")
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf),
        mock_csv_query_result(emptydf),
        entities,
    )
    assert file_status["to_validate"]
    assert file_status["status_list"] == ["VALID"]
    assert file_status["error_list"] == []


def test_diffnametovalidate_check_existing_file_status():
    """
    If name is different from stored name, must re-validate file
    """
    entity = synapseclient.Entity(name="second.txt", id="syn1234", md5="3333")
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf),
        mock_csv_query_result(emptydf),
        entities,
    )
    assert file_status["to_validate"]
    assert file_status["status_list"] == ["VALID"]
    assert file_status["error_list"] == []


def test_twoinvalid_check_existing_file_status():
    """
    Test to make sure both invalid statues get passed back
    """
    validation_statusdf = pd.DataFrame(
        {
            "id": ["syn1234", "syn2345"],
            "status": ["INVALID", "INVALID"],
            "md5": ["3333", "44444"],
            "name": ["first.txt", "second.txt"],
        }
    )
    error_trackerdf = pd.DataFrame(
        {
            "id": ["syn1234", "syn2345"],
            "errors": ["Invalid file format", "Invalid formatting issues"],
        }
    )
    first_entity = synapseclient.Entity(name="first.txt", id="syn1234", md5="3333")
    second_entity = synapseclient.Entity(name="second.txt", id="syn2345", md5="44444")
    entities = [first_entity, second_entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf),
        mock_csv_query_result(error_trackerdf),
        entities,
    )
    assert not file_status["to_validate"]
    assert file_status["status_list"] == ["INVALID", "INVALID"]
    assert file_status["error_list"] == [
        "Invalid file format",
        "Invalid formatting issues",
    ]


def test_error_check_existing_file_status():
    """
    With exception of clinical files, there should never be
    more than 2 files being passed in for validation.
    """
    with pytest.raises(
        ValueError, match="There should never be more than 2 files being validated."
    ):
        entities = ["foo", "doo", "boo"]
        input_to_database.check_existing_file_status(emptydf, emptydf, entities)


def test_valid_validatefile(genie_config):
    """
    Tests the behavior of a file that gets validated that becomes
    valid
    """
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame()
    entity = synapseclient.Entity(
        name="data_clinical_supp_SAGE.txt",
        id="syn1234",
        md5="44444",
        path="/path/to/data_clinical_supp_SAGE.txt",
    )
    entity["modifiedOn"] = "2019-03-24T12:00:00.Z"
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = "333"
    entity.createdBy = "444"
    entities = [entity]
    valid_results = example_filetype_format.ValidationResults(errors="", warnings="")
    message = "Is valid"
    filetype = "clinical"
    # Only a list is returned as opposed a list of lists when there are
    # invalid errors
    status_error_list_results = ([{"entity": entity, "status": "VALIDATED"}], [])
    expected_results = (
        [
            {
                "entity": entity,
                "status": "VALIDATED",
                "fileType": filetype,
                "center": center,
            }
        ],
        [],
        [],
    )
    with patch.object(
        GenieValidationHelper, "determine_filetype", return_value=filetype
    ) as patch_determine_filetype, patch.object(
        input_to_database,
        "check_existing_file_status",
        return_value={"status_list": [], "error_list": [], "to_validate": True},
    ) as patch_check, patch.object(
        GenieValidationHelper,
        "validate_single_file",
        return_value=(valid_results, message),
    ) as patch_validate, patch.object(
        input_to_database,
        "_get_status_and_error_list",
        return_value=status_error_list_results,
    ) as patch_staterror_list, patch.object(
        input_to_database, "_send_validation_error_email"
    ) as patch_send_email:

        validate_results = input_to_database.validatefile(
            syn,
            None,
            entities,
            validation_statusdf,
            error_trackerdf,
            center,
            genie_config=genie_config,
        )

        assert expected_results == validate_results
        patch_validate.assert_called_once_with(
            oncotree_link=oncotree_link, nosymbol_check=False
        )
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities
        )
        patch_determine_filetype.assert_called_once()
        patch_staterror_list.assert_called_once_with(
            valid_results.is_valid(), message, entities
        )
        patch_send_email.assert_not_called()


def test_invalid_validatefile(genie_config):
    """
    Tests the behavior of a file that gets validated that becomes
    invalid
    """
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame(columns=["id"], dtype=str)
    entity = synapseclient.Entity(
        name="data_clinical_supp_SAGE.txt",
        id="syn1234",
        md5="44444",
        path="/path/to/data_clinical_supp_SAGE.txt",
    )
    entity["modifiedOn"] = "2019-03-24T12:00:00.Z"
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = "333"
    entity.createdBy = "444"
    entities = [entity]
    valid_results = example_filetype_format.ValidationResults(
        errors="Is invalid", warnings=""
    )
    message = "Is invalid"
    filetype = "clinical"
    status_error_list_results = (
        [{"entity": entity, "status": "INVALID"}],
        [{"entity": entity, "errors": message}],
    )
    expected_results = (
        [
            {
                "entity": entity,
                "status": "INVALID",
                "fileType": filetype,
                "center": center,
            }
        ],
        [{"entity": entity, "errors": message, "fileType": filetype, "center": center}],
        [(["data_clinical_supp_SAGE.txt"], "Is invalid", ["333", "444"])],
    )

    with patch.object(
        GenieValidationHelper, "determine_filetype", return_value=filetype
    ) as patch_determine_filetype, patch.object(
        input_to_database,
        "check_existing_file_status",
        return_value={"status_list": [], "error_list": [], "to_validate": True},
    ) as patch_check, patch.object(
        GenieValidationHelper,
        "validate_single_file",
        return_value=(valid_results, message),
    ) as patch_validate, patch.object(
        input_to_database,
        "_get_status_and_error_list",
        return_value=status_error_list_results,
    ) as patch_staterror_list:

        validate_results = input_to_database.validatefile(
            syn,
            None,
            entities,
            validation_statusdf,
            error_trackerdf,
            center,
            genie_config=genie_config,
        )

        assert expected_results == validate_results
        patch_validate.assert_called_once_with(
            oncotree_link=oncotree_link, nosymbol_check=False
        )
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities
        )
        patch_determine_filetype.assert_called_once()
        patch_staterror_list.assert_called_once_with(
            valid_results.is_valid(), message, entities
        )


def test_already_validated_validatefile():
    """
    Test already validated files
    """
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame()
    entity = synapseclient.Entity(
        name="data_clinical_supp_SAGE.txt",
        id="syn1234",
        md5="44444",
        path="/path/to/data_clinical_supp_SAGE.txt",
    )
    entity["modifiedOn"] = "2019-03-24T12:00:00.Z"
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = "333"
    entity.createdBy = "444"
    entities = [entity]
    threads = 0
    valid = False
    errors = "Invalid file"
    filetype = "markdown"
    status = "INVALID"
    check_file_status_dict = {
        "status_list": [status],
        "error_list": [errors],
        "to_validate": False,
    }

    status_error_list_results = (
        [{"entity": entity, "status": status}],
        [{"entity": entity, "errors": errors}],
    )
    expected_results = (
        [{"entity": entity, "status": status, "fileType": filetype, "center": center}],
        [{"entity": entity, "errors": errors, "fileType": filetype, "center": center}],
        [],
    )
    with patch.object(
        GenieValidationHelper, "determine_filetype", return_value=filetype
    ) as patch_determine_filetype, patch.object(
        input_to_database,
        "check_existing_file_status",
        return_value=check_file_status_dict,
    ) as patch_check, patch.object(
        GenieValidationHelper, "validate_single_file", return_value=(valid, errors)
    ) as patch_validate, patch.object(
        input_to_database,
        "_get_status_and_error_list",
        return_value=status_error_list_results,
    ) as patch_staterror_list, patch.object(
        input_to_database, "_send_validation_error_email"
    ) as patch_send_email:

        validate_results = input_to_database.validatefile(
            syn,
            None,
            entities,
            validation_statusdf,
            error_trackerdf,
            center,
            threads,
            oncotree_link,
        )

        assert expected_results == validate_results
        patch_validate.assert_not_called()
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities
        )
        patch_determine_filetype.assert_called_once()
        patch_staterror_list.assert_not_called()
        patch_send_email.assert_not_called()


def test_valid__get_status_and_error_list():
    """
    Tests the correct status and error lists received
    when file is valid.
    """
    # modified_on = 1561143558000
    modified_on_string = "2019-06-21T18:59:18.456Z"

    entity = synapseclient.Entity(
        id="syn1234",
        md5="44444",
        path="/path/to/foobar.txt",
        name="data_clinical_supp_SAGE.txt",
    )
    entity.properties.modifiedOn = modified_on_string

    entities = [entity]

    valid = True
    message = "valid"
    # filetype = 'clinical'

    (
        input_status_list,
        invalid_errors_list,
    ) = input_to_database._get_status_and_error_list(valid, message, entities)
    assert input_status_list == [{"entity": entity, "status": "VALIDATED"}]
    assert not invalid_errors_list


def test_invalid__get_status_and_error_list():
    """
    Tests the correct status and error lists received
    when file is invalid.
    """
    # modified_on = 1561143558000
    modified_on_string = "2019-06-21T18:59:18.456Z"
    entity = synapseclient.Entity(
        id="syn1234",
        md5="44444",
        path="/path/to/foobar.txt",
        name="data_clinical_supp_SAGE.txt",
    )
    entity.properties.modifiedOn = modified_on_string

    entities = [entity]
    # filetype = "clinical"
    # This valid variable control the validation status
    valid = False
    errors = "invalid file content"

    (
        input_status_list,
        invalid_errors_list,
    ) = input_to_database._get_status_and_error_list(valid, errors, entities)
    assert input_status_list == [{"entity": entity, "status": "INVALID"}]
    assert invalid_errors_list == [{"entity": entity, "errors": errors}]


def test__send_validation_error_email():
    message = "invalid error message here"
    filenames = ["data_clinical_supp_SAGE.txt"]
    file_user = "333"
    message_objs = [dict(filenames=filenames, messages=message)]
    with patch.object(
        syn, "getUserProfile", return_value={"userName": "trial"}
    ) as patch_syn_getuserprofile, patch.object(
        syn, "sendMessage"
    ) as patch_syn_sendmessage, patch(
        "genie.input_to_database.datetime"
    ) as mock_datetime:
        mock_datetime.datetime.today.return_value = datetime(2019, 11, 1, 00, 00, 00, 0)
        mock_datetime.side_effect = lambda *args, **kw: date(*args, **kw)

        input_to_database._send_validation_error_email(syn, file_user, message_objs)
        error_message = (
            "Dear trial,\n\n"
            "You have invalid files! "
            "Here are the reasons why:\n\n"
            "Filenames: data_clinical_supp_SAGE.txt, "
            f"Errors:\n {message}\n\n"
        )
        patch_syn_sendmessage.assert_called_once_with(
            userIds=[file_user],
            messageBody=error_message,
            messageSubject="GENIE Validation Error - 2019-11-01 00:00:00",
        )
        patch_syn_getuserprofile.assert_called_once_with(file_user)


class emptytable_mock:
    """
    Validation status tablequery dataframe mocking and
    error tracking tablequery dataframe mocking
    This is used because assert_called_once_with has a hard
    time with comparing pandas dataframes
    """

    tableId = "syn555"

    def asDataFrame(self):
        return []


class TestValidation:
    def setup(self):
        valid = [
            [
                sample_clinical_entity.id,
                sample_clinical_entity.path,
                sample_clinical_entity.md5,
                "VALIDATED",
                sample_clinical_entity.name,
                1553428800000,
                "clinical",
                center,
                sample_clinical_entity,
            ]
        ]
        self.empty_validation = pd.DataFrame(
            columns=[
                "id",
                "path",
                "md5",
                "status",
                "name",
                "modifiedOn",
                "fileType",
                "center",
                "entity",
            ]
        )
        self.validation_statusdf = pd.DataFrame(
            valid,
            columns=[
                "id",
                "path",
                "md5",
                "status",
                "name",
                "modifiedOn",
                "fileType",
                "center",
                "entity",
            ],
        )
        error = [
            [
                sample_clinical_entity.id,
                "errors",
                sample_clinical_entity.name,
                "clinical",
                center,
                sample_clinical_entity,
            ]
        ]
        self.errors_df = pd.DataFrame(
            error, columns=["id", "errors", "name", "fileType", "center", "entity"]
        )
        self.empty_errors = pd.DataFrame(
            columns=["id", "errors", "name", "fileType", "center", "entity"]
        )

        self.with_dupsdf = pd.DataFrame(
            {
                "id": [
                    sample_clinical_entity.id,
                    "syn2345",
                    "syn5555",
                    "syn1224",
                    "syn34444",
                ],
                "name": [
                    "first.cbs",
                    "second.seg",
                    "data_clinical_supp_1",
                    "data_clinical_supp_2",
                    "data_clinical_supp_3",
                ],
                "center": ["SAGE"] * 5,
                "fileType": ["type"] * 5,
                "entity": ["entity"] * 5,
            }
        )
        self.duplicateddf = self.with_dupsdf.copy()
        self.duplicateddf["errors"] = input_to_database.DUPLICATED_FILE_ERROR

        self.no_dupsdf = pd.DataFrame(
            {
                "id": [
                    sample_clinical_entity.id,
                    "syn2345",
                    "syn5555",
                    "syn1224",
                    "syn34444",
                ],
                "name": [
                    "cbs.txt",
                    "second.seg",
                    "no_clinical.txt",
                    "data_clinical_supp_2",
                    "data_clinical_supp_3",
                ],
                "center": ["SAGE"] * 5,
                "fileType": ["type"] * 5,
                "entity": ["entity"] * 5,
            }
        )
        self.empty_dup = pd.DataFrame(
            columns=["id", "name", "center", "fileType", "entity", "errors"]
        )

    def test_build_validation_status_table(self):
        input_valid_statuses = [
            {
                "entity": sample_clinical_entity,
                "status": "VALIDATED",
                "fileType": "clinical",
                "center": center,
            }
        ]
        validationdf = input_to_database.build_validation_status_table(
            input_valid_statuses
        )
        assert validationdf.equals(self.validation_statusdf)

    def test_build_validation_status_table__empty(self):
        input_valid_statuses = []
        validationdf = input_to_database.build_validation_status_table(
            input_valid_statuses
        )
        assert validationdf.equals(self.empty_validation)

    def test_build_error_tracking_table(self):
        invalid_errors = [
            {
                "entity": sample_clinical_entity,
                "errors": "errors",
                "fileType": "clinical",
                "center": center,
            }
        ]
        errordf = input_to_database.build_error_tracking_table(invalid_errors)
        assert errordf.equals(self.errors_df)

    def test_build_error_tracking_table__empty(self):
        invalid_errors = []
        errordf = input_to_database.build_error_tracking_table(invalid_errors)
        assert errordf.equals(self.empty_errors)

    def test_update_status_and_error_tables(self):
        """Test updating validation status and error table"""
        validation_status_table = emptytable_mock()
        error_tracker_table = emptytable_mock()

        with patch.object(process_functions, "updateDatabase") as mock_update:
            input_to_database.update_status_and_error_tables(
                syn,
                self.validation_statusdf,
                self.errors_df,
                validation_status_table,
                error_tracker_table,
            )
            assert mock_update.call_count == 2

    def test_dups_get_duplicated_files(self):
        """Test get all duplicates: cbs/seg/clinical"""
        dupsdf = input_to_database.get_duplicated_files(self.with_dupsdf)
        assert dupsdf.equals(self.duplicateddf)

    def test_nodups_get_duplicated_files(self):
        """Test no duplicated"""
        dupsdf = input_to_database.get_duplicated_files(self.no_dupsdf)
        assert dupsdf.equals(self.empty_dup)

    def test__update_tables_content(self):
        """Tests duplicates are added to the tables and errors/statues are
        updated
        """
        errorsdf = self.errors_df.copy()
        errorsdf["errors"] = input_to_database.DUPLICATED_FILE_ERROR

        validationdf = self.validation_statusdf
        validationdf["status"] = "INVALID"

        with patch.object(
            input_to_database, "get_duplicated_files", return_value=self.duplicateddf
        ):
            updated_tables = input_to_database._update_tables_content(
                self.validation_statusdf, self.errors_df
            )
        assert updated_tables["duplicated_filesdf"].equals(self.duplicateddf)
        assert updated_tables["error_trackingdf"].equals(errorsdf)
        assert updated_tables["validation_statusdf"].equals(validationdf)

    def test__update_tables_content__remove_old_duplicates(self):
        """Tests that old duplicates are removed"""
        errorsdf = self.errors_df.copy()
        errorsdf["errors"] = input_to_database.DUPLICATED_FILE_ERROR
        validationdf = self.validation_statusdf
        validationdf["status"] = "INVALID"

        with patch.object(
            input_to_database, "get_duplicated_files", return_value=self.empty_dup
        ):
            updated_tables = input_to_database._update_tables_content(
                validationdf, errorsdf
            )
        assert updated_tables["duplicated_filesdf"].empty
        assert updated_tables["error_trackingdf"].empty
        assert updated_tables["validation_statusdf"].empty

    def test__update_tables_content__remove_old_errors(self):
        """Tests that old errors are removed"""
        errorsdf = self.errors_df.copy()
        errorsdf["errors"] = input_to_database.DUPLICATED_FILE_ERROR
        validationdf = self.validation_statusdf
        validationdf["status"] = "VALIDATED"

        with patch.object(
            input_to_database, "get_duplicated_files", return_value=self.empty_dup
        ):
            updated_tables = input_to_database._update_tables_content(
                validationdf, errorsdf
            )
        assert updated_tables["duplicated_filesdf"].empty
        assert updated_tables["error_trackingdf"].empty
        assert updated_tables["validation_statusdf"].empty

    def test_validation(self, genie_config):
        """Test validation steps"""
        modified_on = 1561143558000
        process = "main"
        entity = synapseclient.Entity(
            id="syn1234",
            md5="44444",
            path="/path/to/foobar.txt",
            name="data_clinical_supp_SAGE.txt",
        )
        entities = [entity]
        filetype = "clinical"
        input_status_list = [
            [
                entity.id,
                entity.path,
                entity.md5,
                "VALIDATED",
                entity.name,
                modified_on,
                filetype,
                center,
            ]
        ]
        invalid_errors_list = []
        messages = []
        new_tables = {
            "validation_statusdf": self.validation_statusdf,
            "error_trackingdf": self.errors_df,
            "duplicated_filesdf": self.empty_dup,
        }
        validationstatus_mock = emptytable_mock()
        errortracking_mock = emptytable_mock()
        valiate_cls = Mock()
        with patch.object(
            syn, "tableQuery", side_effect=[validationstatus_mock, errortracking_mock]
        ) as patch_query, patch.object(
            input_to_database,
            "validatefile",
            return_value=(input_status_list, invalid_errors_list, messages),
        ) as patch_validatefile, patch.object(
            input_to_database,
            "build_validation_status_table",
            return_value=self.validation_statusdf,
        ), patch.object(
            input_to_database, "build_error_tracking_table", return_value=self.errors_df
        ), patch.object(
            input_to_database, "_update_tables_content", return_value=new_tables
        ), patch.object(
            input_to_database, "update_status_and_error_tables"
        ):

            valid_filedf = input_to_database.validation(
                syn,
                "syn123",
                center,
                process,
                entities,
                format_registry={"test": valiate_cls},
                genie_config=genie_config,
            )
            assert patch_query.call_count == 2
            patch_validatefile.assert_called_once_with(
                syn=syn,
                project_id="syn123",
                entities=entity,
                validation_status_table=validationstatus_mock,
                error_tracker_table=errortracking_mock,
                center="SAGE",
                format_registry={"test": valiate_cls},
                genie_config=genie_config,
            )

            assert valid_filedf.equals(
                self.validation_statusdf[["id", "path", "fileType", "name"]]
            )


@pytest.mark.parametrize(
    "process, genieclass, filetype",
    [
        ("main", Mock(), "clinical"),
        ("main", Mock(), "maf"),
    ],
)
def test_main_processfile(genie_config, process, genieclass, filetype):
    validfiles = {
        "id": ["syn1"],
        "path": ["/path/to/data_clinical_supp.txt"],
        "fileType": [filetype],
        "name": ["data_clinical_supp_SAGE.txt"],
    }
    validfilesdf = pd.DataFrame(validfiles)
    center = "SAGE"
    path_to_genie = "./"
    format_registry = {filetype: genieclass}

    input_to_database.processfiles(
        syn,
        validfilesdf,
        center,
        path_to_genie,
        processing=process,
        format_registry=format_registry,
        genie_config=genie_config,
    )
    genieclass.assert_called_once()


def test_mainnone_processfile(genie_config):
    """If file type is None, the processing function is not called"""
    validfiles = {
        "id": ["syn1"],
        "path": ["/path/to/data_clinical_supp_SAGE.txt"],
        "fileType": [None],
        "name": ["data_clinical_supp_SAGE.txt"],
    }
    validfilesdf = pd.DataFrame(validfiles)
    center = "SAGE"
    path_to_genie = "./"
    process_cls = Mock()

    with patch.object(Clinical, "process") as patch_clin:
        input_to_database.processfiles(
            syn,
            validfilesdf,
            center,
            path_to_genie,
            processing="main",
            format_registry={"main": process_cls},
            genie_config=genie_config,
        )
        patch_clin.assert_not_called()


def test_mutation_processfile(genie_config):
    """
    Make sure mutation is called correctly
    """
    validfiles = {
        "id": ["syn1"],
        "path": ["/path/to/data_clinical_supp_SAGE.txt"],
        "fileType": [None],
    }
    validfilesdf = pd.DataFrame(validfiles)
    center = "SAGE"
    path_to_genie = "./"
    process_cls = Mock()

    with patch.object(process_mutation, "process_mutation_workflow") as patch_process:
        input_to_database.processfiles(
            syn,
            validfilesdf,
            center,
            path_to_genie,
            processing="mutation",
            format_registry={"vcf": process_cls},
            genie_config=genie_config,
        )
        patch_process.assert_called_once_with(
            syn=syn,
            center=center,
            validfiles=validfilesdf,
            genie_config=genie_config,
            workdir=path_to_genie,
        )
