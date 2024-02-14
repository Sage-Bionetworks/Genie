"""Tests validate.py"""

from unittest.mock import Mock, patch

import pandas as pd
import pytest
import synapseclient
from synapseclient.core.exceptions import SynapseHTTPError

from genie import example_filetype_format, extract, load, validate, process_functions

CENTER = "SAGE"
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


def test_perfect_determine_filetype(syn):
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


def test_wrongfilename_noerror_determine_filetype(syn):
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


def test_valid_validate_single_file(syn, genie_config):
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


def test_filetype_validate_single_file(syn):
    """
    Tests that if filetype is passed in that an error is thrown
    if it is an incorrect filetype
    """
    entitylist = [WRONG_NAME_ENT]
    expected_error = (
        "----------------ERRORS----------------\n"
        "Your filename is incorrect! Please change your "
        "filename before you run the validator or specify "
        "--filetype if you are running the validator locally. "
        "If specifying filetype, options are: [wrong]\n"
    )

    with patch.object(FileFormat, "validateFilename", side_effect=AssertionError):
        validator = validate.GenieValidationHelper(
            syn, None, CENTER, entitylist, format_registry={"wrong": FileFormat}
        )

        valid_cls, message = validator.validate_single_file()
        assert message == expected_error
        assert not valid_cls.is_valid()


def test_wrongfiletype_validate_single_file(syn):
    """
    Tests that if there is no filetype for the filename passed
    in, an error is thrown
    """
    entitylist = [WRONG_NAME_ENT]
    expected_error = (
        "----------------ERRORS----------------\n"
        "Your filename is incorrect! Please change your "
        "filename before you run the validator or specify "
        "--filetype if you are running the validator locally. "
        "If specifying filetype, options are: [wrong]\n"
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


def test_accepted_chromosomes_value():
    """Testing global accepted_chromosomes before we use it"""
    assert validate.ACCEPTED_CHROMOSOMES == [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "X",
        "Y",
        "MT",
    ]


def test_nopermission__check_parentid_permission_container(syn):
    """Throws error if no permissions to access"""
    parentid = "syn123"
    with patch.object(syn, "get", side_effect=SynapseHTTPError), pytest.raises(
        ValueError,
        match="Provided Synapse id must be your input folder "
        "Synapse id or a Synapse Id of a folder inside "
        "your input directory",
    ):
        validate._check_parentid_permission_container(syn, parentid)


def test_notcontainer__check_parentid_permission_container(syn):
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


def test_valid__check_parentid_permission_container(syn):
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


def test_nonexistentcol__validate_chromosome():
    """Checks that no errors or warnings get thrown if
    chromosome vol doesn't exist"""
    input_df = pd.DataFrame(
        {
            "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID1-1", "ID3-1"],
            "SITE1_HUGO_SYMBOL": ["af", "af", "ff"],
            "SITE2_HUGO_SYMBOL": ["af", "af", "ff"],
        }
    )
    errors, warnings = validate._validate_chromosome(
        df=input_df, col="CHROM", fileformat="FOO", allow_chr=True
    )
    assert (
        errors == "" and warnings == ""
    ), "Error and warnings should be empty when chromosome col doesn't exist!"


def test_valid_nochar__validate_chromosome():
    """Checks that no errors or warnings get thrown if
    chromosome col has all the valid values and chr is not allowed"""
    input_df = pd.DataFrame(
        {
            "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID1-1", "ID3-1"],
            "CHROM": ["1", "19", "MT"],
        }
    )
    errors, warnings = validate._validate_chromosome(
        df=input_df, col="CHROM", fileformat="FOO", allow_chr=False
    )
    assert (
        errors == "" and warnings == ""
    ), "Error and warnings should be empty when chromosome col is valid!"


def test_valid_allowchr__validate_chromosome():
    """Checks that no errors get thrown if
    chromosome col has all the valid values and has chr in the values"""
    input_df = pd.DataFrame(
        {
            "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID1-1", "ID3-1"],
            "CHROM": ["1", "chr19", "chrMT"],
        }
    )
    errors, warnings = validate._validate_chromosome(
        df=input_df, col="CHROM", fileformat="FOO", allow_chr=True
    )
    assert errors == ""
    assert warnings == "FOO: Should not have the chr prefix in front of chromosomes.\n"


def test_invalid_allowchr__validate_chromosome():
    """Checks that errors and warnings get thrown if
    chromosome col has invalid values and chr is allowed in the values"""
    input_df = pd.DataFrame(
        {
            "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID1-1", "ID3-1"],
            "CHROM": ["chr1", "chr9", "chrZ"],
        }
    )
    errors, warnings = validate._validate_chromosome(
        df=input_df, col="CHROM", fileformat="FOO", allow_chr=True
    )
    assert (
        errors
        == "FOO: Please double check your CHROM column.  This column must only be these values: {possible_vals}\n".format(
            possible_vals=", ".join(validate.ACCEPTED_CHROMOSOMES)
        )
    )
    assert warnings == "FOO: Should not have the chr prefix in front of chromosomes.\n"


def test_invalid_nochr__validate_chromosome():
    """Checks that errors and warnings get thrown if
    chromosome col has invalid values and chr is not allowed"""
    input_df = pd.DataFrame(
        {
            "sample_id": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID1-1", "ID3-1"],
            "CHROM": ["chr2", 9, 100],
        }
    )
    errors, warnings = validate._validate_chromosome(
        df=input_df, col="CHROM", fileformat="FOO", allow_chr=False
    )
    assert (
        errors == "FOO: Should not have the chr prefix in front of chromosomes.\n"
        "FOO: Please double check your CHROM column.  This column must only be these values: {possible_vals}\n".format(
            possible_vals=", ".join(validate.ACCEPTED_CHROMOSOMES)
        )
    )
    assert warnings == "", "Warnings should be empty"


@pytest.mark.parametrize(
    "test_na_allowed,expected_val",
    [(True, True), (False, False)],
    ids=[
        "allow_na_is_true",
        "allow_na_is_false",
    ],
)
def test_that__validate_chromosome_calls_check_col_and_values_with_correct_na_allowed_val(
    test_na_allowed, expected_val
):
    input_df = pd.DataFrame({"SITE1_CHROMOSOME": [2, 3, 4]})
    with patch.object(
        process_functions, "check_col_and_values", return_value=("", "")
    ) as check_col_and_values_mock:
        validate._validate_chromosome(
            df=input_df,
            col="SITE1_CHROMOSOME",
            fileformat="Structural Variant",
            allow_na=test_na_allowed,
        )
        check_col_and_values_mock.assert_called_once_with(
            df=input_df,
            col="SITE1_CHROMOSOME",
            possible_values=validate.ACCEPTED_CHROMOSOMES,
            filename="Structural Variant",
            na_allowed=expected_val,
        )


ONCOTREE_ENT = "syn222"


class argparser:
    oncotree_link = "link"
    parentid = "syn3333"
    filetype = None
    center = "try"
    filepath = ["path.csv"]
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


def test_perform_validate(syn, genie_config):
    """Make sure all functions are called"""
    arg = argparser()
    valid = True
    with patch.object(
        validate, "_check_parentid_permission_container", return_value=None
    ) as patch_check_parentid, patch.object(
        extract, "get_genie_config", return_value=genie_config
    ) as patch_get_config, patch.object(
        validate, "_check_center_input"
    ) as patch_check_center, patch.object(
        extract, "_get_oncotreelink"
    ) as patch_get_onco, patch.object(
        validate.GenieValidationHelper,
        "validate_single_file",
        return_value=(valid, "foo"),
    ) as patch_validate, patch.object(
        load, "store_files"
    ) as patch_syn_upload:
        validate._perform_validate(syn, arg)
        patch_check_parentid.assert_called_once_with(syn=syn, parentid=arg.parentid)
        patch_get_config.assert_called_once_with(syn=syn, project_id=arg.project_id)
        patch_check_center.assert_called_once_with(arg.center, ["SAGE", "TEST", "GOLD"])
        patch_get_onco.assert_called_once()
        patch_validate.assert_called_once_with(
            nosymbol_check=arg.nosymbol_check, project_id=arg.project_id
        )
        patch_syn_upload.assert_called_once_with(
            syn=syn, filepaths=arg.filepath, parentid=arg.parentid
        )


@pytest.mark.parametrize(
    "test_nested_list,expected",
    [
        ([[]], {"files": {}, "file_info": {"name": "", "path": []}}),
        (
            [
                [{"name": "test_txt1", "path": "/some_path1.txt"}],
                [{"name": "test_txt2", "path": "/some_path2.txt"}],
            ],
            {"files": {}, "file_info": {"name": "", "path": []}},
        ),
        (
            [
                [
                    {"name": "test_file1_1", "path": "/some_path1_1.txt"},
                    {"name": "test_file1_2", "path": "/some_path1_2.txt"},
                ],
                [{"name": "test_file2", "path": "/some_path2.txt"}],
            ],
            {
                "files": {
                    "test_file1_1": "/some_path1_1.txt",
                    "test_file1_2": "/some_path1_2.txt",
                },
                "file_info": {
                    "name": "test_file1_1,test_file1_2",
                    "path": ["/some_path1_1.txt", "/some_path1_2.txt"],
                },
            },
        ),
    ],
    ids=["empty_nested_list", "no_files_with_substr", "valid_files"],
)
def test_that_parse_file_info_in_nested_list_returns_expected(
    test_nested_list, expected
):
    result = validate.parse_file_info_in_nested_list(
        nested_list=test_nested_list, search_str="test_file1"
    )
    assert result == expected


def test_that_parse_file_info_in_nested_list_returns_expected_with_ignore_case():
    test_nested_list = [
        [
            {"name": "TEST-file1-1", "path": "/some_path1_1.txt"},
            {"name": "test-file1-2", "path": "/some_path1_2.txt"},
        ],
        [{"name": "test_file2", "path": "/some_path2.txt"}],
    ]
    expected = {
        "files": {
            "TEST-file1-1": "/some_path1_1.txt",
            "test-file1-2": "/some_path1_2.txt",
        },
        "file_info": {
            "name": "TEST-file1-1,test-file1-2",
            "path": ["/some_path1_1.txt", "/some_path1_2.txt"],
        },
    }
    result = validate.parse_file_info_in_nested_list(
        nested_list=test_nested_list, search_str="tEsT-file1", ignore_case=True
    )
    assert result == expected


def test_that_parse_file_info_in_nested_list_returns_expected_with_allow_underscore():
    test_nested_list = [
        [
            {"name": "test_file1_1", "path": "/some_path1_1.txt"},
            {"name": "test-file1_2", "path": "/some_path1_2.txt"},
        ],
        [{"name": "test_file2", "path": "/some_path2.txt"}],
    ]
    expected = {
        "files": {
            "test_file1_1": "/some_path1_1.txt",
            "test-file1_2": "/some_path1_2.txt",
        },
        "file_info": {
            "name": "test_file1_1,test-file1_2",
            "path": ["/some_path1_1.txt", "/some_path1_2.txt"],
        },
    }
    result = validate.parse_file_info_in_nested_list(
        nested_list=test_nested_list, search_str="test_file1", allow_underscore=True
    )
    assert result == expected


def test_that_parse_file_info_in_nested_list_returns_expected_with_ignore_case_and_allow_underscore():
    test_nested_list = [
        [
            {"name": "TEST-file1_1", "path": "/some_path1_1.txt"},
            {"name": "test_fiLE1_2", "path": "/some_path1_2.txt"},
        ],
        [{"name": "test_file2", "path": "/some_path2.txt"}],
    ]
    expected = {
        "files": {
            "TEST-file1_1": "/some_path1_1.txt",
            "test_fiLE1_2": "/some_path1_2.txt",
        },
        "file_info": {
            "name": "TEST-file1_1,test_fiLE1_2",
            "path": ["/some_path1_1.txt", "/some_path1_2.txt"],
        },
    }
    result = validate.parse_file_info_in_nested_list(
        nested_list=test_nested_list,
        search_str="TEST_file1",
        ignore_case=True,
        allow_underscore=True,
    )
    assert result == expected


@pytest.mark.parametrize(
    "test_df1,test_df2,expected_errors,expected_warnings,ignore_case,allow_underscore",
    [
        (
            pd.DataFrame({"ID": [1, 2]}),
            pd.DataFrame({"ID2": [1, 2, 3]}),
            "",
            "",
            False,
            False,
        ),
        (
            pd.DataFrame({"ID": [1, 2, 3, 4]}),
            pd.DataFrame({"ID2": [1, 2, 3]}),
            "At least one ID in your test1 file does not exist as a ID2 in your test2 file. "
            "Please update your file(s) to be consistent.\n",
            "",
            False,
            False,
        ),
        (
            pd.DataFrame({"ID": [3, 4, 5]}),
            pd.DataFrame({"ID2": [1, 2, 3]}),
            "At least one ID in your test1 file does not exist as a ID2 in your test2 file. "
            "Please update your file(s) to be consistent.\n",
            "",
            False,
            False,
        ),
        (
            pd.DataFrame({"ID": ["TEST1", "TEST2", "TeSt3"]}),
            pd.DataFrame({"ID2": ["test1", "test2", "TEST3"]}),
            "",
            "",
            True,
            False,
        ),
        (
            pd.DataFrame({"ID": ["TEST_1", "TEST_2", "TEST-3"]}),
            pd.DataFrame({"ID2": ["TEST_1", "TEST-2", "TEST_3"]}),
            "",
            "",
            False,
            True,
        ),
        (
            pd.DataFrame({"ID": ["TEST_1", "TEST_2", "TeSt-3"]}),
            pd.DataFrame({"ID2": ["test_1", "test-2", "TEST_3"]}),
            "",
            "",
            True,
            True,
        ),
        (
            pd.DataFrame({"ID": ["TEST__1", "TEST_2", "TeSt-3"]}),
            pd.DataFrame({"ID2": ["test_1", "test-2", "TEST_3"]}),
            "At least one ID in your test1 file does not exist as a ID2 in your test2 file. "
            "Please update your file(s) to be consistent.\n",
            "",
            True,
            True,
        ),
        (
            pd.DataFrame({"ID": ["TEST1", "TEST2", "TeSt3"]}),
            pd.DataFrame({"ID2": ["test1", "test2", "TEST3"]}),
            "At least one ID in your test1 file does not exist as a ID2 in your test2 file. "
            "Please update your file(s) to be consistent.\n",
            "",
            False,
            False,
        ),
    ],
    ids=[
        "all_match",
        "some_match",
        "no_match",
        "ignore_case",
        "allow_underscore",
        "ignore_case_and_allow_underscore",
        "str_some_match",
        "str_no_match",
    ],
)
def test_that_check_values_between_two_df_returns_expected(
    test_df1,
    test_df2,
    expected_errors,
    expected_warnings,
    ignore_case,
    allow_underscore,
):
    errors, warnings = validate.check_values_between_two_df(
        df1=test_df1,
        df1_filename="test1",
        df1_id_to_check="ID",
        df2=test_df2,
        df2_filename="test2",
        df2_id_to_check="ID2",
        ignore_case=ignore_case,
        allow_underscore=allow_underscore,
    )
    assert errors == expected_errors
    assert warnings == expected_warnings


@pytest.mark.parametrize(
    "test_input,expected_errors,expected_warnings",
    [
        (
            pd.DataFrame(
                {"start_pos": [1, 2], "end_pos": [3, 4], "some_col": ["a", "b"]}
            ),
            "",
            "",
        ),
        (
            pd.DataFrame(
                {"start_pos": [1, 2], "end_pos": [0, 4], "some_col": ["a", "b"]}
            ),
            "",
            (
                "test_file: Your variants file has record(s) that have an end position "
                "value less than the start position value. Please update your file to be consistent. "
                "When we annotate using the genome-nexus-annotation-pipeline, the records with this "
                "position discrepancy will be re-annotated with different end position values.\n"
            ),
        ),
        (
            pd.DataFrame(
                {"start_pos": [1, 2], "end_pos": [0, 1], "some_col": ["a", "b"]}
            ),
            "",
            (
                "test_file: Your variants file has record(s) that have an end position "
                "value less than the start position value. Please update your file to be consistent. "
                "When we annotate using the genome-nexus-annotation-pipeline, the records with this "
                "position discrepancy will be re-annotated with different end position values.\n"
            ),
        ),
        (
            pd.DataFrame(
                {"start_pos": [1, 2], "end_pos": [1, 2], "some_col": ["a", "b"]}
            ),
            "",
            "",
        ),
    ],
    ids=["start_lt_end_pos", "end_lt_start_pos", "end_all_lt_start_pos", "equal_pos"],
)
def test_that_check_variant_start_and_end_positions_returns_expected(
    test_input, expected_errors, expected_warnings
):
    errors, warnings = validate.check_variant_start_and_end_positions(
        input_df=test_input,
        start_pos_col="start_pos",
        end_pos_col="end_pos",
        filename="test_file",
    )
    assert errors == expected_errors
    assert warnings == expected_warnings


@pytest.mark.parametrize(
    "input_str,expected,ignore_case,allow_underscore",
    [
        ("SAGe-1", "sage-1", True, False),
        ("SAGe_1", "SAGe-1", False, True),
        ("SAGe_1", "sage-1", True, True),
        (120, 120, True, True),
    ],
    ids=[
        "ignore_case",
        "allow_underscore",
        "ignore_case_and_allow_underscore",
        "non_str",
    ],
)
def test_that_standardize_string_for_validation_returns_expected(
    input_str, expected, ignore_case, allow_underscore
):
    test_str = validate.standardize_string_for_validation(
        input_string=input_str,
        ignore_case=ignore_case,
        allow_underscore=allow_underscore,
    )
    assert test_str == expected


def get_invalid_allele_rows_test_cases():
    return [
        {
            "name": "correct_alleles",
            "input": pd.DataFrame(
                {
                    "REFERENCE_ALLELE": [
                        "NANANANA",
                        "ACGTN",
                        "A",
                        "C",
                        "T",
                        "G",
                        "-",
                        "N",
                    ]
                }
            ),
            "expected_index": pd.Index([]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "N"],
            "allowed_ind_alleles": ["-"],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "incorrect_alleles",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["@##", "ACGTX", "XXX"]}),
            "expected_index": pd.Index([0, 1, 2]),
            "allowed_comb_alleles": ["A", "T", "C", "G"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "case_ignored",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["acgtg", "acgt", "-", "a"]}),
            "expected_index": pd.Index([]),
            "allowed_comb_alleles": ["A", "T", "C", "G"],
            "allowed_ind_alleles": ["-"],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "case_not_ignored",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["acgt-G"]}),
            "expected_index": pd.Index([0]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "-"],
            "allowed_ind_alleles": [],
            "ignore_case": False,
            "allow_na": True,
        },
        {
            "name": "no_ind_alleles_incorrect",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["ACG-T", "ACGT", "G-CT"]}),
            "expected_index": pd.Index([0, 2]),
            "allowed_comb_alleles": ["A", "T", "C", "G"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "no_ind_alleles_correct",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["ACT", "ACGT", "G"]}),
            "expected_index": pd.Index([]),
            "allowed_comb_alleles": ["A", "T", "C", "G"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "missing_entries_not_allowed",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["ACGT-G", pd.NA, None]}),
            "expected_index": pd.Index([1, 2]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "-"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": False,
        },
        {
            "name": "missing_entries_allowed",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["ACGT-G", pd.NA, None]}),
            "expected_index": pd.Index([]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "-"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "no_specified_alleles_values",
            "input": pd.DataFrame({"REFERENCE_ALLELE": ["ACGT-G", "ACGF", "B"]}),
            "expected_index": pd.Index([]),
            "allowed_comb_alleles": [],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "float_nas_not_allowed",
            "input": pd.DataFrame(
                {"REFERENCE_ALLELE": [1.5, 2.0, float("nan"), 3.5, 4.0]}
            ),
            "expected_index": pd.Index([0, 1, 2, 3, 4]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "-"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": False,
        },
        {
            "name": "float_nas_allowed",
            "input": pd.DataFrame(
                {"REFERENCE_ALLELE": [1.5, 2.0, float("nan"), 3.5, 4.0]}
            ),
            "expected_index": pd.Index([0, 1, 3, 4]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "-"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "all_missing_nas_allowed",
            "input": pd.DataFrame(
                {"REFERENCE_ALLELE": [float("nan"), float("nan"), float("nan")]}
            ),
            "expected_index": pd.Index([]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "-"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": True,
        },
        {
            "name": "all_missing_nas_not_allowed",
            "input": pd.DataFrame(
                {"REFERENCE_ALLELE": [float("nan"), float("nan"), float("nan")]}
            ),
            "expected_index": pd.Index([0, 1, 2]),
            "allowed_comb_alleles": ["A", "T", "C", "G", "-"],
            "allowed_ind_alleles": [],
            "ignore_case": True,
            "allow_na": False,
        },
    ]


@pytest.mark.parametrize(
    "test_cases", get_invalid_allele_rows_test_cases(), ids=lambda x: x["name"]
)
def test_that_get_invalid_allele_rows_returns_expected(test_cases):
    invalid_rows = validate.get_invalid_allele_rows(
        test_cases["input"],
        input_col="REFERENCE_ALLELE",
        allowed_comb_alleles=test_cases["allowed_comb_alleles"],
        allowed_ind_alleles=test_cases["allowed_ind_alleles"],
        ignore_case=test_cases["ignore_case"],
        allow_na=test_cases["allow_na"],
    )
    assert invalid_rows.equals(test_cases["expected_index"])


def get_allele_validation_message_test_cases():
    return [
        {
            "name": "has_invalid_alleles",
            "input_invalid_rows": pd.Index([1, 2, 3]),
            "allowed_comb_alleles": ["A", "C", "T", "G"],
            "allowed_ind_alleles": ["-"],
            "expected_error": (
                "maf: Your REFERENCE_ALLELE column has invalid allele values. "
                "This is the list of accepted allele values that can appear individually "
                "or in combination with each other: A,C,T,G.\n"
                "This is the list of accepted allele values that can only appear individually: -\n"
            ),
            "expected_warning": "",
        },
        {
            "name": "has_no_invalid_alleles",
            "input_invalid_rows": [],
            "allowed_comb_alleles": [],
            "allowed_ind_alleles": [],
            "expected_error": "",
            "expected_warning": "",
        },
        {
            "name": "has_invalid_alleles_empty_ind_alleles",
            "input_invalid_rows": pd.Index([1, 2, 3]),
            "allowed_comb_alleles": ["A", "C", "T", "G"],
            "allowed_ind_alleles": [],
            "expected_error": (
                "maf: Your REFERENCE_ALLELE column has invalid allele values. "
                "This is the list of accepted allele values that can appear individually "
                "or in combination with each other: A,C,T,G.\n"
                "This is the list of accepted allele values that can only appear individually: \n"
            ),
            "expected_warning": "",
        },
    ]


@pytest.mark.parametrize(
    "test_cases", get_allele_validation_message_test_cases(), ids=lambda x: x["name"]
)
def test_that_get_allele_validation_message_returns_expected(test_cases):
    error, warning = validate.get_allele_validation_message(
        invalid_indices=test_cases["input_invalid_rows"],
        invalid_col="REFERENCE_ALLELE",
        allowed_comb_alleles=test_cases["allowed_comb_alleles"],
        allowed_ind_alleles=test_cases["allowed_ind_alleles"],
        fileformat="maf",
    )
    assert error == test_cases["expected_error"]
    assert warning == test_cases["expected_warning"]
