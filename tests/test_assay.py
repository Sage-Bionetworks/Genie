"""Test assay information validation and processing"""

import copy
from unittest.mock import patch

import pandas as pd
import pytest
from genie import extract, process_functions
from genie_registry.assay import Assayinfo

GDC_DATA_DICT = {
    "properties": {
        "library_strategy": {"enum": ["value1", "value2"]},
        "library_selection": {"enum": ["value1", "value2"]},
        "platform": {"enum": ["value1", "value2"]},
        "instrument_model": {"enum": ["value1", "value2"]},
        "target_capture_kit": {"enum": ["value1", "value2"]},
    }
}


@pytest.fixture()
def assay_info(syn):
    return Assayinfo(syn, "SAGE", genie_config={"sample": "syn1234"})


def test_filetype(assay_info):
    assert assay_info._fileType == "assayinfo"


def test_invalidname__validatefilename(assay_info):
    """Assertion error thrown for wrong filename"""
    with pytest.raises(AssertionError):
        assay_info.validateFilename(["foo"])


def test_correct__validatefilename(assay_info):
    """Filetype returned when filename valid"""
    assert assay_info.validateFilename(["assay_information.yaml"]) == "assayinfo"


def test_validinput__validate(assay_info):
    """Valid input should have no errors or warnings"""
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["SAGE-1", "SAGE-3"],
        "is_paired_end": [True, False],
        "library_strategy": ["Targeted Sequencing", "WXS"],
        "library_selection": ["value1", "value2"],
        "platform": ["value1", "value2"],
        "instrument_model": ["value1", "value2"],
        "target_capture_kit": ["value1", "value2"],
        "variant_classifications": ["Frame_Shift_Ins", "Frame_Shift_Ins"],
        "read_length": [22, float("nan")],
        "number_of_genes": [5, 20],
        "gene_padding": [10, None],
        "calling_strategy": ["tumor_only", "tumor_normal"],
        "specimen_tumor_cellularity": [">10%", ">20%"],
        "alteration_types": ["snv;small_indels", "intragenic_cna"],
        "preservation_technique": ["FFPE", "FFPE;fresh_frozen"],
        "coverage": ["hotspot_regions;introns", "introns"],
    }
    uniq_seq_df = pd.DataFrame({"seq": ["SAGE-1", "SAGE-3"]})
    assay_info_df = pd.DataFrame(assay_info_dict)
    test_dict = copy.deepcopy(GDC_DATA_DICT)
    with patch.object(
        extract, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = assay_info._validate(assay_info_df, skip_database_checks=False)
        assert error == ""
        assert warning == ""
        patch_get_gdc.assert_called()


def test_case__validate(assay_info):
    """Valid input with lowercase SEQ_ASSAY_ID, should have no errors or warnings"""
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["sage-1", "SAGE-3"],
        "is_paired_end": [True, False],
        "library_strategy": ["Targeted Sequencing", "WXS"],
        "library_selection": ["value1", "value2"],
        "platform": ["value1", "value2"],
        "instrument_model": ["value1", "value2"],
        "target_capture_kit": ["value1", "value2"],
        "variant_classifications": ["Frame_Shift_Ins", "Frame_Shift_Ins"],
        "read_length": [22, float("nan")],
        "number_of_genes": [5, 20],
        "gene_padding": [10, None],
        "calling_strategy": ["tumor_only", "tumor_normal"],
        "specimen_tumor_cellularity": [">10%", ">20%"],
        "alteration_types": ["snv;small_indels", "intragenic_cna"],
        "preservation_technique": ["FFPE", "FFPE;fresh_frozen"],
        "coverage": ["hotspot_regions;introns", "introns"],
    }
    uniq_seq_df = pd.DataFrame({"seq": ["SAGE-1", "SAGE-3"]})
    assay_info_df = pd.DataFrame(assay_info_dict)
    test_dict = copy.deepcopy(GDC_DATA_DICT)
    with patch.object(
        extract, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = assay_info._validate(assay_info_df, skip_database_checks=False)
        assert error == ""
        assert warning == ""
        patch_get_gdc.assert_called()


def test_underscore__validate(assay_info):
    """Valid input with underscore in SEQ_ASSAY_ID, should have no errors or warnings"""
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["SAGE_1", "SAGE-3"],
        "is_paired_end": [True, False],
        "library_strategy": ["Targeted Sequencing", "WXS"],
        "library_selection": ["value1", "value2"],
        "platform": ["value1", "value2"],
        "instrument_model": ["value1", "value2"],
        "target_capture_kit": ["value1", "value2"],
        "variant_classifications": ["Frame_Shift_Ins", "Frame_Shift_Ins"],
        "read_length": [22, float("nan")],
        "number_of_genes": [5, 20],
        "gene_padding": [10, None],
        "calling_strategy": ["tumor_only", "tumor_normal"],
        "specimen_tumor_cellularity": [">10%", ">20%"],
        "alteration_types": ["snv;small_indels", "intragenic_cna"],
        "preservation_technique": ["FFPE", "FFPE;fresh_frozen"],
        "coverage": ["hotspot_regions;introns", "introns"],
    }
    uniq_seq_df = pd.DataFrame({"seq": ["SAGE-1", "SAGE-3"]})
    assay_info_df = pd.DataFrame(assay_info_dict)
    test_dict = copy.deepcopy(GDC_DATA_DICT)
    with patch.object(
        extract, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = assay_info._validate(assay_info_df, skip_database_checks=False)
        assert error == ""
        assert warning == ""
        patch_get_gdc.assert_called()


def test__missingcols__validate(assay_info):
    """Test missing columns"""
    assay_info_df = pd.DataFrame()
    test_dict = copy.deepcopy(GDC_DATA_DICT)
    with patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = assay_info._validate(assay_info_df, skip_database_checks=False)
    expected_errors = (
        "Assay_information.yaml: Must have SEQ_ASSAY_ID column.\n"
        "Assay_information.yaml: Must have is_paired_end column.\n"
        "Assay_information.yaml: Must have library_selection column.\n"
        "Assay_information.yaml: Must have library_strategy column.\n"
        "Assay_information.yaml: Must have platform column.\n"
        "Assay_information.yaml: Must have instrument_model column.\n"
        "Assay_information.yaml: Must have target_capture_kit column.\n"
        "Assay_information.yaml: Must have read_length column.\n"
        "Assay_information.yaml: Must have number_of_genes column.\n"
        "Assay_information.yaml: Must have calling_strategy column.\n"
        "Assay_information.yaml: "
        "Must have specimen_tumor_cellularity column.\n"
        "Assay_information.yaml: Must have alteration_types column.\n"
        "Assay_information.yaml: Must have preservation_technique column.\n"
        "Assay_information.yaml: Must have coverage column.\n"
    )
    assert error == expected_errors

    expected_warnings = (
        "Assay_information.yaml: Doesn't have variant_classifications column. "
        "This column will be added\n"
        "Assay_information.yaml: gene_padding is "
        "by default 10 if not specified.\n"
    )
    assert warning == expected_warnings
    patch_get_gdc.assert_called()


def test_fillcols__process(assay_info):
    """
    Standardization of SEQ_ASSAY_ID
    Add in CENTER, gene_padding, and variant_classifications if missing
    """

    assay_info_dict = {"SEQ_ASSAY_ID": ["SAGE-Foo_1"], "SEQ_PIPELINE_ID": ["SAGE_Foo"]}
    assay_info_df = pd.DataFrame(assay_info_dict)
    processed_assay_df = assay_info._process(assay_info_df)
    expected_assay_df = pd.DataFrame(
        {
            "SEQ_ASSAY_ID": ["SAGE-FOO-1"],
            "SEQ_PIPELINE_ID": ["SAGE-FOO"],
            "gene_padding": [10],
            "variant_classifications": [float("nan")],
            "CENTER": ["SAGE"],
        }
    )

    assert expected_assay_df.equals(processed_assay_df[expected_assay_df.columns])


def test_default10__process(assay_info):
    """
    gene_padding default 10
    """

    assay_info_dict = {
        "SEQ_ASSAY_ID": ["SAGE-1", "SAGE-2"],
        "SEQ_PIPELINE_ID": ["SAGE-1", "SAGE-2"],
        "gene_padding": [20, float("nan")],
        "variant_classifications": ["test", "test"],
    }
    assay_info_df = pd.DataFrame(assay_info_dict)
    processed_assay_df = assay_info._process(assay_info_df)
    expected_assay_df = pd.DataFrame(
        {
            "SEQ_ASSAY_ID": ["SAGE-1", "SAGE-2"],
            "SEQ_PIPELINE_ID": ["SAGE-1", "SAGE-2"],
            "gene_padding": [20, 10],
            "variant_classifications": ["test", "test"],
            "CENTER": ["SAGE", "SAGE"],
        }
    )
    assert expected_assay_df.equals(processed_assay_df[expected_assay_df.columns])


def test_invalid__validate(assay_info):
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["SAGE-1", "SAG-2"],
        "is_paired_end": [True, "foo"],
        "library_strategy": ["foo", "WXS"],
        "library_selection": ["foo", "PCR"],
        "platform": ["foo", "Illumina"],
        "instrument_model": ["foo", "Illumina HiSeq 4000"],
        "variant_classifications": ["foo", "Frame_Shift_Ins"],
        "target_capture_kit": ["foo", "doo"],
        "read_length": [22, "foo"],
        "number_of_genes": [5, "foo"],
        "gene_padding": [10, "foo"],
        "calling_strategy": ["tumor_ony", "tumor_normal"],
        "specimen_tumor_cellularity": [">10", ">20%"],
        "alteration_types": ["snv;small_indel", "intragenic_cna"],
        "preservation_technique": ["FPE", "FFPE;fresh_frozen"],
        "coverage": ["hotsot_regions;introns", "introns"],
    }
    assay_info_df = pd.DataFrame(assay_info_dict)
    uniq_seq_df = pd.DataFrame({"seq": ["SAGE-1", "SAGE-2"]})

    # Must use deepcopy because a dict.copy is a shallow copy
    # Which just points to reference keys
    test_dict = copy.deepcopy(GDC_DATA_DICT)
    with patch.object(
        extract, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = assay_info._validate(assay_info_df, skip_database_checks=False)
        expected_errors = (
            "Assay_information.yaml: "
            "Please make sure all your SEQ_ASSAY_IDs start with your "
            "center abbreviation.\n"
            "Assay_information.yaml: You are missing SEQ_ASSAY_IDs: SAGE-2\n"
            "Assay_information.yaml: "
            "Please double check your is_paired_end column.  "
            "This column must only be these values: True, False\n"
            "Assay_information.yaml: "
            "Please double check your library_selection column.  "
            "This column must only be these values: value1, value2\n"
            "Assay_information.yaml: "
            "Please double check your library_strategy column.  "
            "This column must only be these values: Targeted Sequencing, WXS\n"
            "Assay_information.yaml: "
            "Please double check your platform column.  "
            "This column must only be these values: value1, value2\n"
            "Assay_information.yaml: "
            "Please double check your instrument_model column.  "
            "This column must only be these values: value1, value2, "
            "Illumina NovaSeq 6000, None\n"
            "Assay_information.yaml: "
            "Please double check your variant_classifications column.  "
            "This column must only be these values: Splice_Site, "
            "Nonsense_Mutation, Frame_Shift_Del, Frame_Shift_Ins, "
            "Nonstop_Mutation, Translation_Start_Site, In_Frame_Ins, "
            "In_Frame_Del, Missense_Mutation, Intron, Splice_Region, "
            "Silent, RNA, 5'UTR, 3'UTR, IGR, 5'Flank, 3'Flank, None\n"
            "Assay_information.yaml: "
            "Please double check your read_length.  "
            "It must be an integer or null.\n"
            "Assay_information.yaml: "
            "Please double check your number_of_genes. "
            "It must be an integer.\n"
            "Assay_information.yaml: "
            "Please double check your gene_padding. "
            "It must be an integer or blank.\n"
            "Assay_information.yaml: "
            "Please double check your calling_strategy column.  This "
            "column must only be these values: tumor_only, tumor_normal, "
            "plasma_normal\n"
            "Assay_information.yaml: "
            "Please double check your specimen_tumor_cellularity. "
            "It must in this format >(num)%. ie. >10%\n"
            "Assay_information.yaml: "
            "Please double check your alteration_types column.  "
            "This column must only be these values: snv, small_indels, "
            "gene_level_cna, intragenic_cna, structural_variants\n"
            "Assay_information.yaml: "
            "Please double check your preservation_technique column.  "
            "This column must only be these values: FFPE, fresh_frozen, NA\n"
            "Assay_information.yaml: "
            "Please double check your coverage column.  "
            "This column must only be these values: "
            "hotspot_regions, coding_exons, introns, promoters\n"
        )

        patch_get_gdc.assert_called_once_with("read_group")
        assert error == expected_errors
        assert warning == ""


@pytest.mark.parametrize(
    "skip_database_checks, all_seq_assays, clinical_db_ids, expected_error, expect_called",
    [
        # 1) skip_database_checks=True -> always "", even if IDs "missing"
        (
            True,
            {"A-1": {}, "B-2": {}},
            ["A-1", "C-3"],
            "",
            False,
        ),
        # 2) skip_database_checks=False AND all assay IDs are present in DB -> ""
        (
            False,
            {"A-1": {}, "B-2": {}},
            ["A-1", "B-2"],
            "",
            True,
        ),
        # 3) skip_database_checks=False AND there are more assay ids than in DB -> ""
        (
            False,
            {"A-1": {}, "B-2": {}, "C-2": {}},
            ["A-1", "B-2"],
            "",
            True,
        ),
        # 4) skip_database_checks=False AND DB has at least one ID not in assay file -> error
        (
            False,
            {"A-1": {}, "B-2": {}},
            ["A-1", "B-2", "C-3"],
            "Assay_information.yaml: You are missing SEQ_ASSAY_IDs: C-3\n",
            True,
        ),
    ],
    ids=[
        "skip_database_check_is_true",
        "seq_assay_ids_match_db",
        "assay_has_extra_seq_assay_ids",
        "assay_has_missing_seq_assay_ids",
    ],
)
def test_validate_all_seq_assay_ids_exist_in_clinical_database_patch_object(
    assay_info,
    skip_database_checks,
    all_seq_assays,
    clinical_db_ids,
    expected_error,
    expect_called,
):
    fake_df = pd.DataFrame({"seq": clinical_db_ids})

    # Patch extract.get_syntabledf right where the function looks it up
    with patch.object(extract, "get_syntabledf", return_value=fake_df) as mocked_get:
        err = assay_info.validate_all_seq_assay_ids_exist_in_clinical_database(
            all_seq_assays=all_seq_assays,
            skip_database_checks=skip_database_checks,
        )

    assert err == expected_error

    if expect_called:
        mocked_get.assert_called_once()
    else:
        mocked_get.assert_not_called()
