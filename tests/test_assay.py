"""Test assay information validation and processing"""
import copy
from unittest.mock import patch, create_autospec

import pandas as pd
import pytest
import synapseclient

from genie_registry.assay import Assayinfo
from genie import process_functions

GDC_DATA_DICT = {
    "properties": {
        "library_strategy": {"enum": ["value1", "value2"]},
        "library_selection": {"enum": ["value1", "value2"]},
        "platform": {"enum": ["value1", "value2"]},
        "instrument_model": {"enum": ["value1", "value2"]},
        "target_capture_kit": {"enum": ["value1", "value2"]},
    }
}

SYN = create_autospec(synapseclient.Synapse)
ASSAY_INFO = Assayinfo(SYN, "SAGE")


def test_filetype():
    assert ASSAY_INFO._fileType == "assayinfo"


def test_invalidname__validatefilename():
    """Assertion error thrown for wrong filename"""
    with pytest.raises(AssertionError):
        ASSAY_INFO.validateFilename(["foo"])


def test_correct__validatefilename():
    """Filetype returned when filename valid"""
    assert ASSAY_INFO.validateFilename(["assay_information.yaml"]) == "assayinfo"


def test_validinput__validate():
    """Valid input should have no errors or warnings"""
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["SAGE-1", "SAGE-3"],
        "is_paired_end": [True, False],
        "library_strategy": ["value1", "value2"],
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
        process_functions, "get_synid_database_mappingdf", return_value="syn123"
    ), patch.object(
        process_functions, "getDatabaseSynId", return_value="syn1234"
    ), patch.object(
        process_functions, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = ASSAY_INFO._validate(assay_info_df, "syn9999")
        assert error == ""
        assert warning == ""
        patch_get_gdc.assert_called()


def test_case__validate():
    """Valid input should have no errors or warnings"""
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["sage-1", "SAGE-3"],
        "is_paired_end": [True, False],
        "library_strategy": ["value1", "value2"],
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
        process_functions, "get_synid_database_mappingdf", return_value="syn123"
    ), patch.object(
        process_functions, "getDatabaseSynId", return_value="syn1234"
    ), patch.object(
        process_functions, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = ASSAY_INFO._validate(assay_info_df, "syn9999")
        assert error == ""
        assert warning == ""
        patch_get_gdc.assert_called()


def test_underscore__validate():
    """Valid input should have no errors or warnings"""
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["SAGE_1", "SAGE-3"],
        "is_paired_end": [True, False],
        "library_strategy": ["value1", "value2"],
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
        process_functions, "get_synid_database_mappingdf", return_value="syn123"
    ), patch.object(
        process_functions, "getDatabaseSynId", return_value="syn1234"
    ), patch.object(
        process_functions, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = ASSAY_INFO._validate(assay_info_df, "syn9999")
        assert error == ""
        assert warning == ""
        patch_get_gdc.assert_called()


def test__missingcols__validate():
    """Test missing columns"""
    assay_info_df = pd.DataFrame()
    test_dict = copy.deepcopy(GDC_DATA_DICT)
    with patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = ASSAY_INFO._validate(assay_info_df, "syn99999")
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


def test_fillcols__process():
    """
    Standardization of SEQ_ASSAY_ID
    Add in CENTER, gene_padding, and variant_classifications if missing
    """

    assay_info_dict = {"SEQ_ASSAY_ID": ["SAGE-Foo_1"], "SEQ_PIPELINE_ID": ["SAGE_Foo"]}
    assay_info_df = pd.DataFrame(assay_info_dict)
    processed_assay_df = ASSAY_INFO._process(assay_info_df)
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


def test_default10__process():
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
    processed_assay_df = ASSAY_INFO._process(assay_info_df)
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


def test_invalid__validate():
    assay_info_dict = {
        "SEQ_ASSAY_ID": ["SAGE-1", "SAG-2"],
        "is_paired_end": [True, "foo"],
        "library_strategy": ["foo", "ChIP-Seq"],
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
        process_functions, "get_synid_database_mappingdf", return_value="syn123"
    ), patch.object(
        process_functions, "getDatabaseSynId", return_value="syn1234"
    ), patch.object(
        process_functions, "get_syntabledf", return_value=uniq_seq_df
    ), patch.object(
        process_functions, "get_gdc_data_dictionary", return_value=test_dict
    ) as patch_get_gdc:
        error, warning = ASSAY_INFO._validate(assay_info_df, "syn9999")
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
            "This column must only be these values: value1, value2\n"
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

        patch_get_gdc.called_once_with("read_group")
        assert error == expected_errors
        assert warning == ""
