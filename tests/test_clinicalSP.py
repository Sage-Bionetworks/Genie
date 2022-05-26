from unittest import mock

import pandas as pd
import pytest
import synapseclient

from genie_registry.clinicalSP import clinicalSP

syn = mock.create_autospec(synapseclient.Synapse)

clinsp_class = clinicalSP(syn, "SAGE")


def test_filetype():
    assert clinsp_class._fileType == "clinicalSP"


def test_incorrect_validatefilename():
    with pytest.raises(AssertionError):
        clinsp_class.validateFilename(["foo"])


def test_correct_validatefilename():
    assert clinsp_class.validateFilename(["nonGENIE_data_clinical.txt"]) == "clinicalSP"


def test_perfect__processing():
    expected_clindf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            SAMPLE_ID=[
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
            SEQ_ASSAY_ID=[
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
            ],
        )
    )
    # TEST patient processing
    clindf = pd.DataFrame(
        dict(
            PATIENT_ID=["ID1", "ID2", "ID3", "ID4", "ID5"],
            SAMPLE_ID=["ID1-1", "ID2-1", "ID3-1", "ID4-1", "ID5-1"],
            SEQ_ASSAY_ID=[
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
            ],
        )
    )

    clindf = clinsp_class._process(clindf)
    assert expected_clindf.equals(clindf[expected_clindf.columns])


def test_perfect__validate():
    clindf = pd.DataFrame(
        dict(
            PATIENT_ID=["ID1", "ID2", "ID3", "ID4", "ID5"],
            SAMPLE_ID=["ID1-1", "ID2-1", "ID3-1", "ID4-1", "ID5-1"],
            SEQ_ASSAY_ID=[
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
            ],
        )
    )

    error, warning = clinsp_class._validate(clindf)
    assert error == ""
    assert warning == ""


def test_missingcols__validate():
    clindf = pd.DataFrame()
    error, warning = clinsp_class._validate(clindf)
    expected_errors = (
        "nonGENIE_data_clinical.txt: File must have SAMPLE_ID column.\n"
        "nonGENIE_data_clinical.txt: File must have SEQ_ASSAY_ID column.\n"
        "nonGENIE_data_clinical.txt: File must have PATIENT_ID column.\n"
    )
    assert error == expected_errors
    assert warning == ""


def test_errors__validate():
    clindf = pd.DataFrame(
        dict(
            PATIENT_ID=[float("nan"), "ID2", "ID3", "ID4", "ID5"],
            SAMPLE_ID=["ID1-1", float("nan"), "ID3-1", "ID4-1", "ID5-1"],
            SEQ_ASSAY_ID=[
                float("nan"),
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
            ],
        )
    )

    error, warning = clinsp_class._validate(clindf)

    expected_errors = (
        "nonGENIE_data_clinical.txt: "
        "There can't be any blank values for SAMPLE_ID\n"
        "nonGENIE_data_clinical.txt: "
        "There can't be any blank values for PATIENT_ID\n"
    )
    expected_warnings = (
        "nonGENIE_data_clinical.txt: "
        "Please double check your SEQ_ASSAY_ID columns, "
        "there are empty rows.\n"
    )

    assert error == expected_errors
    assert warning == expected_warnings
