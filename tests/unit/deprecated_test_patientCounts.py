import synapseclient
import pandas as pd
import mock
import pytest
from genie import patientCounts

json_oncotreeurl = \
    "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"

syn = mock.create_autospec(synapseclient.Synapse)

pc = patientCounts(syn, "SAGE")

onco_map_dict = {
  "AMPCA": {
    'CANCER_TYPE': "Ampullary Cancer",
    'CANCER_TYPE_DETAILED': "Ampullary Carcinoma",
    'ONCOTREE_PRIMARY_NODE': "AMPULLA_OF_VATER",
    'ONCOTREE_SECONDARY_NODE': "AMPCA"},
  "TESTIS": {
    'CANCER_TYPE': "Testicular Cancer, NOS",
    'CANCER_TYPE_DETAILED': "Testis",
    'ONCOTREE_PRIMARY_NODE': "TESTIS",
    'ONCOTREE_SECONDARY_NODE': ''},
  "UCEC": {
    'CANCER_TYPE': "Endometrial Cancer",
    'CANCER_TYPE_DETAILED': "Endometrial Carcinoma",
    'ONCOTREE_PRIMARY_NODE': "UTERUS",
    'ONCOTREE_SECONDARY_NODE': "UCEC"}}


def test_processing():
    expectedpcDf = pd.DataFrame(dict(
        CENTER=['SAGE', 'SAGE', 'SAGE'],
        ONCOTREE_CODE=['AMPCA', 'UCEC', 'TESTIS'],
        NUM_PATIENTS_PD1_PDL1=[1, 2, 3],
        PRIMARY_CODE=['AMPULLA_OF_VATER', 'UTERUS', 'TESTIS']))

    pcDf = pd.DataFrame(dict(
        CENTER=['foo', 'foo', 'foo'],
        ONCOTREE_CODE=['AMPCA', 'UCEC', 'TESTIS'],
        NUM_PATIENTS_PD1_PDL1=[1, 2, 3]))
    with mock.patch(
            "genie.process_functions.get_oncotree_code_mappings",
            return_value=onco_map_dict) as mock_get_onco_map:
        newpcDf = pc._process(pcDf, json_oncotreeurl)
        mock_get_onco_map.called_once_with(json_oncotreeurl)
        assert expectedpcDf.equals(newpcDf[expectedpcDf.columns])


def test_perfect_validateFilename():
    assert pc.validateFilename(["patient_counts.txt"]) == "patientCounts"


def test_invalid_validateFilename():
    with pytest.raises(AssertionError):
        pc.validateFilename(["foo"])


def test_perfect__validate():
    pcDf = pd.DataFrame(dict(
        CENTER=['foo', 'foo', 'foo'],
        ONCOTREE_CODE=['AMPCA', 'UCEC', 'TESTIS'],
        NUM_PATIENTS_PD1_PDL1=[1, 2, 3]))
    with mock.patch(
            "genie.process_functions.get_oncotree_code_mappings",
            return_value=onco_map_dict) as mock_get_onco_map:
        error, warning = pc._validate(pcDf, json_oncotreeurl)
        mock_get_onco_map.called_once_with(json_oncotreeurl)
        assert error == ""
        assert warning == ""


def test_invalid__validate():
    pcDf = pd.DataFrame(dict(
        ONCOTREE_CODE=['AMPCA', 'AMPCA', 'CHOLD', 'AMPCA', 'TESTIS'],
        NUM_PATIENTS_PD1_PDL1=[1, 2, 3, 4, float('nan')]))
    with mock.patch(
            "genie.process_functions.get_oncotree_code_mappings",
            return_value=onco_map_dict):
        error, warning = pc._validate(pcDf, json_oncotreeurl)
        expectedErrors = (
            "Patient Counts: Must not have any duplicated ONCOTREE CODES.\n"
            "Patient Counts: Please double check that all your ONCOTREE CODES "
            "exist in the mapping. You have 1 codes that don't map. These are "
            "the codes that don't map: CHOLD\n"
            "Patient Counts: Must not have any null values, "
            "and must be all integers.\n")
        assert error == expectedErrors
        assert warning == ""


def test_missingcol__validate():
    '''
    Test for missing columns
    '''
    pcdf = pd.DataFrame()
    with mock.patch(
            "genie.process_functions.get_oncotree_code_mappings",
            return_value=onco_map_dict):
        error, warning = pc._validate(pcdf, json_oncotreeurl)
        expectedErrors = (
            "Patient Counts: File must have ONCOTREE_CODE column.\n"
            "Patient Counts: File must have NUM_PATIENTS_PD1_PDL1 column.\n")
        assert error == expectedErrors
        assert warning == ""
