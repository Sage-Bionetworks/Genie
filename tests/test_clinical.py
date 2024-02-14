from collections import Counter
import datetime
import json
from unittest import mock

import pandas as pd
import pytest
import synapseclient

from genie import process_functions, validate
import genie_registry
from genie_registry.clinical import Clinical


def createMockTable(dataframe):
    table = mock.create_autospec(synapseclient.table.CsvFileTable)
    table.asDataFrame.return_value = dataframe
    return table


def table_query_results(*args):
    return table_query_results_map[args]


no_nan = pd.DataFrame(
    dict(
        CODE=[1, 2, 3, 4, 99],
        CBIO_LABEL=["Test", "Why", "foo", "Me", "Unknown"],
        DESCRIPTION=["non", "asdf", "asdf", "asdff", "asdfasdf"],
    )
)

sexdf = pd.DataFrame(
    dict(
        CODE=[1, 2, 99],
        CBIO_LABEL=["Male", "Female", "Unknown"],
        DESCRIPTION=["Male", "Female", "Not coded"],
    )
)

table_query_results_map = {
    ("select * from syn7434222",): createMockTable(sexdf),
    ("select * from syn7434236",): createMockTable(no_nan),
    ("select * from syn7434242",): createMockTable(no_nan),
    ("select * from syn7434273",): createMockTable(no_nan),
}

json_oncotreeurl = (
    "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"
)

onco_map_dict = {
    "AMPCA": {
        "CANCER_TYPE": "Ampullary Cancer",
        "CANCER_TYPE_DETAILED": "Ampullary Carcinoma",
        "ONCOTREE_PRIMARY_NODE": "AMPULLA_OF_VATER",
        "ONCOTREE_SECONDARY_NODE": "AMPCA",
    },
    "TESTIS": {
        "CANCER_TYPE": "Testicular Cancer, NOS",
        "CANCER_TYPE_DETAILED": "Testis",
        "ONCOTREE_PRIMARY_NODE": "TESTIS",
        "ONCOTREE_SECONDARY_NODE": "",
    },
    "UCEC": {
        "CANCER_TYPE": "Endometrial Cancer",
        "CANCER_TYPE_DETAILED": "Endometrial Carcinoma",
        "ONCOTREE_PRIMARY_NODE": "UTERUS",
        "ONCOTREE_SECONDARY_NODE": "UCEC",
    },
}


@pytest.fixture
def clin_class(syn, genie_config):
    syn.tableQuery.side_effect = table_query_results
    return Clinical(syn, "SAGE", genie_config=genie_config)


@pytest.fixture(params=[(["foo"]), (["foo", "data_clinical_supp_sample_SAGE.txt"])])
def filename_fileformat_map(request):
    return request.param


@pytest.fixture
def valid_clinical_df():
    patientdf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            SEX=[1, 2, 1, 2, 99],
            PRIMARY_RACE=[1, 2, 3, 4, 99],
            SECONDARY_RACE=[1, 2, 3, 4, 99],
            TERTIARY_RACE=[1, 2, 3, 4, 99],
            ETHNICITY=[1, 2, 3, 4, 99],
            BIRTH_YEAR=[1222, "Unknown", 1920, 1990, 1990],
            CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
            YEAR_CONTACT=["Unknown", "Not Collected", ">89", "<18", 1990],
            INT_CONTACT=["Unknown", "Not Collected", ">32485", "<6570", 2000],
            YEAR_DEATH=["Unknown", "Not Collected", "Unknown", "Not Applicable", "<18"],
            INT_DOD=["Unknown", "Not Collected", "Unknown", "Not Applicable", "<6570"],
            DEAD=["Unknown", "Not Collected", "Unknown", False, True],
        )
    )

    sampledf = pd.DataFrame(
        dict(
            SAMPLE_ID=[
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            AGE_AT_SEQ_REPORT=[">32485", "Unknown", "<6570", 20000, 100000],
            ONCOTREE_CODE=["AMPCA", "AMPCA", "Unknown", "AMPCA", "AMPCA"],
            SAMPLE_TYPE=[1, 2, 3, 4, 4],
            SEQ_ASSAY_ID=["SAGE-1-1", "SAGE-SAGE-1", "SAGE-1", "SAGE-1", "SAGE-1"],
            SEQ_DATE=["Jan-2013", "ApR-2013", "Jul-2013", "Oct-2013", "release"],
            SAMPLE_CLASS=["Tumor", "cfDNA", "cfDNA", "cfDNA", "cfDNA"],
        )
    )

    clinicaldf = patientdf.merge(sampledf, on="PATIENT_ID")
    yield clinicaldf


def test_filetype(clin_class):
    assert clin_class._fileType == "clinical"


def test_incorrect_validatefilename(clin_class, filename_fileformat_map):
    filepath_list = filename_fileformat_map
    with pytest.raises(AssertionError):
        clin_class.validateFilename(filepath_list)


def test_correct_validatefilename(clin_class):
    assert clin_class.validateFilename(["data_clinical_supp_SAGE.txt"]) == "clinical"
    assert (
        clin_class.validateFilename(
            [
                "data_clinical_supp_sample_SAGE.txt",
                "data_clinical_supp_patient_SAGE.txt",
            ]
        )
        == "clinical"
    )


def test_patient_fillvs__process(clin_class):
    """
    This will be removed once vital status values are required
    - capitalized column headers
    - remapping of values
    - Fill out CENTER column
    - Append GENIE-CENTER-..
    """
    expected_patientdf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            SEX=["Male", "Female", "Male", "Female", "Unknown"],
            PRIMARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            SECONDARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            TERTIARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            ETHNICITY=["Test", "Why", "foo", "Me", "Unknown"],
            BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
        )
    )

    patientdf = pd.DataFrame(
        dict(
            PATIENT_Id=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            sex=[1, 2, 1, 2, 99],
            PRIMARY_RACE=[1, 2, 3, 4, 99],
            Secondary_RACE=[1, 2, 3, 4, 99],
            TERTIARY_RACE=[1, 2, 3, 4, 99],
            ETHNICITY=[1, 2, 3, 4, 99],
            BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
            CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
        )
    )
    patient_cols = [
        "PATIENT_ID",
        "SEX",
        "PRIMARY_RACE",
        "SECONDARY_RACE",
        "TERTIARY_RACE",
        "ETHNICITY",
        "BIRTH_YEAR",
        "CENTER",
        "YEAR_CONTACT",
        "YEAR_DEATH",
        "INT_CONTACT",
        "INT_DOD",
        "DEAD",
    ]
    clinical_template = pd.DataFrame(columns=patient_cols)
    new_patientdf = clin_class._process(patientdf, clinical_template)
    assert new_patientdf.columns.isin(expected_patientdf.columns).all()
    assert expected_patientdf.equals(new_patientdf[expected_patientdf.columns])


def test_patient_lesscoltemplate__process(clin_class):
    """
    Test scope is excluding values.
    Only those value defined by the scope will be written out
    """
    expected_patientdf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            SEX=["Male", "Female", "Male", "Female", "Unknown"],
            PRIMARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            SECONDARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            TERTIARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            ETHNICITY=["Test", "Why", "foo", "Me", "Unknown"],
            BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
        )
    )
    # TEST patient processing
    patientdf = pd.DataFrame(
        dict(
            PATIENT_Id=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            sex=[1, 2, 1, 2, 99],
            PRIMARY_RACE=[1, 2, 3, 4, 99],
            Secondary_RACE=[1, 2, 3, 4, 99],
            TERTIARY_RACE=[1, 2, 3, 4, 99],
            ETHNICITY=[1, 2, 3, 4, 99],
            BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
            CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
        )
    )
    patient_cols = [
        "PATIENT_ID",
        "SEX",
        "PRIMARY_RACE",
        "SECONDARY_RACE",
        "TERTIARY_RACE",
        "ETHNICITY",
        "BIRTH_YEAR",
        "CENTER",
        "YEAR_CONTACT",
        "YEAR_DEATH",
    ]
    clinical_template = pd.DataFrame(columns=patient_cols)
    new_patientdf = clin_class._process(patientdf, clinical_template)

    assert new_patientdf.columns.isin(patient_cols).all()
    assert expected_patientdf[expected_patientdf.columns].equals(
        new_patientdf[expected_patientdf.columns]
    )


def test_patient_vs__process(clin_class):
    """
    Test vital status columns being propogated with same data
    """
    expected_patientdf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            SEX=["Male", "Female", "Male", "Female", "Unknown"],
            PRIMARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            SECONDARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            TERTIARY_RACE=["Test", "Why", "foo", "Me", "Unknown"],
            ETHNICITY=["Test", "Why", "foo", "Me", "Unknown"],
            BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
            YEAR_DEATH=["Unknown", "Not Collected", "Not Applicable", 1990, 1990],
            YEAR_CONTACT=["Unknown", "Not Collected", 1990, 1990, 1990],
            INT_CONTACT=[">32485", "<6570", "Unknown", "Not Collected", 2000],
            INT_DOD=[">32485", "<6570", "Unknown", "Not Collected", "Not Applicable"],
            DEAD=[True, False, "Unknown", "Not Collected", True],
        )
    )
    # TEST patient processing
    # Clinical file headers are capitalized prior to processing
    patientdf = pd.DataFrame(
        dict(
            PATIENT_Id=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            sex=[1, 2, 1, 2, 99],
            PRIMARY_RACE=[1, 2, 3, 4, 99],
            Secondary_RACE=[1, 2, 3, 4, 99],
            TERTIARY_RACE=[1, 2, 3, 4, 99],
            ETHNICITY=[1, 2, 3, 4, 99],
            BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
            CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
            YEAR_DEATH=["Unknown", "Not Collected", "Not Applicable", 1990, 1990],
            YEAR_CONTACT=["Unknown", "Not Collected", 1990, 1990, 1990],
            INT_CONTACT=[">32485", "<6570", "Unknown", "Not Collected", 2000],
            INT_DOD=[">32485", "<6570", "Unknown", "Not Collected", "Not Applicable"],
            DEAD=[True, False, "Unknown", "Not Collected", True],
        )
    )
    patient_cols = [
        "PATIENT_ID",
        "SEX",
        "PRIMARY_RACE",
        "SECONDARY_RACE",
        "TERTIARY_RACE",
        "ETHNICITY",
        "BIRTH_YEAR",
        "CENTER",
        "YEAR_CONTACT",
        "YEAR_DEATH",
        "INT_CONTACT",
        "INT_DOD",
        "DEAD",
    ]
    clinical_template = pd.DataFrame(columns=patient_cols)
    new_patientdf = clin_class._process(patientdf, clinical_template)
    assert expected_patientdf.equals(new_patientdf[expected_patientdf.columns])


def test_sample__process(clin_class):
    """
    Test sample processing
    - column headers are capitalized
    - SEQ_DATE is normalized (Mon-YYYY, Release)
    - Allow UNKNOWN oncotree value
    - Add on GENIE-CENTER-...
    - Remapping of SAMPLE_TYPE/SAMPLE_TYPE_DETAILED value
    - SEQ_YEAR from SEQ_DATE, nan if SEQ_DATE is Release

    """
    expected_sampledf = pd.DataFrame(
        dict(
            SAMPLE_ID=[
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            AGE_AT_SEQ_REPORT=[100000, 100000, 100000, 100000, 100000],
            ONCOTREE_CODE=["AMPCA", "UNKNOWN", "AMPCA", "AMPCA", "AMPCA"],
            SAMPLE_TYPE=["Test", "Why", "foo", "Me", "Me"],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
            SAMPLE_TYPE_DETAILED=["non", "asdf", "asdf", "asdff", "asdff"],
            SEQ_ASSAY_ID=["SAGE-1", "SAGE-1", "SAGE-1", "SAGE-1", "SAGE-1"],
            SEQ_DATE=["Jan-2012", "Apr-2013", "Jul-2014", "Oct-2015", "Release"],
            SEQ_YEAR=[2012, 2013, 2014, 2015, float("nan")],
        )
    )
    sample_cols = [
        "SAMPLE_ID",
        "PATIENT_ID",
        "AGE_AT_SEQ_REPORT",
        "ONCOTREE_CODE",
        "SAMPLE_TYPE",
        "SEQ_ASSAY_ID",
        "SEQ_DATE",
        "SAMPLE_TYPE_DETAILED",
        "SEQ_YEAR",
    ]

    clinical_template = pd.DataFrame(columns=sample_cols)
    # patient = False
    sampledf = pd.DataFrame(
        dict(
            SAMPLE_ID=[
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            Age_AT_SEQ_REPORT=[100000, 100000, 100000, 100000, 100000],
            ONCOTree_CODE=["AMPCA", " UNKNOWN", "AMPCA", "AMPCA", "AMPCA"],
            SAMPLE_TYPE=[1, 2, 3, 4, 4],
            SEQ_ASSAY_ID=["SAGE-1", "SAGE-1", "SAGE-1", "SAGE-1", "SAGE-1"],
            SEQ_DATE=["Jan-2012", "Apr-2013", "JUL-2014", "Oct-2015", "release"],
        )
    )

    new_sampledf = clin_class._process(sampledf, clinical_template)
    assert new_sampledf.columns.isin(expected_sampledf.columns).all()
    assert expected_sampledf.equals(new_sampledf[expected_sampledf.columns])


def test_perfect__validate(clin_class, valid_clinical_df):
    """
    Test perfect validation
    """
    with mock.patch(
        "genie.process_functions.get_oncotree_code_mappings", return_value=onco_map_dict
    ) as mock_get_onco_map:
        error, warning = clin_class._validate(valid_clinical_df)
        mock_get_onco_map.called_once_with(json_oncotreeurl)
        assert error == ""
        assert warning == ""


def test_nonull__validate(clin_class):
    """
    Test that no null values are allowed in the clinical dataframe
    """
    patientdf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            SEX=[1, 2, 1, 2, float("nan")],
            PRIMARY_RACE=[1, 2, 3, 4, float("nan")],
            SECONDARY_RACE=[1, 2, 3, 4, float("nan")],
            TERTIARY_RACE=[1, 2, 3, 4, float("nan")],
            ETHNICITY=[1, 2, 3, 4, float("nan")],
            BIRTH_YEAR=[float("nan"), "Unknown", 1920, 1990, 1990],
            CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
            YEAR_DEATH=[
                "Unknown",
                "Not Collected",
                "Not Applicable",
                1990,
                float("nan"),
            ],
            YEAR_CONTACT=["Unknown", "Not Collected", float("nan"), 1990, 1990],
            INT_CONTACT=["Unknown", "Not Collected", ">32485", float("nan"), 2000],
            INT_DOD=["Unknown", "Not Collected", "Unknown", float("nan"), "<6570"],
            DEAD=["Unknown", "Not Collected", "Unknown", float("nan"), True],
        )
    )

    sampledf = pd.DataFrame(
        dict(
            SAMPLE_ID=[
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            AGE_AT_SEQ_REPORT=[100000, "Unknown", 20000, float("nan"), 100000],
            ONCOTREE_CODE=["AMPCA", "AMPCA", "Unknown", "AMPCA", "AMPCA"],
            SAMPLE_TYPE=[1, 2, 3, 4, float("nan")],
            SEQ_ASSAY_ID=["SAGE-1-1", "SAGE-SAGE-1", "SAGE-1", "SAGE-1", "SAGE-1"],
            SEQ_DATE=["Jan-2013", "ApR-2013", "Jul-2013", "Oct-2013", "release"],
            SAMPLE_CLASS=["Tumor", "cfDNA", "cfDNA", "cfDNA", float("nan")],
        )
    )

    clinicaldf = patientdf.merge(sampledf, on="PATIENT_ID")
    with mock.patch(
        "genie.process_functions.get_oncotree_code_mappings", return_value=onco_map_dict
    ) as mock_get_onco_map:
        error, warning = clin_class._validate(clinicaldf)
        mock_get_onco_map.called_once_with(json_oncotreeurl)
        expected_errors = (
            "Sample Clinical File: Please double check your "
            "AGE_AT_SEQ_REPORT. It must be an integer, 'Unknown', "
            "'>32485', '<6570'.\n"
            "Sample Clinical File: Please double check your SAMPLE_TYPE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Patient Clinical File: Please double check your BIRTH_YEAR "
            "column, it must be an integer in YYYY format <= {year} or "
            "'Unknown', '>89', '<18'.\n"
            "Patient Clinical File: Please double check your YEAR_DEATH "
            "column, it must be an integer in YYYY format <= {year} or "
            "'Unknown', 'Not Collected', 'Not Applicable', "
            "'Not Released', '>89', '<18'.\n"
            "Patient Clinical File: Please double check your YEAR_CONTACT "
            "column, it must be an integer in YYYY format <= {year} or "
            "'Unknown', 'Not Collected', 'Not Released', '>89', '<18'.\n"
            "Patient Clinical File: Please double check your INT_CONTACT "
            "column, it must be an integer, '>32485', '<6570', 'Unknown', "
            "'Not Released' or 'Not Collected'.\n"
            "Patient Clinical File: Please double check your INT_DOD "
            "column, it must be an integer, '>32485', '<6570', 'Unknown', "
            "'Not Collected', 'Not Released' or 'Not Applicable'.\n"
            "Patient Clinical File: Please double check your DEAD column, "
            "it must be True, False, 'Unknown', "
            "'Not Released' or 'Not Collected'.\n"
            "Patient: you have inconsistent redaction values in "
            "YEAR_CONTACT, INT_CONTACT.\n"
            "Patient: you have inconsistent redaction and text values in "
            "YEAR_DEATH, INT_DOD.\n"
            "Sample Clinical File: SAMPLE_CLASS column must be 'Tumor', or 'cfDNA'\n"
            "Patient Clinical File: Please double check your PRIMARY_RACE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Patient Clinical File: Please double check your SECONDARY_RACE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Patient Clinical File: Please double check your TERTIARY_RACE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Patient Clinical File: Please double check your SEX column.  "
            "This column must only be these values: 1, 2, 99\n"
            "Patient Clinical File: Please double check your ETHNICITY "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
        )
        assert error == expected_errors.format(year=datetime.datetime.utcnow().year)
        assert warning == ""


def test_missingcols__validate(clin_class):
    """
    Test for missing column errors
    """
    clinicaldf = pd.DataFrame()
    with mock.patch(
        "genie.process_functions.get_oncotree_code_mappings", return_value=onco_map_dict
    ) as mock_get_onco_map:
        error, warning = clin_class._validate(clinicaldf)
        mock_get_onco_map.called_once_with(json_oncotreeurl)
        expected_errors = (
            "Sample Clinical File: Must have SAMPLE_ID column.\n"
            "Patient Clinical File: Must have PATIENT_ID column.\n"
            "Sample Clinical File: Must have AGE_AT_SEQ_REPORT column.\n"
            "Sample Clinical File: Must have ONCOTREE_CODE column.\n"
            "Sample Clinical File: Must have SAMPLE_TYPE column.\n"
            "Sample Clinical File: Must have SEQ_ASSAY_ID column.\n"
            "Sample Clinical File: Must have SEQ_DATE column.\n"
            "Patient Clinical File: Must have BIRTH_YEAR column.\n"
            "Patient Clinical File: Must have YEAR_DEATH column.\n"
            "Patient Clinical File: Must have YEAR_CONTACT column.\n"
            "Patient Clinical File: Must have INT_CONTACT column.\n"
            "Patient Clinical File: Must have INT_DOD column.\n"
            "Patient Clinical File: Must have DEAD column.\n"
            "Patient Clinical File: Must have SEX column.\n"
        )

        expected_warnings = (
            "Patient Clinical File: Doesn't have PRIMARY_RACE column. "
            "This column will be added\n"
            "Patient Clinical File: Doesn't have SECONDARY_RACE column. "
            "This column will be added\n"
            "Patient Clinical File: Doesn't have TERTIARY_RACE column. "
            "This column will be added\n"
            "Patient Clinical File: Doesn't have ETHNICITY column. "
            "This column will be added\n"
        )
        assert error == expected_errors
        assert warning == expected_warnings


def test_errors__validate(clin_class):
    """
    Test for validation errors
    """
    sampleDf = pd.DataFrame(
        dict(
            SAMPLE_ID=[
                float("nan"),
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            PATIENT_ID=[
                "GENIE-SAGE-ID6",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                float("nan"),
                "GENIE-SAGE-ID5",
            ],
            AGE_AT_SEQ_REPORT=[10, 100000, ">doo", 100000, 100000],
            ONCOTREE_CODE=["AMPCAD", "TESTIS", "AMPCA", "AMPCA", "UCEC"],
            SAMPLE_TYPE=[1, 2, 3, 4, 6],
            SEQ_ASSAY_ID=[float("nan"), "Sage-1", "SAGE-1", "S-SAGE-1", "SAGE-1"],
            SEQ_DATE=["Jane-2013", "Jan-2013", "Jan-2013", "Jan-2013", "Jan-2013"],
            YEAR_DEATH=["Unknown", "Not Collected", "Not Applicable", 19930, 1990],
            YEAR_CONTACT=["Unknown", "Not Collected", 1990, 1990, 19940],
            INT_CONTACT=[">32485", "<6570", 1990, "Not Collected", ">foobar"],
            INT_DOD=[">32485", "<6570", "Unknown", "Not Collected", "<dense"],
            DEAD=[1, False, "Unknown", "Not Collected", "Not Applicable"],
        )
    )

    patientDf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID6",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                float("nan"),
                "GENIE-SAGE-ID5",
            ],
            SEX=[1, 2, 1, 5, float("nan")],
            PRIMARY_RACE=[1, 2, 3, 6, float("nan")],
            SECONDARY_RACE=[1, 2, 3, 6, float("nan")],
            TERTIARY_RACE=[1, 2, 3, 6, float("nan")],
            ETHNICITY=[1, 2, 3, 6, float("nan")],
            BIRTH_YEAR=[1990, 1990, ">90", 1990, 1990],
            CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
        )
    )
    clinicalDf = patientDf.merge(sampleDf, on="PATIENT_ID")
    with mock.patch(
        "genie.process_functions.get_oncotree_code_mappings", return_value=onco_map_dict
    ) as mock_get_onco_map:
        error, warning = clin_class._validate(clinicalDf)
        mock_get_onco_map.called_once_with(json_oncotreeurl)

        expectedErrors = (
            "Sample Clinical File: SAMPLE_ID must start with GENIE-SAGE\n"
            "Patient Clinical File: PATIENT_ID must start with GENIE-SAGE\n"
            "Sample Clinical File: PATIENT_ID's much be contained in the "
            "SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n"
            "Patient Clinical File: All samples must have associated patient "
            "information and no null patient ids allowed. "
            "These samples are missing patient data: GENIE-SAGE-ID4-1\n"
            "Sample Clinical File: Please double check your "
            "AGE_AT_SEQ_REPORT. It must be an integer, 'Unknown', "
            "'>32485', '<6570'.\n"
            "Sample Clinical File: Please double check that all your "
            "ONCOTREE CODES exist in the mapping. You have 1 samples that "
            "don't map. These are the codes that don't map: AMPCAD\n"
            "Sample Clinical File: Please double check your SAMPLE_TYPE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Sample Clinical File: Please double check your SEQ_ASSAY_ID "
            "columns, there are empty rows.\n"
            "Sample Clinical File: Please make sure your SEQ_ASSAY_IDs start "
            "with your center abbreviation: S-SAGE-1.\n"
            "Sample Clinical File: SEQ_DATE must be one of five values- "
            "For Jan-March: use Jan-YEAR. "
            "For Apr-June: use Apr-YEAR. "
            "For July-Sep: use Jul-YEAR. "
            "For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) "
            "For values that don't have SEQ_DATES that you want "
            "released use 'release'.\n"
            "Patient Clinical File: Please double check your BIRTH_YEAR "
            "column, it must be an integer in YYYY format <= {year} or "
            "'Unknown', '>89', '<18'.\n"
            "Patient Clinical File: Please double check your YEAR_DEATH "
            "column, it must be an integer in YYYY format <= {year} or "
            "'Unknown', 'Not Collected', 'Not Applicable', "
            "'Not Released', '>89', '<18'.\n"
            "Patient Clinical File: Please double check your YEAR_CONTACT "
            "column, it must be an integer in YYYY format <= {year} or "
            "'Unknown', 'Not Collected', 'Not Released', '>89', '<18'.\n"
            "Patient Clinical File: Please double check your INT_CONTACT "
            "column, it must be an integer, '>32485', '<6570', 'Unknown', "
            "'Not Released' or 'Not Collected'.\n"
            "Patient Clinical File: Please double check your INT_DOD column, "
            "it must be an integer, '>32485', '<6570', 'Unknown', "
            "'Not Collected', 'Not Released' or 'Not Applicable'.\n"
            "Patient Clinical File: Please double check your DEAD column, "
            "it must be True, False, 'Unknown', "
            "'Not Released' or 'Not Collected'.\n"
            "Patient: you have inconsistent redaction and text values in "
            "YEAR_CONTACT, INT_CONTACT.\n"
            "Patient: you have inconsistent redaction and text values in "
            "YEAR_DEATH, INT_DOD.\n"
            "Patient Clinical File: DEAD value is inconsistent with "
            "INT_DOD for at least one patient.\n"
            "Patient Clinical File: Please double check your PRIMARY_RACE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Patient Clinical File: Please double check your SECONDARY_RACE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Patient Clinical File: Please double check your TERTIARY_RACE "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
            "Patient Clinical File: Please double check your SEX column.  "
            "This column must only be these values: 1, 2, 99\n"
            "Patient Clinical File: Please double check your ETHNICITY "
            "column.  This column must only be these values: 1, 2, 3, 4, 99\n"
        )
        expectedWarnings = (
            "Sample Clinical File: All patients must have associated sample "
            "information. These patients are missing sample data: GENIE-SAGE-ID6\n"
            "Sample Clinical File: Some SAMPLE_IDs have conflicting SEX and "
            "ONCOTREE_CODES: GENIE-SAGE-ID2-1,GENIE-SAGE-ID5-1\n"
        )
        assert error == expectedErrors.format(year=datetime.datetime.utcnow().year)
        assert warning == expectedWarnings


def test_duplicated__validate(clin_class):
    """
    Test for duplicated SAMPLE_ID error and and in the case that
    both sample and patient
    are uploaded, it could be a duplicated PATIENT_ID error
    """
    patientDf = pd.DataFrame(
        dict(
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                float("nan"),
            ],
            SEX=[1, 2, 1, 2, float("nan")],
            PRIMARY_RACE=[1, 2, 3, 4, float("nan")],
            SECONDARY_RACE=[1, 2, 3, 4, float("nan")],
            TERTIARY_RACE=[1, 2, 3, 4, float("nan")],
            ETHNICITY=[1, 2, 3, 4, float("nan")],
            BIRTH_YEAR=["Unknown", 1990, 1990, 1990, float("nan")],
            CENTER=["FOO", "FOO", "FOO", "FOO", float("nan")],
            YEAR_CONTACT=["Unknown", "Not Collected", ">89", "<18", float("nan")],
            INT_CONTACT=["Unknown", "Not Collected", ">32485", "<6570", float("nan")],
            YEAR_DEATH=[
                "Unknown",
                "Not Collected",
                "Unknown",
                "Not Applicable",
                float("nan"),
            ],
            INT_DOD=[
                "Unknown",
                "Not Collected",
                "Unknown",
                "Not Applicable",
                float("nan"),
            ],
            DEAD=["Unknown", "Not Collected", "Unknown", False, float("nan")],
        )
    )

    sampleDf = pd.DataFrame(
        dict(
            SAMPLE_ID=[
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                float("nan"),
            ],
            PATIENT_ID=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                float("nan"),
            ],
            AGE_AT_SEQ_REPORT=[100000, 100000, 100000, float("nan")],
            ONCOTREE_CODE=["AMPCA", "UNKNOWN", "AMPCA", float("nan")],
            SAMPLE_TYPE=[1, 3, 4, float("nan")],
            SEQ_ASSAY_ID=["SAGE-1-1", "SAGE-1", "SAGE-1", float("nan")],
            SEQ_DATE=["Jan-2013", "Jul-2013", "Oct-2013", float("nan")],
        )
    )

    clinicalDf = patientDf.merge(sampleDf, on="PATIENT_ID")
    with mock.patch(
        "genie.process_functions.get_oncotree_code_mappings", return_value=onco_map_dict
    ) as mock_get_onco_map:
        error, warning = clin_class._validate(clinicalDf)
        mock_get_onco_map.called_once_with(json_oncotreeurl)
        expectedErrors = (
            "Clinical file(s): No empty rows allowed.\n"
            "Sample Clinical File: No duplicated SAMPLE_ID allowed.\n"
            "If there are no duplicated SAMPLE_IDs, and both sample and "
            "patient files are uploaded, then please check to make sure no "
            "duplicated PATIENT_IDs exist in the patient clinical file.\n"
        )

        assert error == expectedErrors
        assert warning == ""


class fake_oncotree:
    text = json.dumps(
        {
            "TISSUE": {
                "children": {
                    "AMPCA": {
                        "level": 1,
                        "mainType": "Ampullary Cancer",
                        "name": "Ampullary Carcinoma",
                        "children": {
                            "TESTIS": {
                                "level": 2,
                                "mainType": "Testicular Cancer, NOS",
                                "name": "Testis",
                                "children": [],
                            },
                            "UCEC": {
                                "level": 2,
                                "mainType": "Endometrial Cancer",
                                "name": "Endometrial Carcinoma",
                                "children": [],
                            },
                        },
                    }
                }
            }
        }
    )


expected_onco_mapping = {
    "AMPCA": {
        "CANCER_TYPE": "Ampullary Cancer",
        "CANCER_TYPE_DETAILED": "Ampullary Carcinoma",
        "ONCOTREE_PRIMARY_NODE": "AMPCA",
        "ONCOTREE_SECONDARY_NODE": "",
    },
    "TESTIS": {
        "CANCER_TYPE": "Testicular Cancer, NOS",
        "CANCER_TYPE_DETAILED": "Testis",
        "ONCOTREE_PRIMARY_NODE": "AMPCA",
        "ONCOTREE_SECONDARY_NODE": "TESTIS",
    },
    "UCEC": {
        "CANCER_TYPE": "Endometrial Cancer",
        "CANCER_TYPE_DETAILED": "Endometrial Carcinoma",
        "ONCOTREE_PRIMARY_NODE": "AMPCA",
        "ONCOTREE_SECONDARY_NODE": "UCEC",
    },
}


def test_get_oncotree_code_mappings():
    from genie import process_functions

    with mock.patch(
        "genie.process_functions.retry_get_url", return_value=fake_oncotree
    ) as retry_get_url:
        onco_mapping = process_functions.get_oncotree_code_mappings(json_oncotreeurl)
        retry_get_url.called_once_with(json_oncotreeurl)
        assert onco_mapping == expected_onco_mapping


def test__check_year_no_errors():
    """Tests perfect checking year validation function"""
    perfectdf = pd.DataFrame({"BIRTH_YEAR": ["Unknown", 1990, 1990, 1990, 1990]})
    error = genie_registry.clinical._check_year(
        clinicaldf=perfectdf,
        year_col="BIRTH_YEAR",
        filename="Filename",
        allowed_string_values=["Unknown"],
    )
    assert error == ""


def test__check_year_too_big_year():
    """Tests year can't be greater than current year"""
    year_now = datetime.datetime.utcnow().year

    errordf = pd.DataFrame({"BIRTH_YEAR": ["Unknown", year_now + 1, 1990, 1990, 1990]})
    error = genie_registry.clinical._check_year(
        errordf, "BIRTH_YEAR", "Filename", allowed_string_values=["Unknown"]
    )
    assert error == (
        "Filename: Please double check your BIRTH_YEAR column, "
        f"it must be an integer in YYYY format <= {year_now} or 'Unknown'.\n"
    )


def test__check_year_invalid():
    """Tests year can't be greater than current year"""
    year_now = datetime.datetime.utcnow().year

    errordf = pd.DataFrame({"BIRTH_YEAR": ["Unknown", 1990, 1990, 1990, 1990]})
    error = genie_registry.clinical._check_year(errordf, "BIRTH_YEAR", "Filename")
    assert error == (
        "Filename: Please double check your BIRTH_YEAR column, "
        f"it must be an integer in YYYY format <= {year_now}.\n"
    )


def test_remap_clinical_values_sampletype():
    """Test adding of sample type detailed when remapping clinical values"""
    testdf = pd.DataFrame({"SAMPLE_TYPE": [1, 2, 99]})
    expecteddf = pd.DataFrame(
        {
            "SAMPLE_TYPE": ["Test", "Why", "Unknown"],
            "SAMPLE_TYPE_DETAILED": ["non", "asdf", "asdfasdf"],
        }
    )
    remappeddf = genie_registry.clinical.remap_clinical_values(
        testdf, sexdf, no_nan, sexdf, no_nan
    )
    assert expecteddf.equals(remappeddf)


@pytest.mark.parametrize(
    "col", ["SEX", "PRIMARY_RACE", "SECONDARY_RACE", "TERTIARY_RACE", "ETHNICITY"]
)
def test_remap_clinical_values(col):
    """Test Remapping clinical values"""
    testdf = pd.DataFrame({col: [1, 2, 99]})
    expecteddf = pd.DataFrame({col: ["Male", "Female", "Unknown"]})
    remappeddf = genie_registry.clinical.remap_clinical_values(
        testdf, sexdf, sexdf, sexdf, sexdf
    )
    assert expecteddf.equals(remappeddf)


def test__check_int_year_consistency_valid():
    """Test valid vital status consistency"""
    testdf = pd.DataFrame(
        {
            "INT_2": [1, 2, "Unknown", "Unknown"],
            "YEAR_1": [1, 4, "Unknown", 1],
            "FOO_3": [1, 3, "Unknown", 1],
        }
    )
    error = genie_registry.clinical._check_int_year_consistency(
        clinicaldf=testdf, cols=["INT_2", "YEAR_1"], string_vals=["Unknown"]
    )
    assert error == ""


@pytest.mark.parametrize(
    ("inconsistent_df", "expected_err"),
    [
        (
            pd.DataFrame({"INT_2": [1, ">32485", 4], "YEAR_1": [1, 4, 4]}),
            "Patient: you have inconsistent redaction values in INT_2, YEAR_1.\n",
        ),
        (
            pd.DataFrame(
                {"INT_2": [1, 3, "Unknown"], "YEAR_1": [1, 4, "Not Applicable"]}
            ),
            "Patient: you have inconsistent text values in INT_2, YEAR_1.\n",
        ),
        (
            pd.DataFrame(
                {
                    "INT_2": [1, "Unknown", "Unknown"],
                    "YEAR_1": [1, "Not Applicable", "Unknown"],
                }
            ),
            "Patient: you have inconsistent text values in INT_2, YEAR_1.\n",
        ),
        (
            pd.DataFrame({"INT_2": [1, 2, ">32485"], "YEAR_1": [1, 4, "<18"]}),
            "Patient: you have inconsistent redaction values in INT_2, YEAR_1.\n",
        ),
        (
            pd.DataFrame({"INT_2": [1, 2, "<6570"], "YEAR_1": [1, 4, ">89"]}),
            "Patient: you have inconsistent redaction values in INT_2, YEAR_1.\n",
        ),
        (
            pd.DataFrame(
                {
                    "INT_2": ["<6570", "Unknown", "Unknown"],
                    "YEAR_1": [1, "Not Applicable", "Unknown"],
                }
            ),
            "Patient: you have inconsistent redaction and text values in "
            "INT_2, YEAR_1.\n",
        ),
        (
            pd.DataFrame(
                {
                    "INT_2": ["12345", "Unknown", "Unknown"],
                    "YEAR_1": ["Unknown", "Unknown", "Unknown"],
                }
            ),
            "Patient: you have inconsistent text values in INT_2, YEAR_1.\n",
        ),
    ],
)
def test__check_int_year_consistency_inconsistent(inconsistent_df, expected_err):
    """Test inconsistent vital status values"""
    error = genie_registry.clinical._check_int_year_consistency(
        clinicaldf=inconsistent_df,
        cols=["INT_2", "YEAR_1"],
        string_vals=["Unknown", "Not Applicable"],
    )
    assert error == expected_err


@pytest.mark.parametrize(
    "valid_df",
    [
        pd.DataFrame({"INT_DOD": [11111, "Not Applicable"], "DEAD": [True, False]}),
        pd.DataFrame(
            {"INT_DOD": ["Not Applicable", "Not Applicable"], "DEAD": [False, False]}
        ),
    ],
)
def test__check_int_dead_consistency_valid(valid_df):
    """Test valid vital status consistency"""
    error = genie_registry.clinical._check_int_dead_consistency(clinicaldf=valid_df)
    assert error == ""


@pytest.mark.parametrize(
    "inconsistent_df",
    [
        pd.DataFrame(
            {"INT_DOD": ["Not Applicable", "Not Applicable"], "DEAD": [True, False]}
        ),
        pd.DataFrame({"INT_DOD": [1111, "Not Released"], "DEAD": [True, False]}),
        pd.DataFrame({"INT_DOD": [1111, 11111], "DEAD": [True, False]}),
        pd.DataFrame(
            {"INT_DOD": [1111, "Not Released"], "DEAD": [True, "Not Applicable"]}
        ),
    ],
)
def test__check_int_dead_consistency_inconsistent(inconsistent_df):
    """Test valid vital status consistency"""
    error = genie_registry.clinical._check_int_dead_consistency(
        clinicaldf=inconsistent_df
    )
    assert error == (
        "Patient Clinical File: DEAD value is inconsistent with "
        "INT_DOD for at least one patient.\n"
    )


def get_cross_validate_bed_files_test_cases():
    return [
        {
            "name": "all_match",
            "test_clinical_df": pd.DataFrame(
                {
                    "SEQ_ASSAY_ID": [
                        "SAGE-1-1",
                        "SAGE-SAGE-1",
                        "SAGE-1",
                        "SAGE-1",
                        "SaGe-1",
                    ]
                }
            ),
            "test_ancillary_files": [
                [{"name": "SAGE-SAGE-1.bed", "path": ""}],
                [{"name": "SAGE-1-1.bed", "path": ""}],
                [{"name": "SAGE-1.bed", "path": ""}],
            ],
            "expected_missing_files": [],
        },
        {
            "name": "partial_match",
            "test_clinical_df": pd.DataFrame(
                {"SEQ_ASSAY_ID": ["SAGE-1-1", "SAGE-1-2", "SaGe-1_1"]}
            ),
            "test_ancillary_files": [
                [{"name": "SAGE-SAGE-1.bed", "path": ""}],
                [{"name": "SAGE-1-1.bed", "path": ""}],
                [{"name": "SAGE-1.bed", "path": ""}],
            ],
            "expected_missing_files": ["SAGE-1-2.bed"],
        },
        {
            "name": "no_match",
            "test_clinical_df": pd.DataFrame(
                {"SEQ_ASSAY_ID": ["SAGE-1-2", "SAGE-1-3", "SaGe_1_2"]}
            ),
            "test_ancillary_files": [
                [{"name": "SAGE-SAGE-1.bed", "path": ""}],
                [{"name": "SAGE-1-1.bed", "path": ""}],
                [{"name": "SAGE-1.bed", "path": ""}],
            ],
            "expected_missing_files": ["SAGE-1-2.bed", "SAGE-1-3.bed"],
        },
        {
            "name": "no_bed_files",
            "test_clinical_df": pd.DataFrame(
                {"SEQ_ASSAY_ID": ["SAGE-1-2", "SAGE-1-3", "SAge-1_2"]}
            ),
            "test_ancillary_files": [
                [{"name": "SAGE-1.txt", "path": ""}],
            ],
            "expected_missing_files": ["SAGE-1-2.bed", "SAGE-1-3.bed"],
        },
    ]


@pytest.mark.parametrize(
    "test_cases", get_cross_validate_bed_files_test_cases(), ids=lambda x: x["name"]
)
def test_that_cross_validate_bed_files_exist_returns_correct_msgs(
    clin_class, test_cases
):
    clin_class.ancillary_files = test_cases["test_ancillary_files"]
    missing_files = clin_class._cross_validate_bed_files_exist(
        test_cases["test_clinical_df"]
    )
    assert Counter(test_cases["expected_missing_files"]) == Counter(missing_files)


def test_that_cross_validate_bed_files_exist_calls_expected_methods(clin_class):
    test_clinical_df = pd.DataFrame({"SEQ_ASSAY_ID": ["SAGE-SAGE-1"]})
    test_ancillary_files = (
        [
            [{"name": "SAGE-SAGE-1.bed", "path": ""}],
        ],
    )
    with mock.patch.object(
        validate,
        "parse_file_info_in_nested_list",
        return_value={"files": "SAGE-1-1.bed"},
    ) as patch_parse_file_info:
        clin_class.ancillary_files = test_ancillary_files
        clin_class._cross_validate_bed_files_exist(test_clinical_df)
        patch_parse_file_info.assert_called_once_with(
            nested_list=clin_class.ancillary_files,
            search_str="sage-sage-1.bed",
            ignore_case=True,
            allow_underscore=True,
        )


@pytest.mark.parametrize(
    "missing_files,expected_error,expected_warning",
    [
        (
            [],
            "",
            "",
        ),
        (
            ["test1.bed", "test2.bed"],
            "At least one SEQ_ASSAY_ID in your clinical file does not have an associated BED file. "
            "Please update your file(s) to be consistent.\n"
            "Missing BED files: test1.bed, test2.bed\n",
            "",
        ),
    ],
    ids=["no_missing_files", "missing_files"],
)
def test_that_cross_validate_bed_files_exist_message_returns_correct_msgs(
    clin_class, missing_files, expected_error, expected_warning
):
    errors, warnings = clin_class._cross_validate_bed_files_exist_message(missing_files)
    assert errors == expected_error
    assert warnings == expected_warning


def test_that__cross_validate_calls_expected_methods(clin_class):
    with mock.patch.object(
        Clinical, "_cross_validate_assay_info_has_seq", return_value=("", "")
    ) as patch__cross_validate_assay, mock.patch.object(
        Clinical, "_cross_validate_bed_files_exist", return_value=("", "")
    ) as patch__cross_validate_bed, mock.patch.object(
        Clinical, "_cross_validate_bed_files_exist_message", return_value=("", "")
    ) as patch__cross_validate_bed_msg:
        clin_class._cross_validate(clinicaldf=pd.DataFrame({"something": [1]}))
        patch__cross_validate_assay.assert_called_once()
        patch__cross_validate_bed.assert_called_once()
        patch__cross_validate_bed_msg.assert_called_once()


def test_that__cross_validate_returns_correct_format_for_errors_warnings(clin_class):
    with mock.patch.object(
        Clinical, "_cross_validate_assay_info_has_seq", return_value=("test1", "")
    ) as patch__cross_validate_assay, mock.patch.object(
        Clinical, "_cross_validate_bed_files_exist", return_value=["something_missing"]
    ) as patch__cross_validate_bed, mock.patch.object(
        Clinical,
        "_cross_validate_bed_files_exist_message",
        return_value=("test3\n", ""),
    ) as patch__cross_validate_bed_msg:
        errors, warnings = clin_class._cross_validate(
            clinicaldf=pd.DataFrame({"something": [1]})
        )
        assert errors == "test1test3\n"
        assert warnings == ""


def test_that__cross_validate_assay_info_has_seq_does_not_read_files_if_no_assay_files(
    clin_class,
):
    with mock.patch.object(
        validate,
        "parse_file_info_in_nested_list",
        return_value={"files": {}, "file_info": {"name": "", "path": ""}},
    ) as patch_assay_files, mock.patch.object(
        process_functions,
        "get_assay_dataframe",
    ) as patch_get_df:
        errors, warnings = clin_class._cross_validate_assay_info_has_seq(
            pd.DataFrame({})
        )
        assert warnings == ""
        assert errors == ""
        patch_get_df.assert_not_called()


def test_that__cross_validate_assay_info_has_seq_does_not_call_check_col_exist_if_assay_df_read_error(
    clin_class,
):
    with mock.patch.object(
        validate,
        "parse_file_info_in_nested_list",
        return_value={"files": {"some_file"}, "file_info": {"name": "", "path": ""}},
    ) as patch_assay_files, mock.patch.object(
        process_functions,
        "get_assay_dataframe",
        side_effect=Exception("mocked error"),
    ) as patch_get_df, mock.patch.object(
        process_functions, "checkColExist"
    ) as patch_check_col_exist:
        errors, warnings = clin_class._cross_validate_assay_info_has_seq(
            pd.DataFrame({})
        )
        assert warnings == ""
        assert errors == ""
        patch_check_col_exist.assert_not_called()


def test_that__cross_validate_assay_info_has_seq_does_not_call_check_values_if_id_cols_do_not_exist(
    clin_class,
):
    with mock.patch.object(
        validate,
        "parse_file_info_in_nested_list",
        return_value={"files": {"some_file"}, "file_info": {"name": "", "path": ""}},
    ) as patch_assay_files, mock.patch.object(
        process_functions,
        "get_assay_dataframe",
        return_value=pd.DataFrame({"test_col": [2, 3, 4]}),
    ) as patch_get_df, mock.patch.object(
        process_functions, "checkColExist", return_value=False
    ) as patch_check_col_exist, mock.patch.object(
        validate, "check_values_between_two_df"
    ) as patch_check_values:
        errors, warnings = clin_class._cross_validate_assay_info_has_seq(
            pd.DataFrame({})
        )
        assert warnings == ""
        assert errors == ""
        patch_check_values.assert_not_called()


def test_that__cross_validate_assay_info_has_seq_calls_check_values_between_two_df(
    clin_class,
):
    with mock.patch.object(
        validate,
        "parse_file_info_in_nested_list",
        return_value={"files": {"some_file"}, "file_info": {"name": "", "path": ""}},
    ) as patch_assay_files, mock.patch.object(
        process_functions,
        "get_assay_dataframe",
        return_value=None,
    ) as patch_get_df, mock.patch.object(
        process_functions, "checkColExist", return_value=True
    ) as patch_check_col_exist, mock.patch.object(
        validate, "check_values_between_two_df", return_value=("", "")
    ) as patch_check_values:
        input_df = None
        clin_class._cross_validate_assay_info_has_seq(clinicaldf=input_df)
        patch_check_values.assert_called_once_with(
            df1=None,
            df1_filename="clinical file",
            df1_id_to_check="SEQ_ASSAY_ID",
            df2=None,
            df2_filename="assay information",
            df2_id_to_check="SEQ_ASSAY_ID",
            ignore_case=True,
            allow_underscore=True,
        )


@pytest.mark.parametrize(
    "test_assay_df,expected_error,expected_warning",
    [
        (
            pd.DataFrame(
                dict(
                    SEQ_ASSAY_ID=["SAGE-1", "SAGE-2", "SAGE-3"],
                )
            ),
            "At least one SEQ_ASSAY_ID in your clinical file file does not exist as a SEQ_ASSAY_ID in your assay information file. "
            "Please update your file(s) to be consistent.\n",
            "",
        ),
        (
            pd.DataFrame(
                dict(
                    SEQ_ASSAY_ID=["SAGE-1-1", "SAGE-SAGE-1"],
                )
            ),
            "At least one SEQ_ASSAY_ID in your clinical file file does not exist as a SEQ_ASSAY_ID in your assay information file. "
            "Please update your file(s) to be consistent.\n",
            "",
        ),
        (
            pd.DataFrame(
                dict(
                    SEQ_ASSAY_ID=[
                        "SAGE-1-1",
                        "SAGE-SAGE-1",
                        "SAGE-1",
                        "SAGE-2",
                        "SAGE-3",
                    ],
                )
            ),
            "",
            "",
        ),
        (
            pd.DataFrame(
                dict(
                    SEQ_ASSAY_ID=[
                        "SAGE-1-1",
                        "SAGE-SAGE-1",
                        "SAGE-1",
                    ],
                )
            ),
            "",
            "",
        ),
    ],
    ids=["diff_ids", "semi_match_ids", "match_in_clinical", "exact_match"],
)
def test_that__cross_validate_assay_info_has_seq_returns_expected_msg_if_valid(
    clin_class, valid_clinical_df, test_assay_df, expected_warning, expected_error
):
    with mock.patch.object(
        validate,
        "parse_file_info_in_nested_list",
        return_value={
            "files": {"some_file"},
            "file_info": {"name": "assay_information.yaml", "path": ""},
        },
    ) as patch_assay_files, mock.patch.object(
        process_functions,
        "get_assay_dataframe",
        return_value=test_assay_df,
    ) as patch_get_df:
        errors, warnings = clin_class._cross_validate_assay_info_has_seq(
            valid_clinical_df
        )
        assert warnings == expected_warning
        assert errors == expected_error
