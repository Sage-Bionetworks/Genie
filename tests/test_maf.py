from unittest import mock
from unittest.mock import patch

import pandas as pd
import pytest
import synapseclient

import genie_registry.maf
from genie_registry.maf import maf
from genie_registry.mafSP import mafSP

syn = mock.create_autospec(synapseclient.Synapse)

maf_class = maf(syn, "SAGE")
mafsp_class = mafSP(syn, "SAGE")


def test_invalidname_validateFilename():
    with pytest.raises(AssertionError):
        maf_class.validateFilename(['foo'])


def test_valid_validateFilename():
    assert maf_class.validateFilename([
        "data_mutations_extended_SAGE.txt"]) == "maf"


def test_perfect_validation():
    mafDf = pd.DataFrame(dict(
        CHROMOSOME=[1, 2, 3, 4, 5],
        START_POSITION=[1, 2, 3, 4, 2],
        REFERENCE_ALLELE=["A", "A", "A", "A", "A"],
        TUMOR_SAMPLE_BARCODE=["ID1-1", "ID1-1", "ID1-1", "ID1-1", "ID1-1"],
        T_ALT_COUNT=[1, 2, 3, 4, 3],
        T_DEPTH=[1, 2, 3, 4, 3],
        T_REF_COUNT=[1, 2, 3, 4, 3],
        N_DEPTH=[1, 2, 3, 4, 3],
        N_REF_COUNT=[1, 2, 3, 4, 3],
        N_ALT_COUNT=[1, 2, 3, 4, 3],
        TUMOR_SEQ_ALLELE2=["A", "A", "A", "A", "A"]))

    error, warning = maf_class._validate(mafDf)
    assert error == ""
    assert warning == ""
    error, warning = mafsp_class._validate(mafDf)
    assert error == ""
    assert warning == ""


def test_firstcolumn_validation():
    """Tests if first column isn't correct"""
    mafDf = pd.DataFrame({
        'REFERENCE_ALLELE': ["A", "B", "C", "D", "E"],
        "START_POSITION": [1, 2, 3, 4, 2],
        "CHROMOSOME": ["A", "A", "A", "A", "A"],
        "TUMOR_SAMPLE_BARCODE": ["ID1-1", "ID1-1", "ID1-1", "ID1-1", "ID1-1"],
        "T_ALT_COUNT": [1, 2, 3, 4, 3],
        "T_DEPTH": [1, 2, 3, 4, 3],
        "T_REF_COUNT": [1, 2, 3, 4, 3],
        "N_DEPTH": [1, 2, 3, 4, 3],
        "N_REF_COUNT": [1, 2, 3, 4, 3],
        "N_ALT_COUNT": [1, 2, 3, 4, 3],
        "TUMOR_SEQ_ALLELE2": ["A", "A", "A", "A", "A"]})
    order = ['REFERENCE_ALLELE', 'START_POSITION', 'CHROMOSOME',
             'TUMOR_SAMPLE_BARCODE', 'T_ALT_COUNT', 'T_DEPTH',
             'T_REF_COUNT', 'N_DEPTH', 'N_REF_COUNT', 'N_ALT_COUNT',
             'TUMOR_SEQ_ALLELE2']
    error, warning = maf_class._validate(mafDf[order])
    expectedErrors = ("Mutation File: First column header must be "
                      "one of these: CHROMOSOME, HUGO_SYMBOL, "
                      "TUMOR_SAMPLE_BARCODE.\n")
    assert error == expectedErrors
    assert warning == ""


def test_missingcols_validation():
    """Tests missing columns"""
    emptydf = pd.DataFrame()
    error, warning = maf_class._validate(emptydf)
    expected_errors = ("Mutation File: Must at least have these headers: "
                       "CHROMOSOME,START_POSITION,REFERENCE_ALLELE,"
                       "TUMOR_SAMPLE_BARCODE,T_ALT_COUNT,"
                       "TUMOR_SEQ_ALLELE2. "
                       "If you are writing your maf file with R, please make"
                       "sure to specify the 'quote=FALSE' parameter.\n"
                       "Mutation File: If you are missing T_DEPTH, "
                       "you must have T_REF_COUNT!\n")
    expected_warnings = ("Mutation File: Does not have the column headers "
                         "that can give extra information to the processed "
                         "mutation file: T_REF_COUNT, N_DEPTH, N_REF_COUNT, "
                         "N_ALT_COUNT.\n")
    assert error == expected_errors
    assert warning == expected_warnings


def test_errors_validation():
    mafDf = pd.DataFrame(dict(
        CHROMOSOME=[1, "chr2", "WT", 4, 5],
        START_POSITION=[1, 2, 3, 4, 2],
        REFERENCE_ALLELE=["NA", float('nan'), "A", "A", "A"],
        TUMOR_SAMPLE_BARCODE=["ID1-1", "ID1-1", "ID1-1", "ID1-1", "ID1-1"],
        TUMOR_SEQ_ALLELE2=["B", "B", "B", "B", "B"],
        T_ALT_COUNT=[1, 2, 3, 4, 3],
        T_DEPTH=[1, 2, 3, 4, 3],
        N_REF_COUNT=[1, 2, 3, 4, 3],
        N_ALT_COUNT=[1, 2, 3, 4, 3]))

    error, warning = maf_class._validate(mafDf)

    expectedErrors = (
        "Mutation File: "
        "REFERENCE_ALLELE can't have any blank or null values.\n"
        "Mutation File: "
        "CHROMOSOME column cannot have any values that start "
        "with 'chr' or any 'WT' values.\n"
    )
    expectedWarnings = ("Mutation File: "
                        "Does not have the column headers that can give "
                        "extra information to the processed mutation file: "
                        "T_REF_COUNT, N_DEPTH.\n"
                        "Mutation File: "
                        "REFERENCE_ALLELE column contains 'NA' values, "
                        "which cannot be placeholders for blank values.  "
                        "Please put in empty strings for blank values.\n")

    assert error == expectedErrors
    assert warning == expectedWarnings


def test_invalid_validation():
    mafDf = pd.DataFrame(dict(
        CHROMOSOME=[1, 2, 3, 4, 2, 4],
        T_ALT_COUNT=[1, 2, 3, 4, 3, 4],
        START_POSITION=[1, 2, 3, 4, 2, 4],
        REFERENCE_ALLELE=["A", "A", "A", "A", "A", "A"],
        TUMOR_SAMPLE_BARCODE=["ID1-1", "ID1-1", "ID1-1", "ID1-1", "ID1-1", "ID1-1"],
        N_DEPTH=[1, 2, 3, 4, 3, 4],
        N_REF_COUNT=[1, 2, 3, 4, 3, 4],
        N_ALT_COUNT=[1, 2, 3, 4, 3, 4],
        TUMOR_SEQ_ALLELE2=["NA", float('nan'), "A", "A", "A", "A"]))

    with patch.object(genie_registry.maf, "_check_tsa1_tsa2",
                      return_value="") as check_tsa1_tsa2:
        error, warning = maf_class._validate(mafDf)
        check_tsa1_tsa2.assert_called_once_with(mafDf)
    expectedErrors = (
        "Mutation File: Should not have duplicate rows\n"
        "Mutation File: "
        "If you are missing T_DEPTH, you must have T_REF_COUNT!\n"
        "Mutation File: "
        "TUMOR_SEQ_ALLELE2 can't have any blank or null values.\n"
    )
    expectedWarnings = (
        "Mutation File: TUMOR_SEQ_ALLELE2 column contains 'NA' values, "
        "which cannot be placeholders for blank values.  "
        "Please put in empty strings for blank values.\n"
        "Mutation File: Does not have the column headers that can give "
        "extra information to the processed mutation file: T_REF_COUNT.\n")
    assert error == expectedErrors
    assert warning == expectedWarnings


@pytest.mark.parametrize("col", ["temp", "REFERENCE_ALLELE"])
def test_noerror__check_allele_col(col):
    """Test error and warning is an empty string if REF col isn't passed in"""
    df = pd.DataFrame(dict(
        REFERENCE_ALLELE=["A", "A"]
    ))
    error, warning = genie_registry.maf._check_allele_col(df, col)
    assert error == ""
    assert warning == ""


def test_warning__check_allele_col():
    """Test warning occurs when 'NA' string is passed in"""
    df = pd.DataFrame(dict(
        TEMP=["NA", "A"]
    ))
    error, warning = genie_registry.maf._check_allele_col(df, "TEMP")
    assert error == ""
    assert warning == (
        "Mutation File: "
        "TEMP column contains 'NA' values, "
        "which cannot be placeholders for blank values.  "
        "Please put in empty strings for blank values.\n"
    )


def test_error__check_allele_col():
    """Test error occurs when blank allele is passed in"""
    df = pd.DataFrame(dict(
        TEMP=[float('nan'), "A"]
    ))
    error, warning = genie_registry.maf._check_allele_col(df, "TEMP")
    assert error == (
        "Mutation File: TEMP can't have any blank or null values.\n"
    )
    assert warning == ""


def test_invalid__check_tsa1_tsa2():
    """Test the scenario in which maf file has TSA1 and TSA2 and fails"""
    df = pd.DataFrame(dict(
        REFERENCE_ALLELE=["A", "A", "A"],
        TUMOR_SEQ_ALLELE1=["B", "B", "B"],
        TUMOR_SEQ_ALLELE2=["C", "C", "C"]
    ))
    error = genie_registry.maf._check_tsa1_tsa2(df)
    assert error == (
        "Mutation File: Contains both "
        "TUMOR_SEQ_ALLELE1 and TUMOR_SEQ_ALLELE2 columns. "
        "The values in TUMOR_SEQ_ALLELE1 must be the same as "
        "all the values in REFERENCE_ALELLE OR TUMOR_SEQ_ALLELE2."
    )


@pytest.mark.parametrize(
    "df",
    [pd.DataFrame(dict(REFERENCE_ALLELE=["A", "A", "A"],
                       TUMOR_SEQ_ALLELE1=["C", "C", "C"],
                       TUMOR_SEQ_ALLELE2=["C", "C", "C"])),
     pd.DataFrame(dict(REFERENCE_ALLELE=["C", "C", "C"],
                       TUMOR_SEQ_ALLELE1=["C", "C", "C"],
                       TUMOR_SEQ_ALLELE2=["A", "A", "A"]))]
)
def test_valid__check_tsa1_tsa2(df):
    """Test valid TSA1 and TSA2"""
    error = genie_registry.maf._check_tsa1_tsa2(df)
    assert error == ""
