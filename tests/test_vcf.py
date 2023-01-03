import pandas as pd
import pytest

from genie_registry.vcf import vcf


@pytest.fixture
def vcf_class(syn):
    return vcf(syn, "SAGE")


def test_processing():
    pass


@pytest.mark.parametrize("filepath_list", [["foo"], ["GENIE-SAGE-ID1.bed"]])
def test_incorrect_validatefilename(vcf_class, filepath_list):
    with pytest.raises(AssertionError):
        vcf_class.validateFilename(filepath_list)


def test__validate_filename(vcf_class):
    assert vcf_class.validateFilename(["GENIE-SAGE-ID1.vcf"]) == "vcf"


def test_validation_valid_no_samples(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == ""
    assert warning == ""


def test_validation_valid_one_sample_tumor(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "FORMAT": ["AG", "AG", "AG"],
            "TUMOR": ["AG", "AG", "AG"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == ""
    assert warning == ""


def test_validation_valid_one_sample(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "FORMAT": ["AG", "AG", "AG"],
            "GENIE-SAGE-1-1": ["AG", "AG", "AG"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == ""
    assert warning == ""


def test_validation_missing_format_col(vcf_class):
    """vcf must have FORMAT column is greater than 8 columns"""
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "GENIE-SAGE-1-1": ["AG", "AG", "AG"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == "vcf: Must have FORMAT header if sample columns exist.\n"
    assert warning == ""


def test_validation_invalid_one_sample(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "FORMAT": ["AG", "AG", "AG"],
            "GENIE-SAE-1-1": ["AG", "AG", "AG"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == (
        "vcf: tumor sample column must start with GENIE-SAGE if vcf represents "
        "a single sample and TUMOR is not the sample column header.\n"
    )
    assert warning == ""


def test_validation_valid_two_samples(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "FORMAT": ["AG", "AG", "AG"],
            "GENIE-SAGE-1-1": ["AG", "AG", "AG"],
            "GENIE-SAGE-1-2": ["AG", "AG", "AG"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == ""
    assert warning == ""


def test_validation_invalid_two_samples_tumor(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "FORMAT": ["AG", "AG", "AG"],
            "GENIE-SAE-1-1": ["AG", "AG", "AG"],
            "GENIE-SAGE-1-2": ["AG", "AG", "AG"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == "vcf: tumor sample column must start with GENIE-SAGE\n"
    assert warning == ""


def test_validation_invalid_two_samples_normal(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "FORMAT": ["AG", "AG", "AG"],
            "GENIE-SAGE-1-1": ["AG", "AG", "AG"],
            "GENIE-SAE-1-2": ["AG", "AG", "AG"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    assert error == "vcf: normal sample column must start with GENIE-SAGE\n"
    assert warning == ""


def test_validation_invalid_white_space(vcf_class):
    vcfDf = pd.DataFrame(
        {
            "#CHROMM": ["2", "9", "12"],
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AA ED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    expectedError = (
        "vcf: Must have these headers: CHROM, POS, ID, REF, "
        "ALT, QUAL, FILTER, INFO.\n"
    )
    expectedWarning = "vcf: Should not have any white spaces in any of the columns.\n"
    assert error == expectedError
    assert warning == expectedWarning


def test_validation_invalid_content(vcf_class):
    """CHeck for duplicated variants and invalid chromosome values"""
    vcfDf = pd.DataFrame(
        {
            "#CHROM": ["chr2", "chrM", float("nan"), "chr2"],
            "POS": [69688533, 99401860, 53701241, 69688533],
            "ID": ["AAK1", "AAED1", "AAAS", "AAK1"],
            "REF": ["AAK1", "AAED1", "AAAS", "AAK1"],
            "ALT": ["AAK1", "AAED1", "AAAS", "AAK1"],
            "QUAL": ["AAK1", "AAED1", "AAAS", "AAK1"],
            "FILTER": ["AAK1", "AAED1", "AAAS", "AAK1"],
            "INFO": ["AAK1", "AAED1", "AAAS", "AAK1"],
        }
    )
    error, warning = vcf_class._validate(vcfDf)
    expectedError = (
        "vcf: Must not have duplicate variants.\n"
        "vcf: May contain rows that are "
        "space delimited instead of tab delimited.\n"
        "vcf: Please double check your #CHROM column.  This column must only be these values: "
        "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT\n"
    )
    expectedWarning = "vcf: Should not have the chr prefix in front of chromosomes.\n"
    assert error == expectedError
    assert warning == expectedWarning


def test_validation_more_than_11_cols(vcf_class):
    """Test that a vcf with more than 11 columns fails and required headers"""
    larger_df = pd.DataFrame(columns=list(range(0, 12)))
    error, warning = vcf_class._validate(larger_df)
    assert error == (
        "vcf: Must have these headers: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.\n"
        "vcf: Should not have more than 11 columns. "
        "Only single sample or matched tumor normal vcf files are accepted.\n"
    )
    assert warning == ""
