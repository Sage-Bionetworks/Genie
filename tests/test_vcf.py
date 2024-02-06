import pandas as pd
import pytest
from unittest.mock import mock_open, patch

from genie import transform
from genie_registry.vcf import vcf

OPEN_BUILTIN = "builtins.open"


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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
            "REF": ["AANT", "AACG", "AAAN"],
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
        "vcf: Your REF column has invalid allele values. "
        "This is the list of accepted allele values that can appear individually "
        "or in combination with each other: A,T,C,G,N.\n"
        "This is the list of accepted allele values that can only appear individually: \n"
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


def test_that__get_dataframe_throws_value_error_if_no_headers(vcf_class):
    file = "CHROM\tALT\tREF\n" "TEST\t3845\tNA"
    with patch(OPEN_BUILTIN, mock_open(read_data=file)):
        with pytest.raises(
            ValueError, match="Your vcf must start with the header #CHROM"
        ):
            vcf_class._get_dataframe(["some_path"])


def test_that__get_dataframe_reads_in_correct_nas(vcf_class):
    file = (
        "#CHROM\tALT\tREF\n"
        "TEST\t3845\tNA\n"
        "TEST\tNA\tnan\n"
        "TEST\t3846\tN/A\n"
        "NA\tnan\tNaN"
    )
    with patch(OPEN_BUILTIN, mock_open(read_data=file)):
        expected = pd.DataFrame(
            {
                "#CHROM": ["TEST", "TEST", "TEST", None],
                "ALT": ["3845", None, "3846", None],
                "REF": ["NA", "nan", None, "NaN"],
            }
        )
        vcf_df = vcf_class._get_dataframe(["some_path"])
        pd.testing.assert_frame_equal(vcf_df, expected)


@pytest.mark.parametrize(
    "input,expected_columns",
    [
        (
            pd.DataFrame(
                {
                    "#CHROM": ["TEST"],
                    "ALT": ["3845"],
                    "REF": ["NA"],
                }
            ),
            ["#CHROM", "ALT"],
        ),
        (
            pd.DataFrame(
                {
                    "#CHROM": ["TEST"],
                    "ALT": ["3845"],
                    "rEf": ["NA"],
                }
            ),
            ["#CHROM", "ALT", "rEf"],
        ),
    ],
    ids=["with_allele_col", "no_allele_col"],
)
def test_that__get_dataframe_uses_correct_columns_to_replace(
    vcf_class, input, expected_columns
):
    file = "#CHROM\tALT\tref\n" "TEST\t3845\tNA"
    with patch(OPEN_BUILTIN, mock_open(read_data=file)), patch.object(
        pd, "read_csv", return_value=input
    ), patch.object(transform, "_convert_values_to_na") as patch_convert_to_na:
        vcf_class._get_dataframe(["some_path"])
        patch_convert_to_na.assert_called_once_with(
            input_df=input,
            values_to_replace=["NA", "nan", "NaN"],
            columns_to_convert=expected_columns,
        )
