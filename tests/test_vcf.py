from unittest import mock

import pandas as pd
import pytest
import synapseclient

from genie_registry.vcf import vcf

syn = mock.create_autospec(synapseclient.Synapse)

vcfClass = vcf(syn, "SAGE")


def test_processing():
    pass


@pytest.mark.parametrize("filepath_list", [["foo"], ["GENIE-SAGE-ID1.bed"]])
def test_incorrect_validatefilename(filepath_list):
    with pytest.raises(AssertionError):
        vcfClass.validateFilename(filepath_list)


def test_validation():
    assert vcfClass.validateFilename(["GENIE-SAGE-ID1.vcf"]) == "vcf"

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

    error, warning = vcfClass._validate(vcfDf)
    assert error == ""
    assert warning == ""

    vcfDf = pd.DataFrame(
        {
            "POS": [69688533, 99401860, 53701241],
            "ID": ["AAK1", "AAED1", "AAAS"],
            "REF": ["AAK1", "AAED1", "AAAS"],
            "ALT": ["AAK1", "AAED1", "AAAS"],
            "QUAL": ["AAK1", "AAED1", "AAAS"],
            "FILTER": ["AAK1", "AAED1", "AAAS"],
            "INFO": ["AAK1", "AAED1", "AAAS"],
            "FOO": ["AAK1", "AAED1", "AAAS"],
            "DOO": ["AAK1", "AA ED1", "AAAS"],
        }
    )

    error, warning = vcfClass._validate(vcfDf)
    expectedError = (
        "vcf: Must have these headers: CHROM, POS, ID, REF, "
        "ALT, QUAL, FILTER, INFO.\n"
        "vcf: Must have FORMAT header if genotype columns exist.\n"
    )
    expectedWarning = (
        "vcf: Should not have any white spaces " "in any of the columns.\n"
    )
    assert error == expectedError
    assert warning == expectedWarning

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

    error, warning = vcfClass._validate(vcfDf)
    expectedError = (
        "vcf: Must not have duplicate variants.\n"
        "vcf: May contain rows that are "
        "space delimited instead of tab delimited.\n"
        "vcf: Must not have variants on chrM.\n"
    )
    expectedWarning = (
        "vcf: Should not have the chr prefix " "in front of chromosomes.\n"
    )
    assert error == expectedError
    assert warning == expectedWarning
