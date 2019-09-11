import mock
import pytest

import pandas as pd
import synapseclient

from genie.vcf import vcf

syn = mock.create_autospec(synapseclient.Synapse)

vcfClass = vcf(syn, "SAGE")


def test_processing():
    pass


@pytest.fixture(params=[
    (["foo"]),
    (["GENIE-SAGE-ID1.bed"])
    ])
def filename_fileformat_map(request):
    return request.param


def test_incorrect_validatefilename(filename_fileformat_map):
    filepath_list = filename_fileformat_map
    with pytest.raises(AssertionError):
        vcfClass.validateFilename(filepath_list)


def test_validation():
    assert vcfClass.validateFilename(["GENIE-SAGE-ID1.vcf"]) == "vcf"

    vcfDf = pd.DataFrame({
        "#CHROM": ['2', '9', '12'],
        "POS": [69688533, 99401860, 53701241],
        "ID": ['AAK1', 'AAED1', 'AAAS'],
        "REF": ['AAK1', 'AAED1', 'AAAS'],
        "ALT": ['AAK1', 'AAED1', 'AAAS'],
        "QUAL": ['AAK1', 'AAED1', 'AAAS'],
        "FILTER": ['AAK1', 'AAED1', 'AAAS'],
        "INFO": ['AAK1', 'AAED1', 'AAAS']})

    error, warning = vcfClass._validate(vcfDf)
    assert error == ""
    assert warning == ""

    vcfDf = pd.DataFrame({
        "POS": [69688533, 99401860, 53701241],
        "ID": ['AAK1', 'AAED1', 'AAAS'],
        "REF": ['AAK1', 'AAED1', 'AAAS'],
        "ALT": ['AAK1', 'AAED1', 'AAAS'],
        "QUAL": ['AAK1', 'AAED1', 'AAAS'],
        "FILTER": ['AAK1', 'AAED1', 'AAAS'],
        "INFO": ['AAK1', 'AAED1', 'AAAS'],
        "FOO": ['AAK1', 'AAED1', 'AAAS'],
        "DOO": ['AAK1', 'AA ED1', 'AAAS']})

    error, warning = vcfClass._validate(vcfDf)
    expectedError = (
        "Your vcf file must have these headers: CHROM, POS, ID, REF, "
        "ALT, QUAL, FILTER, INFO.\n"
        "Your vcf file must have FORMAT header if genotype columns exist.\n")
    expectedWarning = (
        "Your vcf file should not have any white spaces "
        "in any of the columns.\n")
    assert error == expectedError
    assert warning == expectedWarning

    vcfDf = pd.DataFrame({
        "#CHROM": ['chr2', 'chrM', '12'],
        "POS": [69688533, 99401860, 53701241],
        "ID": ['AAK1', 'AAED1', 'AAAS'],
        "REF": ['AAK1', 'AAED1', 'AAAS'],
        "ALT": ['AAK1', 'AAED1', 'AAAS'],
        "QUAL": ['AAK1', 'AAED1', 'AAAS'],
        "FILTER": ['AAK1', 'AAED1', 'AAAS'],
        "INFO": ['AAK1', 'AAED1', 'AAAS']})

    error, warning = vcfClass._validate(vcfDf)
    expectedError = ("Your vcf file must not have variants on chrM.\n")
    expectedWarning = (
        "Your vcf file should not have the chr prefix "
        "in front of chromosomes.\n")
    assert error == expectedError
    assert warning == expectedWarning
