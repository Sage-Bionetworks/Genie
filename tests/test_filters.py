"""
Test GENIE filters
"""
import datetime
import os
import sys
from unittest.mock import patch

import pandas as pd

from genie.process_functions import seqDateFilter
from genie import database_to_staging
from genie.database_to_staging import (
    seq_assay_id_filter,
    redact_phi,
    no_genepanel_filter,
    _to_redact_interval,
    _redact_year,
    _to_redact_difference,
)
from genie.consortium_to_public import commonVariantFilter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, "../../genie"))


# syn = synapseclient.login()
def test_seqassayfilter():
    # SEQ ASSAY ID filter, only seq assay ids with more than
    sagetest = ["SAGE-TEST-1"] * 51
    sagetest.extend(["SAGE-TEST-2"] * 10)

    clinicalDf = pd.DataFrame(sagetest)
    clinicalDf.rename(columns={0: "SEQ_ASSAY_ID"}, inplace=True)
    clinicalDf["SAMPLE_ID"] = range(0, len(sagetest))
    # Make sure that only samples of seq assay ids that have 50 or
    # more samples are returned
    samples = seq_assay_id_filter(clinicalDf).reset_index()
    # Must make sure that seq assays have the same index to do comparisons
    del samples["index"]
    samples = samples["SAMPLE_ID"]
    assert all(samples == pd.Series(range(51, 61)))


def test_seqdatefilter():
    # Check that seq date filter works correctly
    SEQ_DATE = "Jan-2018"
    processingDate = datetime.datetime.strptime(SEQ_DATE, "%b-%Y")
    clinicalDf = pd.DataFrame(
        [
            "SAGE-TEST-1",
            "SAGE-TEST-2",
            "SAGE-TEST-3",
            "SAGE-TEST-4",
            "SAGE-TEST-5",
            "SAGE-TEST-6",
            "SAGE-TEST-7",
        ]
    )
    clinicalDf.rename(columns={0: "SAMPLE_ID"}, inplace=True)
    expected = [
        "SAGE-TEST-1",
        "SAGE-TEST-2",
        "SAGE-TEST-3",
        "SAGE-TEST-4",
        "SAGE-TEST-5",
        "SAGE-TEST-6",
    ]
    clinicalDf["SEQ_DATE"] = [
        "Jan-2018",
        "Apr-2018",
        "Oct-2017",
        "Jul-2017",
        "Apr-2017",
        "Jan-2017",
        "Oct-2016",
    ]
    samples = seqDateFilter(clinicalDf, processingDate, 366)
    assert all(samples == expected)

    # Half year filter
    samples = seqDateFilter(clinicalDf, processingDate, 184)
    expected = ["SAGE-TEST-1", "SAGE-TEST-2", "SAGE-TEST-3", "SAGE-TEST-4"]
    assert all(samples == expected)

    # Test leap year
    SEQ_DATE = "Jan-2017"
    processingDate = datetime.datetime.strptime(SEQ_DATE, "%b-%Y")
    clinicalDf = pd.DataFrame(
        [
            "SAGE-TEST-1",
            "SAGE-TEST-2",
            "SAGE-TEST-3",
            "SAGE-TEST-4",
            "SAGE-TEST-5",
            "SAGE-TEST-6",
            "SAGE-TEST-7",
        ]
    )
    clinicalDf.rename(columns={0: "SAMPLE_ID"}, inplace=True)
    clinicalDf["SEQ_DATE"] = [
        "Jan-2017",
        "Apr-2017",
        "Oct-2016",
        "Jul-2016",
        "Apr-2016",
        "Jan-2016",
        "Oct-2015",
    ]
    expected = [
        "SAGE-TEST-1",
        "SAGE-TEST-2",
        "SAGE-TEST-3",
        "SAGE-TEST-4",
        "SAGE-TEST-5",
        "SAGE-TEST-6",
    ]
    samples = seqDateFilter(clinicalDf, processingDate, 366)
    assert all(samples == expected)

    # Check half year release date during leap year
    SEQ_DATE = "Jul-2016"
    processingDate = datetime.datetime.strptime(SEQ_DATE, "%b-%Y")
    samples = seqDateFilter(clinicalDf, processingDate, 184)
    assert all(samples == expected)


def test__to_redact_interval():
    """Tests the correct boolean vectors are returned for phi and pediatric
    redaction"""
    values = pd.Series([32850, 32485, 6570, 6569, "<foo", ">testing"])
    to_redact, to_redact_peds = _to_redact_interval(values)
    expected_redact = [True, False, False, False, False, True]
    expected_redact_peds = [False, False, False, True, True, False]
    assert to_redact.to_list() == expected_redact
    assert to_redact_peds.to_list() == expected_redact_peds


def test__redact_year():
    """Tests redaction of birth year based on < and >"""
    values = pd.Series([1923, 2003, "<foo", ">testing"])
    redacted = _redact_year(values)
    expected_redact = [1923, 2003, "withheld", "cannotReleaseHIPAA"]
    assert redacted.to_list() == expected_redact


def test___to_redact_difference():
    """Tests if a difference between two year columns is >89, redact"""
    year1 = pd.Series([1923, 2000, float("nan")])
    year2 = pd.Series([1926, 2100, 2000])
    redacted = _to_redact_difference(year1, year2)
    expected_redact = [False, True, False]
    assert redacted.to_list() == expected_redact


def test_redact_phi():
    """Redacts PHI interval and years"""
    return_bools = ([True, False, False], [False, False, True])
    clinicaldf = pd.DataFrame(["SAGE-TEST-1", "SAGE-TEST-2", "SAGE-TEST-3"])
    clinicaldf.rename(columns={0: "SAMPLE_ID"}, inplace=True)
    clinicaldf["AGE_AT_SEQ_REPORT"] = [32850, 32485, 6570]
    clinicaldf["BIRTH_YEAR"] = [1900, "<1900", 1902]
    # These are the years that are returned by the _redact_year
    # Use these against YEAR_CONTACT and YEAR_DEATH to calculate
    # expected_birth
    return_year = pd.Series([2000, 1903, 1904])
    clinicaldf["YEAR_CONTACT"] = [2100, 1904, float("nan")]
    clinicaldf["YEAR_DEATH"] = [2001, float("nan"), 2000]

    expected_age = pd.Series([">32485", 32485, "<6570"])
    expected_birth = pd.Series(["cannotReleaseHIPAA", 1903, "cannotReleaseHIPAA"])
    with patch.object(
        database_to_staging, "_to_redact_interval", return_value=return_bools
    ), patch.object(database_to_staging, "_redact_year", return_value=return_year):
        redacted_clin = redact_phi(
            clinicaldf, interval_cols_to_redact=["AGE_AT_SEQ_REPORT"]
        )
        assert all(redacted_clin["AGE_AT_SEQ_REPORT"] == expected_age)
        assert all(redacted_clin["BIRTH_YEAR"] == expected_birth)


# def test_MAFinBED():
#   syn = mock.create_autospec(synapseclient.Synapse)
#   # MAF in BED filter (Make sure that the only Hugo symbol that
# passes through this filter is the ERRFI1 gene)
#   databaseSynIdMapping = syn.tableQuery('select * from syn11600968')
#   databaseSynIdMappingDf = databaseSynIdMapping.asDataFrame()
# CENTER_MAPPING = syn.tableQuery(
#     'SELECT * FROM syn11601248 where stagingSynId is not null')
#   CENTER_MAPPING_DF = CENTER_MAPPING.asDataFrame()
#   mafPath = os.path.join(SCRIPT_DIR,"testingmaf.txt")
# temp = runMAFinBED(
#     syn, mafPath, CENTER_MAPPING_DF, databaseSynIdMappingDf, test=True)
#   maf = pd.read_csv(mafPath,sep="\t")
#   assert all(maf['Hugo_Symbol'][maf['inBED'] ==  True].unique() == "ERRFI1")
#   os.unlink(mafPath)

# def test_MutationInCis():
#   def createMockTable(dataframe):
#       table = mock.create_autospec(synapseclient.table.CsvFileTable)
#       table.asDataFrame.return_value= dataframe
#       return(table)

#   def table_query_results(*args):
#       return(table_query_results_map[args])

#   databaseMapping = pd.DataFrame(dict(Database=['mutationsInCis'],
#                                       Id=['syn7765462']))

#   centerMapping = pd.DataFrame(dict(center=['SAGE'],
#                                     inputSynId=['syn11601335'],
#                                     stagingSynId=['syn11601337']))

#   mutCisDf = pd.DataFrame(dict(Flag=['mutationsInCis'],
#                                  Center=['syn7765462'],
#                                  Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
#                                  Variant_Classification=["Nonsense_Mutation"],
#                                  Hugo_Symbol=["AKT1"],
#                                  HGVSp_Short=["p.1234"],
#                                  Chromosome=["1"],
#                                  Start_Position=[324234],
#                                  Reference_Allele=["AGC"],
#                                  Tumor_Seq_Allele2=["GCCCT"],
#                                  t_alt_count_num=[3],
#                                  t_depth=[234]))

#   #This is the gene positions that all bed
#   dataframe will be processed against
# table_query_results_map = {
#     ('select * from syn11600968',):
#         createMockTable(databaseMapping),
#     ('SELECT * FROM syn11601248 where staginsgSynId is not null',):
#         createMockTable(centerMapping),
#     ("select * from syn7765462 where Center = 'SAGE'",):
#         createMockTable(mutCisDf)
# }

#   syn = mock.create_autospec(synapseclient.Synapse)
#   syn.tableQuery.side_effect=table_query_results
#   # Mutation in cis filter check (There is one TOSS sample,
#   it should be GENIE-TEST-0-1)
# remove_samples = mutation_in_cis_filter(
#     syn, True, None, None, "syn11601206", None)
#   assert  all(remove_samples == "GENIE-TEST-0-1")


def test_commonvariantfilter():
    """
    Test common variant filter, make sure none
    of the common_variant rows are kept
    """
    mutationDf = pd.DataFrame(
        [
            "common_variant",
            "test1",
            "asdfasd:common_variant",
            "fo;common_variant",
            "test2",
            "test3",
        ]
    )
    mutationDf.rename(columns={0: "FILTER"}, inplace=True)
    maf = commonVariantFilter(mutationDf)
    expected = ["test1", "test2", "test3"]
    assert all(maf["FILTER"] == expected)


def test_no_genepanel_filter():
    """
    Tests the filter to remove samples
    with no bed files
    """
    sagetest = ["SAGE-TEST-1"] * 51
    beddf = pd.DataFrame(sagetest)
    beddf.rename(columns={0: "SEQ_ASSAY_ID"}, inplace=True)
    sagetest.extend(["SAGE-TEST-2"] * 10)

    clinicaldf = pd.DataFrame(sagetest)
    clinicaldf.rename(columns={0: "SEQ_ASSAY_ID"}, inplace=True)
    clinicaldf["SAMPLE_ID"] = sagetest

    remove_samples = no_genepanel_filter(clinicaldf, beddf)
    assert all(remove_samples == ["SAGE-TEST-2"] * 10)
