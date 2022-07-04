from unittest import mock
from unittest.mock import patch

import pandas as pd
import pytest
import synapseclient

from genie_registry.mutationsInCis import mutationsInCis


def createMockTable(dataframe):
    table = mock.create_autospec(synapseclient.table.CsvFileTable)
    table.asDataFrame.return_value = dataframe
    return table


def table_query_results(*args):
    return table_query_results_map[args]


databaseMapping = pd.DataFrame(dict(Database=["mutationsInCis"], Id=["syn7765462"]))
mutCisDf = pd.DataFrame(
    dict(
        Flag=["mutationsInCis"],
        Center=["SAGE"],
        Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
        Variant_Classification=["Nonsense_Mutation"],
        Hugo_Symbol=["AKT1"],
        HGVSp_Short=["p.1234"],
        Chromosome=["1"],
        Start_Position=[324234],
        Reference_Allele=["AGC"],
        Tumor_Seq_Allele2=["GCCCT"],
        t_alt_count_num=[3],
        t_depth=[234],
    )
)

# This is the gene positions that all bed dataframe will be processed against
table_query_results_map = {
    ("SELECT * FROM syn10967259",): createMockTable(databaseMapping),
    ("select * from syn11601206 where Center = 'SAGE'",): createMockTable(mutCisDf),
}

syn = mock.create_autospec(synapseclient.Synapse)
syn.tableQuery.side_effect = table_query_results


@pytest.fixture
def mutCis(genie_config):
    return mutationsInCis(syn, "SAGE", genie_config)


ENTITY = synapseclient.Project("testing", annotations={"dbMapping": ["syn10967259"]})


def test_validation__perfect(mutCis):
    with pytest.raises(AssertionError):
        mutCis.validateFilename(["foo"])
    assert (
        mutCis.validateFilename(["mutationsInCis_filtered_samples.csv"])
        == "mutationsInCis"
    )

    mutCisDf = pd.DataFrame(
        dict(
            Flag=["mutationsInCis"],
            Center=["SAGE"],
            Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
            Variant_Classification=["Nonsense_Mutation"],
            Hugo_Symbol=["AKT1"],
            HGVSp_Short=["p.1234"],
            Chromosome=["1"],
            Start_Position=[324234],
            Reference_Allele=["AGC"],
            Tumor_Seq_Allele2=["GCCCT"],
            t_alt_count_num=[3],
            t_depth=[234],
        )
    )
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = mutCis._validate(mutCisDf)
        assert error == ""
        assert warning == ""


def test_validation__missing_col(mutCis):
    mutCisDf = pd.DataFrame(
        dict(
            Flag=["mutationsInCis"],
            Center=["SAGE"],
            Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
            Hugo_Symbol=["AKT1"],
            HGVSp_Short=["p.1234"],
            Chromosome=["1"],
            Reference_Allele=["AGC"],
            Tumor_Seq_Allele2=["GCCCT"],
            t_alt_count_num=[3],
            t_depth=[234],
        )
    )
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = mutCis._validate(mutCisDf)
        expectedErrors = (
            "Mutations In Cis Filter File: Must at least have these headers: "
            "Variant_Classification,Start_Position.\n"
        )
        assert error == expectedErrors


def test_validation__invalid(mutCis):
    mutCisDf = pd.DataFrame(
        dict(
            Flag=["mutationsInCis"],
            Center=["SAGE"],
            Tumor_Sample_Barcode=["ID1-1"],
            Variant_Classification=["Nonsense_Mutation"],
            Hugo_Symbol=["AKT1"],
            HGVSp_Short=["foo"],
            Chromosome=["1"],
            Start_Position=[324234],
            Reference_Allele=["AGC"],
            Tumor_Seq_Allele2=["GCCCT"],
            t_alt_count_num=[3],
            t_depth=[234],
        )
    )
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = mutCis._validate(mutCisDf)
        expectedErrors = (
            "Mutations In Cis Filter File: All variants must come from "
            "the original mutationInCis_filtered_samples.csv file in each "
            "institution's staging folder.\n"
            "Mutations In Cis Filter File: TUMOR_SAMPLE_BARCODE must "
            "start with GENIE-SAGE\n"
        )
        assert error == expectedErrors
