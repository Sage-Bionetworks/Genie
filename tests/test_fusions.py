from unittest import mock
from unittest.mock import patch

import pandas as pd
import pytest
import synapseclient

from genie_registry.fusions import fusions


def createMockTable(dataframe):
    table = mock.create_autospec(synapseclient.table.CsvFileTable)
    table.asDataFrame.return_value = dataframe
    return(table)


def table_query_results(*args):
    return(table_query_results_map[args])


databaseMapping = pd.DataFrame(dict(Database=['bed'],
                                    Id=['syn8457748']))
symbols = pd.DataFrame(dict(Hugo_Symbol=['AAK1', 'AAED1', 'AAAS'],
                            ID=['AAK1', 'AAED', 'AAAS']))

# This is the gene positions that all bed dataframe will be processed against
table_query_results_map = {
    ("SELECT * FROM syn10967259",): createMockTable(databaseMapping),
    ("select Hugo_Symbol, ID from syn8457748 where CENTER = 'SAGE'",):
        createMockTable(symbols),
}
ENTITY = synapseclient.Project("testing",
                               annotations={'dbMapping': ["syn10967259"]})

syn = mock.create_autospec(synapseclient.Synapse)
syn.tableQuery.side_effect = table_query_results
fusionClass = fusions(syn, "SAGE")


def test_processing():
    expectedFusionDf = pd.DataFrame({
        "HUGO_SYMBOL": ['AAED1', 'AAK1', 'AAAS5'],
        "ENTREZ_GENE_ID": [0, 0, 0],
        "CENTER": ["SAGE", "SAGE", "SAGE"],
        "TUMOR_SAMPLE_BARCODE": ["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID2-1",
                                 "GENIE-SAGE-ID1-3"],
        "FUSION": ["AAED1-AAK1", "AAAS-AAK1", "AAAS-AAK1"],
        "DNA_SUPPORT": ["foo", 'foo', 'foo'],
        "RNA_SUPPORT": ["foo", 'foo', 'foo'],
        "METHOD": ["foo", 'foo', 'foo'],
        "FRAME": ["foo", 'foo', 'foo'],
        "ID": ['AAED', 'AAK1', 'AAAS5']})

    fusionDf = pd.DataFrame({
        "HUGO_SYMBOL": ['AAED', 'AAK1', 'AAAS5'],
        "ENTREZ_GENE_ID": [0, 0, float('nan')],
        "CENTER": ["SAGE", "SAGE", "SAGE"],
        "TUMOR_SAMPLE_BARCODE": ["ID1-1", "ID2-1", "ID1-3"],
        "FUSION": ["AAED-AAK1", "AAAS-AAK1", "AAAS-AAK1"],
        "DNA_SUPPORT": ["foo", 'foo', 'foo'],
        "RNA_SUPPORT": ["foo", 'foo', 'foo'],
        "METHOD": ["foo", 'foo', 'foo'],
        "FRAME": ["foo", 'foo', 'foo']})

    newFusionDf = fusionClass._process(fusionDf, databaseMapping)
    assert expectedFusionDf.equals(newFusionDf[expectedFusionDf.columns])


def test_validation():
    with pytest.raises(AssertionError):
        fusionClass.validateFilename(['foo'])
    assert fusionClass.validateFilename(["data_fusions_SAGE.txt"]) == "fusions"

    fusionDf = pd.DataFrame({
        "HUGO_SYMBOL": ['AAED', 'AAK1', 'AAAS'],
        "ENTREZ_GENE_ID": [0, 0, 0],
        "CENTER": ["SAGE", "SAGE", "SAGE"],
        "TUMOR_SAMPLE_BARCODE": ["ID1-1", "ID2-1", "ID1-3"],
        "FUSION": ["AAED-AAK1", "AAAS-AAK1", "AAAS-AAK1"],
        "DNA_SUPPORT": ["foo", 'foo', 'foo'],
        "RNA_SUPPORT": ["foo", 'foo', 'foo'],
        "METHOD": ["foo", 'foo', 'foo'],
        "FRAME": ["foo", 'foo', 'foo']})
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = fusionClass._validate(fusionDf, False, "syn1234")
        assert error == ""
        assert warning == ""

    fusionDf = pd.DataFrame({
        "HUGO_SYMBOL": [float('nan'), 'AAK1', 'AAAS'],
        "CENTER": ["SAGE", "SAGE", "SAGE"],
        "TUMOR_SAMPLE_BARCODE": ["ID1-1", "ID2-1", "ID1-3"],
        "FUSION": ["AAED-AAK1", "AAAS-AAK1", "AAAS-AAK1"],
        "DNA_SUPPORT": ["foo", 'foo', 'foo'],
        "RNA_SUPPORT": ["foo", 'foo', 'foo'],
        "METHOD": ["foo", 'foo', 'foo'],
        "FRAME": ["foo", 'foo', 'foo']})
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = fusionClass._validate(fusionDf, False, "syn1234")
        expectedErrors = (
            "Your fusion file must at least have these headers: ENTREZ_GENE_ID.\n"
            "Your fusion file should not have any NA/blank Hugo Symbols.\n")

        assert error == expectedErrors
        assert warning == ""
    # Do not check symbols
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = fusionClass._validate(fusionDf, True, "syn1234")
        expectedErrors = (
            "Your fusion file must at least have these headers: ENTREZ_GENE_ID.\n")
        assert error == expectedErrors
        assert warning == ""
