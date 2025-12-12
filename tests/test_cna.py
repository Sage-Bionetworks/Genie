from types import SimpleNamespace
from unittest import mock
from unittest.mock import patch

import pytest
import pandas as pd
import synapseclient

from genie_registry.cna import cna


def createMockTable(dataframe):
    table = mock.create_autospec(synapseclient.table.CsvFileTable)
    table.asDataFrame.return_value = dataframe
    return table


def table_query_results(*args):
    return table_query_results_map[args]


database_mapping = pd.DataFrame(dict(Database=["bed"], Id=["syn8457748"]))
symbols = pd.DataFrame(
    dict(
        Hugo_Symbol=["AAK1", "AAED1", "AAAS", "AAED1"],
        ID=["AAK1", "AAED", "AAAS", "AAD"],
    )
)
# This is the gene positions that all bed dataframe will be processed against
table_query_results_map = {
    ("select Hugo_Symbol, ID from syn11600834 where CENTER = 'SAGE'",): createMockTable(
        symbols
    ),
}

ENTITY = synapseclient.Project("testing", annotations={"dbMapping": ["syn10967259"]})


@pytest.fixture
def cna_class(syn, genie_config):
    syn.tableQuery.side_effect = table_query_results
    return cna(syn, "SAGE", genie_config)


def test_processing(cna_class):
    order = ["Hugo_Symbol", "Entrez_gene_id", "GENIE-SAGE-Id1-1", "GENIE-SAGE-Id2-1"]

    expected_cnadf = pd.DataFrame(
        {
            "Hugo_Symbol": ["AAED1", "AAK1", "AAAS"],
            "GENIE-SAGE-Id1-1": [-0.5, 2.0, 0.5],
            "GENIE-SAGE-Id2-1": [1.0, 1.5, -1.5],
        }
    )

    cnadf = pd.DataFrame(
        {
            "Hugo_Symbol": ["AAED", "AAK1", "AAAS"],
            "Entrez_gene_id": [0, 0, 0],
            "GENIE-SAGE-Id1-1": [-0.5, 2, 0.5],
            "GENIE-SAGE-Id2-1": [1, 1.5, -1.5],
        }
    )
    cnadf = cnadf[order]
    new_cnadf = cna_class._process(cnadf)
    assert expected_cnadf.equals(new_cnadf[expected_cnadf.columns])

    order = ["Hugo_Symbol", "Entrez_gene_id", "GENIE-SAGE-Id1-1", "GENIE-SAGE-Id2-1"]

    expectedCnaDf = pd.DataFrame(
        {
            "Hugo_Symbol": ["AAAS", "AAED1"],
            "GENIE-SAGE-Id1-1": [float("nan"), 1.0],
            "GENIE-SAGE-Id2-1": [float("nan"), 2.0],
        }
    )

    cnaDf = pd.DataFrame(
        {
            "Hugo_Symbol": ["AAED", "AAED1", "foo", "AAAS", "AAD"],
            "Entrez_gene_id": [0, 0, 0, 0, 0],
            "GENIE-SAGE-Id1-1": [1, 1, 0, float("nan"), 1],
            "GENIE-SAGE-Id2-1": [2, 0, -1, float("nan"), float("nan")],
        }
    )
    cnaDf = cnaDf[order]
    newCnaDf = cna_class._process(cnaDf)
    newCnaDf.reset_index(inplace=True, drop=True)
    pd.testing.assert_frame_equal(expectedCnaDf, newCnaDf[expectedCnaDf.columns])


def test_validation(syn, cna_class):
    with pytest.raises(AssertionError):
        cna_class.validateFilename(["foo"])
    assert cna_class.validateFilename(["data_CNA_SAGE.txt"]) == "cna"

    order = ["Hugo_Symbol", "Entrez_gene_id", "GENIE-SAGE-ID1-1", "GENIE-SAGE-ID2-1"]
    cnaDf = pd.DataFrame(
        {
            "Hugo_Symbol": ["AAED", "AAK1", "AAAS"],
            "Entrez_gene_id": [0, 0, 0],
            "GENIE-SAGE-ID1-1": [-2, 0, 1],
            "GENIE-SAGE-ID2-1": [0.5, 1.5, -1.5],
            "GENIE-SAGE-ID3-1": [float("nan"), 2, "NA"],
        }
    )
    cnaDf = cnaDf[order]
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = cna_class._validate(cnaDf, False)
        assert error == ""
        assert warning == ""

    cnaDf = pd.DataFrame(
        {
            "Hugo_Symbol": ["foo", "AAED", "AAED1", "AAD"],
            "GENIE-SAGE-ID1-1": [1, 2, float("nan"), 1],
            "GENIE-SAGE-ID2-1": [2, 1, -1, 2],
        }
    )
    cnaDf.sort_values("Hugo_Symbol", inplace=True)
    cnaDf = cnaDf[["GENIE-SAGE-ID1-1", "Hugo_Symbol", "GENIE-SAGE-ID2-1"]]

    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = cna_class._validate(cnaDf, False)
        expectedErrors = (
            "Your cnv file's first column must be Hugo_Symbol\n"
            "Your CNA file has duplicated Hugo_Symbols (After "
            "remapping of genes): AAD,AAED,AAED1 -> AAED1,AAED1,AAED1.\n"
        )

    assert error == expectedErrors
    assert warning == ""
    order = ["Hugo_Symbol", "GENIE-SAGE-ID1-1", "GENIE-SAGE-ID2-1"]

    cnaDf = pd.DataFrame(
        {
            "Hugo_Symbol": ["AAK1", "AAAS", "AAED1"],
            "GENIE-SAGE-ID1-1": [1, 2, 1],
            "GENIE-SAGE-ID2-1": [2, 1, 3],
        }
    )
    cnaDf = cnaDf[order]
    with patch.object(syn, "get", return_value=ENTITY):
        error, warning = cna_class._validate(cnaDf, False)
        expectedErrors = (
            "All values must be NA/blank, -2, -1.5, -1, -0.5, 0, "
            "0.5, 1, 1.5, or 2.\n"
        )
        assert error == expectedErrors
        assert warning == ""


@pytest.mark.parametrize(
    "skip_database_checks, cnv_hugo_symbols, remap_dict, expected_error, expect_tablequery_called",
    [
        # 1) skip_database_checks=False, no duplicates after remap -> ""
        (
            False,
            ["TP53", "EGFR", "BRCA1"],
            {"TP53": "TP53", "EGFR": "EGFR", "BRCA1": "BRCA1"},
            "",
            True,
        ),
        # 2) skip_database_checks=False, duplicates after remap -> error
        (
            False,
            ["TP53", "TP53_ALT", "EGFR"],
            {"TP53": "TP53", "TP53_ALT": "TP53", "EGFR": "EGFR"},
            "Your CNA file has duplicated Hugo_Symbols (After remapping of genes): "
            "TP53,TP53_ALT -> TP53,TP53.\n",
            True,
        ),
        # 3) skip_database_checks=True, duplicates exist but check skipped -> ""
        (
            True,
            ["TP53", "TP53_ALT", "EGFR"],
            {"TP53": "TP53", "TP53_ALT": "TP53", "EGFR": "EGFR"},
            "",
            False,
        ),
    ],
    ids=["no_dups", "dups", "dups_but_skip_database_checks"],
)
def test_validate_no_dup_symbols_after_remapping(
    cna_class,
    skip_database_checks,
    cnv_hugo_symbols,
    remap_dict,
    expected_error,
    expect_tablequery_called,
):
    cnvDF = pd.DataFrame({"HUGO_SYMBOL": cnv_hugo_symbols})

    # mock bed table returned by Synapse tableQuery
    bedDf = pd.DataFrame(
        {"Hugo_Symbol": list(remap_dict.values()), "ID": range(len(remap_dict))}
    )
    mock_bed_query_result = SimpleNamespace(asDataFrame=lambda: bedDf)

    def mock_validateSymbol(sym, bedDf):
        return remap_dict.get(sym)

    with patch.object(
        cna_class.syn, "tableQuery", return_value=mock_bed_query_result
    ) as mocked_tq, patch(
        "genie_registry.cna.validateSymbol", side_effect=mock_validateSymbol
    ) as mocked_vs:
        err = cna_class.validate_no_dup_symbols_after_remapping(
            cnvDF, skip_database_checks
        )

    assert err == expected_error

    if expect_tablequery_called:
        mocked_tq.assert_called_once()
        mocked_vs.assert_called()
    else:
        mocked_tq.assert_not_called()
        mocked_vs.assert_not_called()
