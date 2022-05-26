"""Test GENIE Bed class"""
import tempfile
import shutil
from unittest import mock
from unittest.mock import patch

import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
import synapseclient

import genie_registry.bed
from genie_registry.bed import bed
from genie_registry.bedSP import bedSP

if not shutil.which("bedtools"):
    pytest.skip("bedtools is not found, skipping bed tests", allow_module_level=True)

GENE_GTF_TEXT = '2\tprotein_coding\tgene\t69688532\t69901481\t.\t-\t.\tgene_id "ENSG00000115977"; gene_name "AAK1"; gene_source "ensembl_havana"; gene_biotype "protein_coding";\n9\tprotein_coding\tgene\t99401859\t99417585\t.\t-\t.\tgene_id "ENSG00000158122"; gene_name "AAED1"; gene_source "ensembl_havana"; gene_biotype "protein_coding";\n12\tprotein_coding\tgene\t53701240\t53718648\t.\t-\t.\tgene_id "ENSG00000094914"; gene_name "AAAS"; gene_source "ensembl_havana"; gene_biotype "protein_coding";\n19\tprotein_coding\tgene\t44047192\t44084625\t.\t-\t.\tgene_id "ENSG00000073050"; gene_name "XRCC1"; gene_source "ensembl_havana"; gene_biotype "protein_coding";\n19\tprotein_coding\tgene\t44080952\t44088116\t.\t+\t.\tgene_id "ENSG00000234465"; gene_name "PINLYP"; gene_source "ensembl_havana"; gene_biotype "protein_coding";'
EXON_GTF_TEXT = '2\tprocessed_transcript\texon\t69688432\t69689532\t.\t-\t.\tgene_id "ENSG00000115977"; transcript_id "ENST00000492192"; exon_number "2"; gene_name "AAK1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AAK1-009"; transcript_source "havana"; exon_id "ENSE00001882560";\n9\tprotein_coding\texon\t99416987\t99417030\t.\t-\t.\tgene_id "ENSG00000158122"; transcript_id "ENST00000411939"; exon_number "1"; gene_name "AAED1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AAED1-003"; transcript_source "havana"; exon_id "ENSE00001768346"; tag "cds_start_NF"; tag "mRNA_start_NF";\n12\tretained_intron\texon\t53702509\t53702599\t.\t-\t.\tgene_id "ENSG00000094914"; transcript_id "ENST00000550033"; exon_number "4"; gene_name "AAAS"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "AAAS-019"; transcript_source "havana"; exon_id "ENSE00003694270";\n19\tprotein_coding\texon\t44084517\t44084625\t.\t-\t.\tgene_id "ENSG00000073050"; transcript_id "ENST00000598165"; exon_number "1"; gene_name "XRCC1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "XRCC1-008"; transcript_source "havana"; exon_id "ENSE00003137784"; tag "cds_end_NF"; tag "mRNA_end_NF";\n19\tprotein_coding\texon\t44084696\t44084739\t.\t+\t.\tgene_id "ENSG00000234465"; transcript_id "ENST00000562255"; exon_number "1"; gene_name "PINLYP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PINLYP-001"; transcript_source "havana"; tag "CCDS"; ccds_id "CCDS58667"; exon_id "ENSE00002599477";'
GENE_TEMP = tempfile.NamedTemporaryFile()
EXON_TEMP = tempfile.NamedTemporaryFile()

with open(GENE_TEMP.name, "w") as gene:
    gene.write(GENE_GTF_TEXT)
with open(EXON_TEMP.name, "w") as exon:
    exon.write(EXON_GTF_TEXT)


def create_mock_table(dataframe):
    table = mock.create_autospec(synapseclient.table.CsvFileTable)
    table.asDataFrame.return_value = dataframe
    return table


def table_query_results(*args):
    return table_query_results_map[args]


symbols = pd.DataFrame(
    dict(
        hgnc_symbol=["AAK1", "AAED1", "AAAS", "PINLYP", "XRCC1"],
        chromosome_name=["2", "9", "12", "19", "19"],
        start_position=[69688532, 99401859, 53701240, 44080952, 44047192],
        end_position=[69901481, 99417585, 53718648, 44088116, 44084625],
    )
)

# This is the gene positions that all bed dataframe will be processed against
table_query_results_map = {
    ("SELECT * FROM syn11806563",): create_mock_table(symbols),
}

syn = mock.create_autospec(synapseclient.Synapse)
syn.tableQuery.side_effect = table_query_results

seq_assay_id = "SAGE-Test"
new_path = "new.bed"
parentid = "synTest"
bed_class = bed(syn, "SAGE")
bedsp_class = bedSP(syn, "SAGE")


def test_perfect___process():
    """Process perfect bed file"""
    expected_beddf = pd.DataFrame(
        dict(
            Chromosome=["2", "9", "12", "19", "19"],
            Start_Position=[69688533, 99401860, 53701241, 44084466, 44084466],
            End_Position=[69901480, 99417584, 53718647, 44084638, 44084638],
            Hugo_Symbol=["AAK1", "AAED1", "AAAS", "XRCC1", "PINLYP"],
            includeInPanel=[True, True, True, True, True],
            clinicalReported=[True, True, False, False, True],
            ID=["AAK1", "AAED1", "AAAS", "XRCC1", "foo"],
            SEQ_ASSAY_ID=[
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
                "SAGE-TEST",
            ],
            Feature_Type=["exon", "exon", "exon", "exon", "exon"],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
        )
    )

    expected_beddf.sort_values("ID", inplace=True)
    expected_beddf.reset_index(drop=True, inplace=True)
    beddf = pd.DataFrame(
        {
            0: ["2", "9", "12", "19", "19"],
            1: [69688533, 99401860, 53701241, 44084466, 44084466],
            2: [69901480, 99417584, 53718647, 44084638, 44084638],
            3: ["AAK1", "AAED1", "AAAS", "XRCC1", "foo"],
            4: [True, True, True, 1, 1],
            5: [True, True, False, 0, 1],
        }
    )
    with patch.object(
        genie_registry.bed, "create_gtf", return_value=(EXON_TEMP.name, GENE_TEMP.name)
    ):
        new_beddf = bed_class._process(
            beddf, seq_assay_id, new_path, parentid, create_panel=False
        )
        new_beddf.sort_values("ID", inplace=True)
        new_beddf.reset_index(drop=True, inplace=True)
        assert_frame_equal(
            expected_beddf, new_beddf[expected_beddf.columns], check_dtype=False
        )


def test_includeinpanel___process():
    """
    Make sure includeInPanel column is propogated and
    intergenic region is captured
    """

    expected_beddf = pd.DataFrame(
        dict(
            Chromosome=["2", "9", "12", "19"],
            Start_Position=[69688432, 1111, 53700240, 44080953],
            End_Position=[69689532, 1111, 53719548, 44084624],
            Hugo_Symbol=["AAK1", float("nan"), "AAAS", float("nan")],
            includeInPanel=[True, True, False, True],
            clinicalReported=[float("nan")] * 4,
            ID=["foo", "bar", "baz", "boo"],
            SEQ_ASSAY_ID=["SAGE-TEST", "SAGE-TEST", "SAGE-TEST", "SAGE-TEST"],
            Feature_Type=["exon", "intergenic", "exon", "exon"],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE"],
        )
    )

    expected_beddf.sort_values("Chromosome", inplace=True)
    expected_beddf.reset_index(drop=True, inplace=True)

    # symbols that can't be map should be null,
    # includeInPanel column should be included if it exists
    beddf = pd.DataFrame(
        {
            0: ["2", "9", "12", "19"],
            1: [69688432, 1111, 53700240, 44080953],
            2: [69689532, 1111, 53719548, 44084624],
            3: ["foo", "bar", "baz", "boo"],
            4: [True, True, 0, 1],
        }
    )
    with patch.object(
        genie_registry.bed, "create_gtf", return_value=(EXON_TEMP.name, GENE_TEMP.name)
    ):
        new_beddf = bedsp_class._process(
            beddf, seq_assay_id, new_path, parentid, create_panel=False
        )
        new_beddf.sort_values("Chromosome", inplace=True)
        new_beddf.reset_index(drop=True, inplace=True)
        assert_frame_equal(
            expected_beddf, new_beddf[expected_beddf.columns], check_dtype=False
        )


def test_clinicalreport___process():
    expected_beddf = pd.DataFrame(
        dict(
            Chromosome=["2", "9", "12", "19"],
            Start_Position=[69688432, 1111, 53700240, 44080953],
            End_Position=[69689532, 1111, 53719548, 44084624],
            Hugo_Symbol=["AAK1", float("nan"), "AAAS", float("nan")],
            includeInPanel=[True, True, False, True],
            clinicalReported=[float("nan")] * 4,
            ID=["foo", "bar", "baz", "boo"],
            SEQ_ASSAY_ID=["SAGE-TEST", "SAGE-TEST", "SAGE-TEST", "SAGE-TEST"],
            Feature_Type=["exon", "intergenic", "exon", "exon"],
            CENTER=["SAGE", "SAGE", "SAGE", "SAGE"],
        )
    )

    expected_beddf.sort_values("Chromosome", inplace=True)
    expected_beddf.reset_index(drop=True, inplace=True)

    # symbols that can't be map should be null,
    # includeInPanel column should be included if it exists
    beddf = pd.DataFrame(
        {
            0: ["2", "9", "12", "19"],
            1: [69688432, 1111, 53700240, 44080953],
            2: [69689532, 1111, 53719548, 44084624],
            3: ["foo", "bar", "baz", "boo"],
            4: [True, True, False, True],
            5: [True, float("nan"), False, True],
        }
    )
    with patch.object(
        genie_registry.bed, "create_gtf", return_value=(EXON_TEMP.name, GENE_TEMP.name)
    ):
        new_beddf = bedsp_class._process(
            beddf, seq_assay_id, new_path, parentid, create_panel=False
        )
        new_beddf.sort_values("Chromosome", inplace=True)
        new_beddf.reset_index(drop=True, inplace=True)
        assert_frame_equal(
            expected_beddf, new_beddf[expected_beddf.columns], check_dtype=False
        )


def test_filetype():
    assert bed_class._fileType == "bed"
    assert bedsp_class._fileType == "bedSP"


@pytest.fixture(
    params=[(["foo"]), (["SAGE-test.txt"]), (["foo"]), (["nonGENIE_SAGE-test.txt"])]
)
def filename_fileformat_map(request):
    return request.param


def test_incorrect_validatefilename(filename_fileformat_map):
    filepath_list = filename_fileformat_map
    with pytest.raises(AssertionError):
        bed_class.validateFilename(filepath_list)


def test_correct_validatefilename():
    assert bed_class.validateFilename(["SAGE-test.bed"]) == "bed"
    assert bedsp_class.validateFilename(["nonGENIE_SAGE-test.bed"]) == "bedSP"


def test_perfect__validate():
    bedDf = pd.DataFrame(
        dict(
            a=["2", "9", "12"],
            b=[69688533, 99401860, 53701241],
            c=[69901480, 99417584, 53718647],
            d=["AAK1", "AAED1", "AAAS"],
            e=[True, True, False],
            f=[True, True, False],
        )
    )

    error, warning = bed_class._validate(bedDf)
    assert error == ""
    assert warning == ""


def test_90percent_boundary__validate():
    bedDf = pd.DataFrame(
        dict(
            a=["2", "9", "12"],
            b=[69688432, 99416585, 53700240],
            c=[69689532, 99417685, 53719548],
            d=["AAK1", "AAED1", "AAAS"],
            e=[True, True, False],
            f=[True, True, False],
        )
    )
    error, warning = bed_class._validate(bedDf)
    assert error == ""
    assert warning == ""


def test_missingcols_failure__validate():
    emptydf = pd.DataFrame()
    error, warning = bed_class._validate(emptydf)
    expected_errors = (
        "BED file: Must at least have five columns in this order: "
        "Chromosome, Start_Position, End_Position, Hugo_Symbol, "
        "includeInPanel. Make sure there are no headers.\n"
    )
    assert error == expected_errors
    assert warning == ""


def test_hugosymbol_failure__validate():
    bedDf = pd.DataFrame(
        dict(
            a=["2", "9", "12"],
            b=[69688533, 99401860, 53701241],
            c=[69901480, 99417584, 53718647],
            d=["+", float("nan"), "AAAS"],
            e=[True, True, False],
        )
    )
    error, warning = bed_class._validate(bedDf)
    expected_errors = (
        "BED file: You cannot submit any null symbols.\n"
        "BED file: Fourth column must be the Hugo_Symbol column, "
        "not the strand column\n"
    )
    assert error == expected_errors
    assert warning == ""


def test_badinputs_failure__validate():
    bedDf = pd.DataFrame(
        dict(
            a=["2", "9", "12"],
            b=["69688533", 99401860, 53701241],
            c=[69901480, "99417584", 53718647],
            d=["AAK1", "AAED1", "AAAS"],
            e=[True, False, "foobar"],
        )
    )

    error, warning = bed_class._validate(bedDf)
    expected_errors = (
        "BED file: The Start_Position column must only be integers. "
        "Make sure there are no headers.\n"
        "BED file: The End_Position column must only be integers. "
        "Make sure there are no headers.\n"
        "BED file: Please double check your includeInPanel column.  "
        "This column must only be these values: True, False\n"
    )
    assert error == expected_errors
    assert warning == ""


def test_90percentboundary_failure__validate():
    # Test 90% boundary failure boundary, with incorrect gene names
    bedDf = pd.DataFrame(
        dict(
            a=["2", "9", "12"],
            b=[69901381, 4345, 11111],
            c=[69911481, 99417590, 11113],
            d=["foo", "foo", "AAAS"],
            e=[True, True, False],
        )
    )

    error, warning = bedsp_class._validate(bedDf)
    expected_errors = (
        "BED file: You have no correct gene symbols. "
        "Make sure your gene symbol column (4th column) is formatted like so: "
        "SYMBOL(;optionaltext).  Optional text can be semi-colon separated.\n"
    )
    expected_warnings = (
        "BED file: Any gene names that can't be remapped will be null.\n"
    )
    assert error == expected_errors
    assert warning == expected_warnings


def test_overlapping__validate():
    """
    Check if genes that overlap have no errors
    - Submitted bed start position in a gene
    - Submitted bed end position in a gene
    - Submitted bed region surrounds a gene
    """
    beddf = pd.DataFrame(
        dict(
            a=["2", "9", "2"],
            b=[1111, 4345, 69901281],
            c=[69880186, 99417590, 70000000],
            d=["AAK1", "AAED1", "AAK1"],
            e=[True, False, True],
        )
    )
    error, warning = bedsp_class._validate(beddf)
    assert error == ""
    assert warning == ""


def test_symbolnull_failure__validate():
    # Test 2 gene symbols returned NULL
    bedDf = pd.DataFrame(
        dict(a=["19"], b=[44080953], c=[44084624], d=["AAK1"], e=[False])
    )

    error, warning = bedsp_class._validate(bedDf)
    expected_errors = (
        "BED file: You have no correct gene symbols. "
        "Make sure your gene symbol column (4th column) is formatted like so: "
        "SYMBOL(;optionaltext).  Optional text can be semi-colon separated.\n"
    )
    expected_warnings = (
        "BED file: Any gene names that can't be remapped will be null.\n"
    )
    assert error == expected_errors
    assert warning == expected_warnings
