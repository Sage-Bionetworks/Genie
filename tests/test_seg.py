import pandas as pd
import pytest

from genie_registry.seg import seg


@pytest.fixture
def seg_class(syn):
    return seg(syn, "SAGE")


def test_processing(seg_class):
    expectedSegDf = pd.DataFrame(
        {
            "ID": [
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            "CHROM": ["1", "2", "3", "4", "5"],
            "LOCSTART": [1, 2, 3, 4, 3],
            "LOCEND": [1, 2, 3, 4, 2],
            "NUMMARK": [1, 2, 3, 4, 3],
            "SEGMEAN": [1, 2, 3.9, 4, 3],
            "CENTER": ["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
        }
    )

    segDf = pd.DataFrame(
        {
            "ID": [
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            "CHROM": ["chr1", 2, 3, 4, 5],
            "LOC.START": [1, 2, 3, 4, 3],
            "LOC.END": [1, 2, 3, 4, 2],
            "NUM.MARK": [1, 2, 3, 4, 3],
            "SEG.MEAN": [1, 2, 3.9, 4, 3],
        }
    )

    newSegDf = seg_class._process(segDf)
    assert expectedSegDf.equals(newSegDf[expectedSegDf.columns])


def test_validation_filename(seg_class):
    with pytest.raises(AssertionError):
        seg_class.validateFilename(["foo"])
    assert seg_class.validateFilename(["genie_data_cna_hg19_SAGE.seg"]) == "seg"


def test_validation_perfect(seg_class):
    segDf = pd.DataFrame(
        {
            "ID": [
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            "CHROM": [1, 2, 3, 4, 5],
            "LOC.START": [1, 2, 3, 4, 3],
            "LOC.END": [1, 2, 3, 4, 3],
            "NUM.MARK": [1, 2, 3, 4, 3],
            "SEG.MEAN": [1, 2, 3, 4, 3],
        }
    )

    error, warning = seg_class._validate(segDf)
    assert error == ""
    assert warning == ""


def test_valdation_invalid(seg_class):
    segDf = pd.DataFrame(
        {
            "ID": ["ID1", "ID2", "ID3", "ID4", "ID5"],
            "CHROM": [1, 2, float("nan"), 4, 5],
            "LOC.START": [1, 2, 3, 4, float("nan")],
            "LOC.END": [1, 2, 3, float("nan"), 3],
            "NUM.MARK": [1, 2, 3, 4, 3],
        }
    )
    expectedErrors = (
        "Your seg file is missing these headers: SEG.MEAN.\n"
        "Seg: No null or empty values allowed in column(s): "
        "CHROM, LOC.END, LOC.START.\n"
        "Seg: ID must start with GENIE-SAGE\n"
    )
    error, warning = seg_class._validate(segDf)
    assert error == expectedErrors
    assert warning == ""

    segDf = pd.DataFrame(
        {
            "ID": ["ID1", "ID2", "ID3", "ID4", "ID5"],
            "CHROM": [1, 2, 3, 4, 5],
            "LOC.START": [1, 2, 3, 4.3, 3],
            "LOC.END": [1, 2, 3.4, 4, 3],
            "NUM.MARK": [1, 2, 3, 33.3, 3],
            "SEG.MEAN": [1, 2, "f.d", 4, 3],
        }
    )
    error, warning = seg_class._validate(segDf)
    expectedErrors = (
        "Seg: Only integars allowed in these column(s): "
        "LOC.END, LOC.START, NUM.MARK.\n"
        "Seg: Only numerical values allowed in SEG.MEAN.\n"
        "Seg: ID must start with GENIE-SAGE\n"
    )
    assert error == expectedErrors
    assert warning == ""
