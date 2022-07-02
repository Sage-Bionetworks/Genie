from unittest import mock

import pandas as pd
import pytest
import synapseclient

from genie_registry.sampleRetraction import sampleRetraction
from genie_registry.patientRetraction import patientRetraction

syn = mock.create_autospec(synapseclient.Synapse)

sr = sampleRetraction(syn, "SAGE")
pr = patientRetraction(syn, "SAGE")


def test_processing():

    expectedsrDf = pd.DataFrame(
        dict(
            genieSampleId=[
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ],
            retractionDate=[
                1523039400000,
                1523039400000,
                1523039400000,
                1523039400000,
                1523039400000,
            ],
            center=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
        )
    )

    srDf = pd.DataFrame(
        {
            0: [
                "GENIE-SAGE-ID1-1",
                "GENIE-SAGE-ID2-1",
                "GENIE-SAGE-ID3-1",
                "GENIE-SAGE-ID4-1",
                "GENIE-SAGE-ID5-1",
            ]
        }
    )

    newsrDf = sr._process(srDf, "2018-04-06T18:30:00")
    assert expectedsrDf.equals(newsrDf[expectedsrDf.columns])

    expectedprDf = pd.DataFrame(
        dict(
            geniePatientId=[
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ],
            retractionDate=[
                1523125800000,
                1523125800000,
                1523125800000,
                1523125800000,
                1523125800000,
            ],
            center=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
        )
    )

    prDf = pd.DataFrame(
        {
            0: [
                "GENIE-SAGE-ID1",
                "GENIE-SAGE-ID2",
                "GENIE-SAGE-ID3",
                "GENIE-SAGE-ID4",
                "GENIE-SAGE-ID5",
            ]
        }
    )

    newprDf = pr._process(prDf, "2018-04-07T18:30:00")
    assert expectedprDf.equals(newprDf[expectedprDf.columns])


def test_validation():
    with pytest.raises(AssertionError):
        sr.validateFilename(["foo"])
    assert sr.validateFilename(["sampleRetraction.csv"]) == "sampleRetraction"

    with pytest.raises(AssertionError):
        pr.validateFilename(["foo"])
    assert pr.validateFilename(["patientRetraction.csv"]) == "patientRetraction"
