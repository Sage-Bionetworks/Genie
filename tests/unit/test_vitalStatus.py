import synapseclient
import pandas as pd
import mock
from genie.vitalStatus import vitalStatus
import pytest

syn = mock.create_autospec(synapseclient.Synapse)

vs = vitalStatus(syn, "SAGE")


def test_processing():

    expectedvsDf = pd.DataFrame(dict(
        PATIENT_ID=["GENIE-SAGE-ID1", "GENIE-SAGE-ID2", "GENIE-SAGE-ID3",
                    "GENIE-SAGE-ID4", "GENIE-SAGE-ID5"],
        YEAR_DEATH=[1999, 2000, 3000, 1000, 1999],
        YEAR_CONTACT=[1999, 2000, 3000, 1000, 1999],
        INT_CONTACT=[1, 2, 3, 4, 3],
        INT_DOD=[1, 2, 3, 4, 3],
        DEAD=[True, False, True, False, True],
        CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"]))

    vsDf = pd.DataFrame(dict(
        PATIENT_ID=["ID1", "ID2", "ID3", "ID4", "ID5"],
        YEAR_DEATH=[1999, 2000, 3000, 1000, 1999],
        YEAR_CONTACT=[1999, 2000, 3000, 1000, 1999],
        INT_CONTACT=[1, 2, 3, 4, 3],
        INT_DOD=[1, 2, 3, 4, 3],
        DEAD=[True, False, True, False, True]))

    newvsDf = vs._process(vsDf)
    assert expectedvsDf.equals(newvsDf[expectedvsDf.columns])


def test_validation():
    with pytest.raises(AssertionError):
        vs.validateFilename(["foo"])

    assert vs.validateFilename(["vital_status.txt"]) == "vitalStatus"

    vsDf = pd.DataFrame(dict(
        PATIENT_ID=["ID1", "ID2", "ID3", "ID4", "ID5"],
        YEAR_DEATH=[1999, 2000, 3000, 1000, float('nan')],
        YEAR_CONTACT=[1999, 2000, float('nan'), 2342, 1999],
        INT_CONTACT=[1, 2, ">32485", float('nan'), "<6570"],
        INT_DOD=[1, ">32485", 3, "<6570", float('nan')],
        DEAD=[True, False, True, False, float('nan')]))

    error, warning = vs._validate(vsDf)
    assert error == ""
    assert warning == ""

    vsDf = pd.DataFrame()
    error, warning = vs._validate(vsDf)
    expectedErrors = ("Vital status file: Must have PATIENT_ID column.\n"
                      "Vital status file: Must have YEAR_DEATH column.\n"
                      "Vital status file: Must have YEAR_CONTACT column.\n"
                      "Vital status file: Must have INT_CONTACT column.\n"
                      "Vital status file: Must have INT_DOD column.\n"
                      "Vital status file: Must have DEAD column.\n")
    assert error == expectedErrors
    assert warning == ""

    vsDf = pd.DataFrame(dict(
        PATIENT_ID=["ID1", "ID2", "ID3", "ID4", "ID5"],
        YEAR_DEATH=[1999, 2000, 3000, 3, 'asdfs'],
        YEAR_CONTACT=[1999, 2000, 3000, 2, 'assdfd'],
        INT_CONTACT=[1, 2, 3, 4, 'string'],
        INT_DOD=[1, 2, 3, 4, 'string'],
        DEAD=[True, False, True, 'boobar', True]))

    error, warning = vs._validate(vsDf)
    expectedErrors = (
        "Vital status file: Please double check your YEAR_DEATH column, "
        "it must be an integer in YYYY format or an empty string.\n"
        "Vital status file: Please double check your YEAR_CONTACT column, "
        "it must be an integer in YYYY format or an empty string.\n"
        "Vital status file: Please double check your INT_CONTACT column, "
        "it must be an integer, an empty string, >32485, or <6570.\n"
        "Vital status file: Please double check your INT_DOD column, "
        "it must be an integer, an empty string, >32485, or <6570.\n"
        "Vital status file: Please double check your DEAD column, "
        "it must be a boolean value or an empty string.\n")
    assert error == expectedErrors
    assert warning == ""
