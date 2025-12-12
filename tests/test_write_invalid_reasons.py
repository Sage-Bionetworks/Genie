"""Test write invalid reasons module"""
import pandas as pd
import pytest
import synapseclient
from unittest import mock
from unittest.mock import patch

from genie import write_invalid_reasons

@pytest.fixture
def mock_centers_mapping():
    df = pd.DataFrame(
        {
            "center": ["A", "B"],
            "errorsSynId": ["synErrA", "synErrB"],
            "stagingSynId":["synStageA", "synStageB"],
            "inputSynId":["synInputA", "synInputB"],
        },
    )
    df.index = ["0_1", "1_1"]
    return df
    
CENTER_ERRORSDF = pd.DataFrame(
    {
        "id": ["syn1234", "syn2345"],
        "errors": ["my errors|n", "errors here|f"],
        "center": ["SAGE", "TEST"],
    }
)
ENT1 = synapseclient.File(id="syn123", name="test", parentId="syn3333")


class QueryResponse:
    @staticmethod
    def asDataFrame():
        return CENTER_ERRORSDF


def test__combine_center_file_errors(syn):
    """Test combining each center's file errors"""
    expected_error = (
        f"\t{ENT1.name} ({ENT1.id}):\n\nmy errors\nn\n\n"
        f"\t{ENT1.name} ({ENT1.id}):\n\nerrors here\nf\n\n"
    )
    calls = [
        mock.call("syn1234", downloadFile=False),
        mock.call("syn2345", downloadFile=False),
    ]
    with patch.object(syn, "get", return_value=ENT1) as patch_synget:
        center_errors = write_invalid_reasons._combine_center_file_errors(
            syn, CENTER_ERRORSDF
        )
        assert center_errors == expected_error
        patch_synget.assert_has_calls(calls)


def test_get_center_invalid_errors(syn):
    """Test getting all center invalid errors"""
    with patch.object(
        syn, "tableQuery", return_value=QueryResponse
    ) as patch_query, patch.object(
        write_invalid_reasons, "_combine_center_file_errors", return_value="errors"
    ) as patch_combine:
        center_invalid = write_invalid_reasons.get_center_invalid_errors(syn, "syn3333")
        assert center_invalid == {"SAGE": "errors", "TEST": "errors"}
        patch_query.assert_called_once_with("SELECT * FROM syn3333")
        assert patch_combine.call_count == 2


def test_write_writes_no_errors_and_correct_errors_and_uses_correct_parent_ids(mock_centers_mapping):
    syn = mock.Mock()
    center_errors = {"A": "A had errors"}  # no errors for center B

    with mock.patch.object(write_invalid_reasons.extract, "get_syntabledf", return_value=mock_centers_mapping), \
         mock.patch.object(write_invalid_reasons, "get_center_invalid_errors", return_value=center_errors), \
         mock.patch.object(write_invalid_reasons.synapseclient, "File") as m_file_cls, \
         mock.patch.object(write_invalid_reasons.os, "remove") as m_remove:

        # Make open() return a different handle per file so we can assert per-center writes
        file_handles = {}

        def _open_side_effect(filename, mode):
            fh = mock.Mock(name=f"handle_{filename}")
            ctx = mock.Mock()
            ctx.__enter__ = mock.Mock(return_value=fh)
            ctx.__exit__ = mock.Mock(return_value=False)
            file_handles[filename] = fh
            return ctx

        with mock.patch("builtins.open", side_effect=_open_side_effect) as m_open:
            # Make File() return different objects so we can assert store calls too (optional)
            ent_a = mock.Mock(name="ent_a")
            ent_b = mock.Mock(name="ent_b")
            m_file_cls.side_effect = [ent_a, ent_b]

            write_invalid_reasons.write(
                syn=syn,
                center_mapping_synid="synCenterMap",
                error_tracker_synid="synErrTrack",
            )

    # assertions: open + writes
    m_open.assert_any_call("A_validation_errors.txt", "w")
    m_open.assert_any_call("B_validation_errors.txt", "w")

    file_handles["A_validation_errors.txt"].write.assert_called_once_with("A had errors")
    file_handles["B_validation_errors.txt"].write.assert_called_once_with("No errors!")

    # assertions: correct Synapse folder IDs used (parentId)
    m_file_cls.assert_any_call("A_validation_errors.txt", parentId="synErrA")
    m_file_cls.assert_any_call("B_validation_errors.txt", parentId="synErrB")

    # Synapse store + cleanup
    assert syn.store.call_args_list == [mock.call(ent_a), mock.call(ent_b)]
    m_remove.assert_any_call("A_validation_errors.txt")
    m_remove.assert_any_call("B_validation_errors.txt")
