"""Test write invalid reasons module"""
from unittest import mock
from unittest.mock import create_autospec, patch

from genie import write_invalid_reasons
import pandas as pd
import synapseclient


SYN = create_autospec(synapseclient.Synapse)
CENTER_ERRORSDF = pd.DataFrame(
    {
        "id": ["syn1234", "syn2345"],
        "errors": ["my errors|n", "errors here|f"],
        "center": ["SAGE", "TEST"],
    }
)
ENT1 = synapseclient.File(id="syn123", name="test", parentId="syn3333")


class QueryResponse:
    def asDataFrame():
        return CENTER_ERRORSDF


def test__combine_center_file_errors():
    """Test combining each center's file errors"""
    expected_error = (
        f"\t{ENT1.name} ({ENT1.id}):\n\nmy errors\nn\n\n"
        f"\t{ENT1.name} ({ENT1.id}):\n\nerrors here\nf\n\n"
    )
    calls = [
        mock.call("syn1234", downloadFile=False),
        mock.call("syn2345", downloadFile=False),
    ]
    with patch.object(SYN, "get", return_value=ENT1) as patch_synget:
        center_errors = write_invalid_reasons._combine_center_file_errors(
            SYN, CENTER_ERRORSDF
        )
        assert center_errors == expected_error
        patch_synget.assert_has_calls(calls)


def test_get_center_invalid_errors():
    """Test getting all center invalid errors"""
    with patch.object(
        SYN, "tableQuery", return_value=QueryResponse
    ) as patch_query, patch.object(
        write_invalid_reasons, "_combine_center_file_errors", return_value="errors"
    ) as patch_combine:
        center_invalid = write_invalid_reasons.get_center_invalid_errors(SYN, "syn3333")
        assert center_invalid == {"SAGE": "errors", "TEST": "errors"}
        patch_query.assert_called_once_with("SELECT * FROM syn3333")
        assert patch_combine.call_count == 2
