import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from unittest import mock

from genie import dashboard_table_updater as dash_update


@pytest.fixture
def mock_syn(tmp_path):
    """Fixture for a mock synapse client and test file."""
    syn = mock.MagicMock()

    # create a fake clinical file
    clinical_path = tmp_path / "clinical.txt"
    clinical_path.write_text("SAMPLE_ID\nS1\nS2\nS3\n")

    # entity returned by syn.get for the clinical file
    clinical_entity = mock.MagicMock()
    clinical_entity.path = str(clinical_path)

    # syn.get() will return this entity
    syn.get.return_value = clinical_entity

    # Default columns of existing table
    syn.getTableColumns.return_value = [{"name": "SAMPLE_ID"}]

    return syn


@pytest.mark.parametrize(
    "input_string_time, expected_output_time",
    [
        ("2018-10-25T20:16:07.959Z", 1540498567000),
        ("2018-10-25T20:16:07.959", 1540498567000),
        ("2018-04-06T18:30:00", 1523039400000),
    ],
    ids=["utc_time", "local_time_zone", "no_milliseconds"],
)
def test_that_string_to_unix_epoch_time_milliseconds_gives_expected_time(
    input_string_time, expected_output_time
):
    output = dash_update.string_to_unix_epoch_time_milliseconds(input_string_time)
    assert output == expected_output_time


def test_that_update_samples_in_release_table_adds_column_and_calls_update(
    tmp_path, mock_syn
):
    """Test that a new release column is added and load._update_table gets correct data."""
    file_mapping = {"sample": "syn123"}
    release = "5.3-consortium"
    samples_in_release_synid = "syn999"

    clinical_df = pd.DataFrame({"SAMPLE_ID": ["S1", "S2", "S3"]})
    existing_df = pd.DataFrame({"SAMPLE_ID": ["S1"], release: [1]})

    with (
        mock.patch.object(pd, "read_csv", return_value=clinical_df) as mock_read,
        mock.patch.object(
            dash_update.extract, "get_syntabledf", return_value=existing_df
        ) as mock_extract,
        mock.patch.object(dash_update.load, "_update_table") as mock_update,
    ):
        dash_update.update_samples_in_release_table(
            syn=mock_syn,
            file_mapping=file_mapping,
            release=release,
            samples_in_release_synid=samples_in_release_synid,
        )

        # assertions on Synapse calls
        mock_syn.get.assert_has_calls(
            [
                mock.call(file_mapping["sample"], followLink=True),
                mock.call(samples_in_release_synid),
            ]
        )
        mock_syn.getTableColumns.assert_called_once_with(samples_in_release_synid)

        # assertions on pd.read_csv
        mock_read.assert_called_once()
        assert list(mock_read.return_value.columns) == ["SAMPLE_ID"]

        # assertions on extract.get_syntabledf
        mock_extract.assert_called_once_with(
            syn=mock_syn,
            query_string=f'SELECT SAMPLE_ID, "{release}" FROM {samples_in_release_synid}',
        )

        # assertions on load._update_table inputs
        mock_update.assert_called_once()
        args, kwargs = mock_update.call_args

        # extract arguments
        _, _, samples_in_releasedf, synid_arg, key_cols = args

        # Ensure correct Synapse ID and key columns
        assert synid_arg == samples_in_release_synid
        assert key_cols == ["SAMPLE_ID"]

        # check that new samples and old samples are in expected order
        assert_frame_equal(
            samples_in_releasedf.reset_index(drop=True),
            pd.DataFrame(
                {"SAMPLE_ID": ["S2", "S3", "S1"], "5.3-consortium": [1, 1, 1]}
            ),
        )


def test_that_update_samples_in_release_table_existing_column_calls_update_directly(
    mock_syn,
):
    """If the release column already exists, we skip creating new column but still call _update_table."""
    release = "5.3-consortium"
    file_mapping = {"sample": "syn123"}
    samples_in_release_synid = "syn999"

    # pre-xisting release column
    mock_syn.getTableColumns.return_value = [{"name": "SAMPLE_ID"}, {"name": release}]

    clinical_df = pd.DataFrame({"SAMPLE_ID": ["S1", "S2"]})
    existing_df = pd.DataFrame({"SAMPLE_ID": ["S2", "S1"], release: [1, 1]})

    with (
        mock.patch.object(pd, "read_csv", return_value=clinical_df) as mock_read,
        mock.patch.object(
            dash_update.extract, "get_syntabledf", return_value=existing_df
        ) as mock_extract,
        mock.patch.object(dash_update.load, "_update_table") as mock_update,
    ):
        dash_update.update_samples_in_release_table(
            syn=mock_syn,
            file_mapping=file_mapping,
            release=release,
            samples_in_release_synid=samples_in_release_synid,
        )
        # assert that samples_in_release schema is not retrieved by syn.get
        mock_syn.get.assert_has_calls(
            [
                mock.call(file_mapping["sample"], followLink=True),
            ]
        )
        args, kwargs = mock_update.call_args

        # extract arguments
        _, _, samples_in_releasedf, synid_arg, key_cols = args

        # ensure correct Synapse ID and key columns
        assert synid_arg == samples_in_release_synid
        assert key_cols == ["SAMPLE_ID"]

        # check that new samples and old samples are in expected order
        assert_frame_equal(
            samples_in_releasedf.reset_index(drop=True),
            pd.DataFrame({"SAMPLE_ID": ["S1", "S2"], "5.3-consortium": [1, 1]}),
        )
