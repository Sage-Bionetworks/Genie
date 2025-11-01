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


def test_update_oncotree_code_tables_calls_update_with_expected_dataframes():
    """Ensure both _update_table calls receive correctly formatted data."""
    syn = mock.MagicMock()

    db_map = pd.DataFrame(
        {
            "Database": ["oncotree", "oncotreeLink", "primaryCode"],
            "Id": ["syn111", "syn222", "syn333"],
        }
    )

    # clinicaldf returned from first extract.get_syntabledf
    clinicaldf = pd.DataFrame(
        {
            "SAMPLE_ID": ["A1", "A2", "A3", "A4"],
            "CENTER": ["DFCI", "DFCI", "MSK", "MSK"],
            "ONCOTREE_CODE": ["BRCA", "LUNG", "LUNG", "SKIN"],
        }
    )

    # Mock oncotree mapping from external URL
    oncotree_mapping = {
        "BRCA": {"ONCOTREE_PRIMARY_NODE": "BREAST"},
        "LUNG": {"ONCOTREE_PRIMARY_NODE": "LUNG"},
        "SKIN": {"ONCOTREE_PRIMARY_NODE": "SKIN"},
    }

    # mock entity returned by syn.get for oncotreeLinkSynId
    oncotree_ent = mock.MagicMock()
    oncotree_ent.externalURL = "http://mock-oncotree.org"

    with (
        mock.patch.object(dash_update.extract, "get_syntabledf") as mock_extract,
        mock.patch.object(dash_update.load, "_update_table") as mock_update,
        mock.patch.object(
            dash_update.process_functions,
            "get_oncotree_code_mappings",
            return_value=oncotree_mapping,
        ),
        mock.patch.object(syn, "get", return_value=oncotree_ent),
    ):
        # Configure get_syntabledf to return clinicaldf first, then mock DB snapshots later
        mock_extract.side_effect = [
            clinicaldf,  # first call: select * from syn7517674
            pd.DataFrame(
                columns=["Oncotree_Code", "DFCI", "MSK", "Total"]
            ),  # second: oncotree DB
            pd.DataFrame(
                columns=["Oncotree_Code", "DFCI", "MSK", "Total"]
            ),  # third: primaryCode DB
        ]

        dash_update.update_oncotree_code_tables(syn, db_map)

        # Two calls to load._update_table
        assert mock_update.call_count == 2

        # First call = oncotree_code_distributiondf update
        args1, kwargs1 = mock_update.call_args_list[0]
        (
            passed_syn1,
            existing_df1,
            new_df1,
            synid1,
            key_cols1,
        ) = args1

        assert synid1 == "syn111"
        assert key_cols1 == ["Oncotree_Code"]
        assert passed_syn1 is syn

        # expected oncotree_code_distributiondf
        expected_df1 = pd.DataFrame(
            {
                "DFCI": [1, 1, 0],
                "MSK": [0, 1, 1],
                "Total": [1, 2, 1],
                "Oncotree_Code": ["BRCA", "LUNG", "SKIN"],
            },
            index=["BRCA", "LUNG", "SKIN"],
        )

        # sort by Oncotree_Code to ensure deterministic order
        assert_frame_equal(
            new_df1,
            expected_df1,
        )

        # Second call = primary_code_distributiondf update
        args2, kwargs2 = mock_update.call_args_list[1]
        (
            passed_syn2,
            existing_df2,
            new_df2,
            synid2,
            key_cols2,
        ) = args2

        assert synid2 == "syn333"
        assert key_cols2 == ["Oncotree_Code"]
        assert passed_syn2 is syn

        # expected primary_code_distributiondf
        expected_df2 = pd.DataFrame(
            {
                "DFCI": [1, 1, 0],
                "MSK": [0, 1, 1],
                "Total": [1, 2, 1],
                "Oncotree_Code": ["BREAST", "LUNG", "SKIN"],
            },
            index=["BREAST", "LUNG", "SKIN"],
        )

        assert_frame_equal(
            new_df2,
            expected_df2,
        )

        # Verify _update_table was called with to_delete=True
        for _, kwargs in mock_update.call_args_list:
            assert kwargs["to_delete"] is True
