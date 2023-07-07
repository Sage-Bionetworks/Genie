"""Test genie.transform module"""
from io import BytesIO
from unittest.mock import patch

import pandas as pd
from pandas.api.types import (
    is_object_dtype,
)
import pytest

from genie import transform
from unittest import mock


class TestConvertCols:
    @pytest.mark.parametrize(
        "test_input, expected",
        [
            (pd.DataFrame({"some_col": [10.0, float("nan")]}), ["10.0", float("nan")]),
            (pd.DataFrame({"some_col": [1, None]}), ["1.0", None]),
            (
                pd.DataFrame({"some_col": ["Val1", float("nan")]}),
                ["Val1", float("nan")],
            ),
            (
                pd.DataFrame({"some_col": ["Val1", float("nan"), 1, "1"]}),
                ["Val1", float("nan"), "1", "1"],
            ),
        ],
        ids=["float_w_na", "int_w_na", "string_w_na", "mixed_w_na"],
    )
    def test_that__convert_col_with_nas_to_str_keep_na_for_any_data_type(
        self, test_input, expected
    ):
        result = transform._convert_col_with_nas_to_str(test_input, "some_col")
        assert result[0] == expected[0]
        assert pd.isna(result[1])

    @pytest.mark.parametrize(
        "test_input, expected",
        [
            (pd.DataFrame({"some_col": [10.0, 11.2]}), ["10.0", "11.2"]),
            (
                pd.DataFrame({"some_col": ["Val1", "Val2"]}),
                ["Val1", "Val2"],
            ),
            (
                pd.DataFrame({"some_col": ["Val1", 1, "2"]}),
                ["Val1", "1", "2"],
            ),
        ],
        ids=["float_no_na", "string_no_na", "mixed_no_na"],
    )
    def test_that__convert_col_with_nas_to_str_returns_correct_vals_with_no_na_data(
        self, test_input, expected
    ):
        result = transform._convert_col_with_nas_to_str(test_input, "some_col")
        assert result == expected

    def test_that__convert_float_col_with_nas_to_int(self):
        test_input = pd.DataFrame({"some_col": [10.0, float("nan")]})
        result = transform._convert_float_col_with_nas_to_int(test_input, "some_col")
        assert result[0] == 10
        assert pd.isna(result[1])

    @pytest.mark.parametrize(
        "test_input, expected",
        [
            (pd.DataFrame({"some_col": [10.0, 11.2]}), [10.0, 11.2]),
            (
                pd.DataFrame({"some_col": ["Val1", "Val2"]}),
                ["Val1", "Val2"],
            ),
            (pd.DataFrame({"some_col": [10, 11]}), [10, 11]),
        ],
        ids=["float_no_na", "string_no_na", "int_no_na"],
    )
    def test_that__convert_float_col_with_nas_to_int_does_nothing_with_no_na_data(
        self, test_input, expected
    ):
        result = transform._convert_float_col_with_nas_to_int(test_input, "some_col")
        assert result == expected

    def test_that__convert_float_col_with_nas_to_int_does_nothing_with_str_data(self):
        test_input = pd.DataFrame({"some_col": ["Val1", float("nan")]})
        result = transform._convert_float_col_with_nas_to_int(test_input, "some_col")
        assert result[0] == "Val1"
        assert pd.isna(result[1])


def read_csv_side_effect(**kwargs):
    if kwargs["low_memory"] == True:
        raise pd.errors.DtypeWarning("some_warning")


class TestConvertMixedDtypes:
    @pytest.fixture(scope="session")
    def test_input(self):
        yield pd.DataFrame({"some_col": [1, "Val2", "1"]})

    @pytest.fixture(scope="session")
    def bytes_io_output(self, test_input):
        output = BytesIO()
        test_input.to_csv(output, index=False)
        output.seek(0)
        yield output

    def test_that__convert_df_with_mixed_dtypes_gets_expected_output(self, test_input):
        # Create your in memory BytesIO file.
        output = BytesIO()
        test_input.to_csv(output, index=False)
        output.seek(0)  # Contains the CSV in memory file.

        df = transform._convert_df_with_mixed_dtypes(
            {"filepath_or_buffer": output, "index_col": False}
        )
        pd.testing.assert_frame_equal(
            df.reset_index(drop=True), pd.DataFrame({"some_col": ["1", "Val2", "1"]})
        )
        assert is_object_dtype(df["some_col"])

    def test_that__convert_df_with_mixed_dtypes_calls_read_csv_once_if_no_exception(
        self, bytes_io_output
    ):
        with mock.patch.object(pd, "read_csv") as mock_read_csv:
            transform._convert_df_with_mixed_dtypes(
                {"filepath_or_buffer": bytes_io_output, "index_col": False}
            )
            mock_read_csv.assert_called_once_with(
                filepath_or_buffer=bytes_io_output,
                index_col=False,
                low_memory=True,
            )

    def test_that__convert_df_with_mixed_dtypes_catches_mixed_dtype_exception(
        self, bytes_io_output
    ):
        # Create your in memory BytesIO file.
        with mock.patch.object(
            pd, "read_csv", side_effect=read_csv_side_effect
        ) as mock_read_csv:
            transform._convert_df_with_mixed_dtypes(
                {"filepath_or_buffer": bytes_io_output, "index_col": False}
            )
            mock_read_csv.assert_has_calls(
                [
                    mock.call(
                        filepath_or_buffer=bytes_io_output,
                        index_col=False,
                        low_memory=True,
                    ),
                    mock.call(
                        filepath_or_buffer=bytes_io_output,
                        index_col=False,
                        low_memory=False,
                        engine="c",
                    ),
                ],
                any_order=False,
            )
