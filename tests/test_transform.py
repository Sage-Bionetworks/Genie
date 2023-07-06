"""Test genie.transform module"""
from io import BytesIO
from unittest.mock import patch

import pandas as pd
from pandas.api.types import (
    is_float_dtype,
    is_integer_dtype,
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


@pytest.mark.parametrize(
    "test_input,expected_output,expected_dtype", 
    [
        (
            pd.DataFrame({"some_col": [1, "Val2", "1"]}),
            pd.DataFrame({"some_col": ["1", "Val2", "1"]}),
            is_object_dtype,
        ),
        (
            pd.DataFrame({"some_col": ["Val1", "Val2", "Val3"]}),
            pd.DataFrame({"some_col": ["Val1", "Val2", "Val3"]}),
            is_object_dtype,
        ),
        (
            pd.DataFrame({"some_col": [1.0, "Val2", 3]}),
            pd.DataFrame({"some_col": ["1.0", "Val2", "3"]}),
            is_object_dtype,
        ),
        (
            pd.DataFrame({"some_col": ["2", "3", 1]}),
            pd.DataFrame({"some_col": [2, 3, 1]}),
            is_integer_dtype,
        ),
        (
            pd.DataFrame({"some_col": [2, None, 1]}),
            pd.DataFrame({"some_col": [2, None, 1]}),
            is_float_dtype,
        ),
        (
            pd.DataFrame({"some_col": [float("nan"), None, "1"]}),
            pd.DataFrame({"some_col": [float("nan"), None, 1]}),
            is_float_dtype,
        ),
        (
            pd.DataFrame({"some_col": [float("nan"), "Val1", "1"]}),
            pd.DataFrame({"some_col": [float("nan"), "Val1", "1"]}),
            is_object_dtype,
        ),
        (
            pd.DataFrame({"some_col": [1, 1.0, 3.0]}),
            pd.DataFrame({"some_col": [1, 1.0, 3.0]}),
            is_float_dtype,
        ),
    ],
    ids=[
        "mixed_dtype_str_int",
        "all_str",
        "mixed_dtype_float_int_str",
        "mixed_dtype_all_int",
        "int_none",
        "int_nan",
        "mixed_dtype_nan_str",
        "int_float",
    ],
)
class TestConvertMixedDtypes:
    def test_that__convert_df_with_mixed_dtypes_gets_expected_output(
        self, test_input, expected_output, expected_dtype
    ):
        # Create your in memory BytesIO file.
        output = BytesIO()
        test_input.to_csv(output, index=False)
        output.seek(0)  # Contains the CSV in memory file.

        df = transform._convert_df_with_mixed_dtypes(
            {"filepath_or_buffer": output, "index_col": False}
        )
        pd.testing.assert_frame_equal(
            df.reset_index(drop=True), expected_output.reset_index(drop=True)
        )
        assert expected_dtype(df["some_col"])

    def test_that__convert_df_with_mixed_dtypes_catches_pandas_exception(
        self, test_input, expected_output, expected_dtype
    ):
        # Create your in memory BytesIO file.
        output = BytesIO()
        test_input.to_csv(output, index=False)
        output.seek(0)  # Contains the CSV in memory file.

        with mock.patch.object(pd, "read_csv") as mock_read_csv:
            transform._convert_df_with_mixed_dtypes(
                {"filepath_or_buffer": output, "index_col": False}
            )
            mock_read_csv.assert_called_once_with(
                filepath_or_buffer=output,
                index_col=False,
                low_memory=True,
            )
