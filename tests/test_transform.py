"""Test genie.transform module"""

import os
from io import BytesIO
from unittest import mock

import pandas as pd
from pandas.api.types import (
    is_object_dtype,
)
import pytest

from genie import transform


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
    if kwargs["low_memory"]:
        raise pd.errors.DtypeWarning("some_warning")


class TestConvertMixedDtypes:
    @pytest.fixture(scope="session")
    def test_input(self):
        input = pd.DataFrame({"some_col": [1, "Val2", "1"]})
        output = BytesIO()
        input.to_csv(output, index=False)
        output.seek(0)
        yield output

    @pytest.fixture(scope="session")
    def test_mixed_dtype_input(self):
        input = pd.DataFrame(
            {
                "some_col": ([1.0] * 100000 + ["a"] * 100000 + [float("nan")] * 100000),
                "some_col2": (
                    [1.0] * 100000 + ["b"] * 100000 + [float("nan")] * 100000
                ),
            }
        )
        input.to_csv("test_mixed_dtype_input.csv", index=False)
        yield "test_mixed_dtype_input.csv"
        os.remove("test_mixed_dtype_input.csv")

    def test_that__convert_df_with_mixed_dtypes_gets_expected_output_with_normal_input(
        self, test_input
    ):
        df = transform._convert_df_with_mixed_dtypes(
            {"filepath_or_buffer": test_input, "index_col": False}
        )
        pd.testing.assert_frame_equal(
            df.reset_index(drop=True), pd.DataFrame({"some_col": ["1", "Val2", "1"]})
        )
        assert is_object_dtype(df["some_col"])

    def test_that__convert_df_with_mixed_dtypes_gets_expected_output_with_large_mixed_dtype_input(
        self, test_mixed_dtype_input
    ):
        df = transform._convert_df_with_mixed_dtypes(
            {"filepath_or_buffer": test_mixed_dtype_input, "index_col": False}
        )
        pd.testing.assert_frame_equal(
            df.reset_index(drop=True),
            pd.DataFrame(
                {
                    "some_col": (
                        ["1.0"] * 100000 + ["a"] * 100000 + [float("nan")] * 100000
                    ),
                    "some_col2": (
                        ["1.0"] * 100000 + ["b"] * 100000 + [float("nan")] * 100000
                    ),
                }
            ),
        )
        assert is_object_dtype(df["some_col"])

    def test_that__convert_df_with_mixed_dtypes_calls_read_csv_once_if_no_exception(
        self, test_input
    ):
        with mock.patch.object(pd, "read_csv") as mock_read_csv:
            transform._convert_df_with_mixed_dtypes(
                {"filepath_or_buffer": test_input, "index_col": False}
            )
            mock_read_csv.assert_called_once_with(
                filepath_or_buffer=test_input,
                index_col=False,
                low_memory=True,
            )

    def test_that__convert_df_with_mixed_dtypes_catches_mixed_dtype_exception(
        self, test_input
    ):
        with mock.patch.object(
            pd, "read_csv", side_effect=read_csv_side_effect
        ) as mock_read_csv:
            transform._convert_df_with_mixed_dtypes(
                {"filepath_or_buffer": test_input, "index_col": False}
            )
            mock_read_csv.assert_has_calls(
                [
                    mock.call(
                        filepath_or_buffer=test_input,
                        index_col=False,
                        low_memory=True,
                    ),
                    mock.call(
                        filepath_or_buffer=test_input,
                        index_col=False,
                        low_memory=False,
                        engine="c",
                    ),
                ],
                any_order=False,
            )


def get__convert_values_to_na_test_cases():
    return [
        {
            "name": "single_col",
            "input_df": pd.DataFrame({"col1": ["N/A", "SAGE-SAGE-1", "-nan"]}),
            "values_to_replace": ["N/A", "-nan"],
            "columns_to_convert": ["col1"],
            "expected_df": pd.DataFrame({"col1": [None, "SAGE-SAGE-1", None]}),
        },
        {
            "name": "multi_col",
            "input_df": pd.DataFrame(
                {"col1": ["N/A", "-nan", "1"], "col2": ["-nan", "N/A", "2"]}
            ),
            "values_to_replace": ["N/A", "-nan"],
            "columns_to_convert": ["col1", "col2"],
            "expected_df": pd.DataFrame(
                {"col1": [None, None, "1"], "col2": [None, None, "2"]}
            ),
        },
        {
            "name": "empty_cols",
            "input_df": pd.DataFrame(columns=["col1", "col2"]),
            "values_to_replace": ["NaN", "NA", "nan"],
            "columns_to_convert": ["col1", "col2"],
            "expected_df": pd.DataFrame(columns=["col1", "col2"]),
        },
        {
            "name": "empty df",
            "input_df": pd.DataFrame(),
            "values_to_replace": ["NaN"],
            "columns_to_convert": ["col3"],
            "expected_df": pd.DataFrame(),
        },
    ]


@pytest.mark.parametrize(
    "test_cases", get__convert_values_to_na_test_cases(), ids=lambda x: x["name"]
)
def test_that__convert_values_to_na_returns_expected(test_cases):
    result = transform._convert_values_to_na(
        test_cases["input_df"],
        test_cases["values_to_replace"],
        test_cases["columns_to_convert"],
    )
    pd.testing.assert_frame_equal(result, test_cases["expected_df"])


def test_that__convert_values_to_na_throws_key_error_if_nonexistent_cols():
    input_df = pd.DataFrame({"col1": ["N/A", "SAGE-SAGE-1", "-nan"]})
    values_to_replace = "NaN"
    columns_to_convert = "col3"
    with pytest.raises(KeyError):
        transform._convert_values_to_na(input_df, values_to_replace, columns_to_convert)
