"""Test genie.transform module"""
from unittest.mock import patch

import pandas as pd
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
        ],
        ids=["float_w_na", "int_w_na", "string_w_na"],
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
        ],
        ids=["float_no_na", "string_no_na"],
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
