from unittest.mock import patch

import pandas as pd
import pytest

from genie.example_filetype_format import FileTypeFormat


@pytest.fixture
def filetype_format_class(syn):
    yield FileTypeFormat(syn, "SAGE", ancillary_files=[["mocked"]])


@pytest.mark.parametrize("invalid_ancillary", [[], None])
def test_that_validate_no_cross_file_without_invalid_ancillary(invalid_ancillary, syn):
    filetype_cls = FileTypeFormat(syn, "SAGE", ancillary_files=invalid_ancillary)
    with patch.object(
        FileTypeFormat,
        "_validate",
        return_value=("", "some_warning\n"),
    ) as patch_validate, patch.object(
        FileTypeFormat,
        "read_file",
        return_value=pd.DataFrame(),
    ), patch.object(
        FileTypeFormat,
        "_cross_validate",
    ) as patch_cross_validate:
        result_cls = filetype_cls.validate(filePathList=["something.txt"])
        patch_validate.assert_called_once()
        patch_cross_validate.assert_not_called()
        assert result_cls.warnings == "some_warning\n"
        assert result_cls.errors == ""


def test_that_validate_returns_expected_msg_if__validate_fails(filetype_format_class):
    with patch.object(
        FileTypeFormat,
        "_validate",
        return_value=("some_error", ""),
    ) as patch_validate, patch.object(
        FileTypeFormat,
        "read_file",
        return_value=pd.DataFrame(),
    ), patch.object(
        FileTypeFormat,
        "_cross_validate",
        return_value=("some_cross_error", "some_cross_warning"),
    ) as patch_cross_validate:
        result_cls = filetype_format_class.validate(filePathList=["something.txt"])
        patch_validate.assert_called_once()
        patch_cross_validate.assert_not_called()
        assert result_cls.warnings == ""
        assert result_cls.errors == "some_error"


def test_that_validate_returns_expected_msg_if__validate_passes(filetype_format_class):
    with patch.object(
        FileTypeFormat,
        "_validate",
        return_value=("", "some_warning\n"),
    ) as patch_validate, patch.object(
        FileTypeFormat,
        "read_file",
        return_value=pd.DataFrame(),
    ), patch.object(
        FileTypeFormat,
        "_cross_validate",
        return_value=("some_cross_error", "some_cross_warning"),
    ) as patch_cross_validate:
        result_cls = filetype_format_class.validate(filePathList=["something.txt"])
        patch_validate.assert_called_once()
        patch_cross_validate.assert_called_once()
        assert result_cls.warnings == "some_warning\nsome_cross_warning"
        assert result_cls.errors == "some_cross_error"


def test_that_validate_throws_exception_if_file_read_error(filetype_format_class):
    with patch.object(
        FileTypeFormat,
        "_validate",
    ) as patch_validate, patch.object(
        FileTypeFormat,
        "read_file",
        side_effect=Exception("mocked error"),
    ), patch.object(FileTypeFormat, "_cross_validate") as patch_cross_validate:
        result_cls = filetype_format_class.validate(filePathList=["something.txt"])
        assert result_cls.warnings == ""
        assert result_cls.errors == (
            "The file(s) (['something.txt']) cannot be read. "
            "Original error: mocked error"
        )
        patch_validate.assert_not_called()
        patch_cross_validate.assert_not_called()
