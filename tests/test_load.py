from unittest.mock import patch

import synapseclient
from synapseclient.core.exceptions import SynapseTimeoutError

from genie import load, __version__


def test_store_file(syn):
    """Test storing of file"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234")
        patch_store.assert_called_once_with(
            synapseclient.File("full/path", parentId="syn1234", versionComment=None),
            used=None,
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_file_version_comment(syn):
    """Test storing of file"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234", version_comment="test")
        patch_store.assert_called_once_with(
            synapseclient.File("full/path", parentId="syn1234", versionComment="test"),
            used=None,
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_file_version_provenance(syn):
    """Test storing of file"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234", used="temp")
        patch_store.assert_called_once_with(
            synapseclient.File("full/path", parentId="syn1234", versionComment=None),
            used="temp",
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_file_annotations(syn):
    """Test storing of file with annotations"""
    with patch.object(syn, "store") as patch_store:
        load.store_file(syn, "full/path", "syn1234", annotations={"test": "test"})
        patch_store.assert_called_once_with(
            synapseclient.File(
                "full/path",
                parentId="syn1234",
                versionComment=None,
                annotations={"test": "test"},
            ),
            used=None,
            executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}",
        )


def test_store_table(syn):
    """Test storing of table"""
    with patch.object(syn, "store") as patch_store:
        load.store_table(syn, "full/path", "syn1234")
        patch_store.assert_called_once()


def test_store_table_error(syn):
    """Test storing of table catches and passes error"""
    with patch.object(syn, "store", side_effect=SynapseTimeoutError) as patch_store:
        load.store_table(syn, "full/path", "syn1234")
        patch_store.assert_called_once()
