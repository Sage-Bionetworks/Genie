"""Test genie.extract module"""

from unittest.mock import patch, MagicMock

import pandas as pd
import pytest
import synapseclient

from genie import extract


class Argparser:
    def asDataFrame(self):
        database_dict = {"Database": ["centerMapping"], "Id": ["syn123"]}
        databasetosynid_mappingdf = pd.DataFrame(database_dict)
        return databasetosynid_mappingdf


ENTITY = synapseclient.Project("foo", annotations={"dbMapping": ["syn1234"]})


@pytest.mark.parametrize(
    "test,staging,synid",
    [
        (False, False, "syn10967259"),
        (False, True, "syn12094210"),
        (True, False, "syn11600968"),
    ],
)
def test_get_synid_database_mappingdf(syn, test, staging, synid):
    """
    Tests getting database mapping config
    no flags
    staging flag
    test flag
    """
    arg = Argparser()
    with patch.object(syn, "get", return_value=ENTITY), patch.object(
        extract, "get_syntabledf", return_value=arg.asDataFrame()
    ) as patch_gettabledf:
        df = extract._get_synid_database_mappingdf(syn, project_id=None)
        patch_gettabledf.assert_called_once_with(
            syn, "SELECT * FROM {}".format(ENTITY.dbMapping[0])
        )
        assert df.equals(arg.asDataFrame())


def test_get_syntabledf(syn):
    """
    Test get_syntabledf queries synapse tables and returns a dataframe
    """
    arg = Argparser()
    with patch.object(syn, "tableQuery", return_value=arg) as patch_syn_tablequery:
        querystring = "select * from foo"
        df = extract.get_syntabledf(syn, querystring)
        patch_syn_tablequery.assert_called_once_with(querystring)
        assert df.equals(arg.asDataFrame())


def test_notnone_get_oncotree_link(syn, genie_config):
    """Test link passed in by user is used"""
    url = "https://www.synapse.org"
    link = extract._get_oncotreelink(syn, genie_config, oncotree_link=url)
    assert link == url


def test_none__getoncotreelink(syn, genie_config):
    """Test oncotree link is gotten"""
    url = "https://www.synapse.org"
    link = synapseclient.File("foo", parentId="foo", externalURL=url)
    with patch.object(syn, "get", return_value=link) as patch_synget:
        oncolink = extract._get_oncotreelink(syn, genie_config)
        patch_synget.assert_called_once_with(genie_config["oncotreeLink"])
        assert oncolink == url


def test_that_get_genie_config_has_expected_calls(syn):
    mock_project = MagicMock()
    mock_project.annotations = {"dbMapping": ["db_syn_id"]}
    mock_db_mapping_config = {
        "errorTracker": "syn11601244",
        "centerMapping": "syn11601248",
        "processTracker": "syn11604890",
    }
    mock_center_config = pd.DataFrame(
        {
            "center": ["TEST"],
            "inputSynId": ["syn11601340"],
            "stagingSynId": ["syn11601342"],
            "errorsSynId": ["syn53239081"],
            "release": [True],
            "mutationInCisFilter": ["ON"],
        }
    )
    expected_genie_config = {
        "errorTracker": "syn11601244",
        "centerMapping": "syn11601248",
        "processTracker": "syn11604890",
        "center_config": {
            "TEST": {
                "center": "TEST",
                "inputSynId": "syn11601340",
                "stagingSynId": "syn11601342",
                "errorsSynId": "syn53239081",
                "release": True,
                "mutationInCisFilter": "ON",
            },
        },
    }

    with patch.object(syn, "get", return_value=mock_project) as patch_get, patch.object(
        extract, "_get_database_mapping_config", return_value=mock_db_mapping_config
    ) as patch_get_db_mapping, patch.object(
        extract, "get_syntabledf", return_value=mock_center_config
    ) as patch_get_center_mapping:
        result = extract.get_genie_config(syn, project_id="syn7208886")

        patch_get.assert_called_once_with("syn7208886")
        patch_get_db_mapping.assert_called_once_with(syn=syn, synid="db_syn_id")
        patch_get_center_mapping.assert_called_once_with(
            syn=syn, query_string=f"SELECT * FROM syn11601248 where release is true"
        )
        assert result == expected_genie_config
