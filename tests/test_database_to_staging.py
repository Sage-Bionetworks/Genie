"""Tests database to staging functions"""
import os
from unittest import mock
from unittest.mock import patch

import pandas as pd
import synapseclient

from genie import database_to_staging, extract, load

FILEVIEW_SYNID = "syn12345"
GENIE_VERSION = "vTEST"
CONSORTIUM_SYNID = "syn2222"


class Tablequerydf:
    """tableQuery.asDataFrame() class"""

    def __init__(self, df):
        self.df = df

    def asDataFrame(self):
        return self.df


def test_store_gene_panel_files(syn):
    data_gene_panel = pd.DataFrame({"mutations": ["PANEL1"]})
    gene_paneldf = pd.DataFrame({"id": ["syn3333"]})

    with mock.patch.object(
        syn, "tableQuery", return_value=Tablequerydf(gene_paneldf)
    ) as patch_syn_table_query, mock.patch.object(
        load, "store_file", return_value=synapseclient.Entity()
    ) as patch_storefile, mock.patch.object(
        syn,
        "get",
        return_value=synapseclient.Entity(path="/foo/bar/PANEL1.txt", versionNumber=2),
    ) as patch_syn_get, mock.patch.object(
        os, "rename"
    ) as patch_os_rename:
        database_to_staging.store_gene_panel_files(
            syn,
            FILEVIEW_SYNID,
            GENIE_VERSION,
            data_gene_panel,
            CONSORTIUM_SYNID,
            ["TEST"],
        )
        patch_syn_table_query.assert_has_calls(
            [
                mock.call("select * from syn12345 limit 1"),
                mock.call(
                    "select id from %s where cBioFileFormat = 'genePanel' "
                    "and fileStage = 'staging' and "
                    "name not in ('data_gene_panel_TEST.txt')" % FILEVIEW_SYNID
                ),
            ]
        )
        assert patch_syn_table_query.call_count == 2
        patch_storefile.assert_called_once_with(
            syn=syn,
            filepath=os.path.join(database_to_staging.GENIE_RELEASE_DIR, "PANEL1.txt"),
            parentid=CONSORTIUM_SYNID,
            version_comment=GENIE_VERSION,
            name="PANEL1.txt",
            annotations={"cBioFileFormat": "genePanel"},
            used="syn3333.2",
        )

        patch_syn_get.assert_called_once_with("syn3333")
        patch_os_rename.assert_called_once_with(
            "/foo/bar/PANEL1.txt",
            os.path.join(database_to_staging.GENIE_RELEASE_DIR, "PANEL1.txt"),
        )


def test_store_assay_info_files(syn):
    """Tests storing of assay information file"""
    assay_infodf = pd.DataFrame({"library_strategy": ["WXS"], "SEQ_ASSAY_ID": ["A"]})
    clinicaldf = pd.DataFrame({"SEQ_ASSAY_ID": ["A"]})
    database_to_staging.GENIE_RELEASE_DIR = "./"
    path = os.path.join(database_to_staging.GENIE_RELEASE_DIR, "assay_information.txt")
    with patch.object(
        syn, "create_snapshot_version", return_value=2
    ) as patch_create_version, patch.object(
        extract, "get_syntabledf", return_value=assay_infodf
    ) as patch_table_query, patch.object(
        load, "store_file", return_value=synapseclient.Entity()
    ) as patch_storefile:
        wes_ids = database_to_staging.store_assay_info_files(
            syn, GENIE_VERSION, FILEVIEW_SYNID, clinicaldf, CONSORTIUM_SYNID
        )
        patch_create_version.assert_called_once_with(
            FILEVIEW_SYNID, comment=GENIE_VERSION
        )
        patch_table_query.assert_called_once_with(
            syn, f"select * from {FILEVIEW_SYNID} where SEQ_ASSAY_ID in ('A')"
        )
        patch_storefile.assert_called_once_with(
            syn=syn,
            filepath=path,
            parentid=CONSORTIUM_SYNID,
            version_comment=GENIE_VERSION,
            name="assay_information.txt",
            used=f"{FILEVIEW_SYNID}.2",
        )
        assert wes_ids == ["A"]
