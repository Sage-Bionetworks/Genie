# import pytest
# import sys
import mock
import os

from genie import database_to_staging
import pandas as pd
import synapseclient

syn = synapseclient.Synapse()
fileviewSynId = "syn12345"
genieVersion = "vTEST"
consortiumReleaseSynId = "syn2222"


def test_store_gene_panel_files():
    current_release_staging = False

    data_gene_panel = pd.DataFrame({'mutations': ['PANEL1']})

    class gene_panel():
        def asDataFrame(self):
            gene_paneldf = pd.DataFrame({'id': ['syn3333']})
            return gene_paneldf

    with mock.patch.object(
            syn, "tableQuery",
            return_value=gene_panel()) as patch_syn_table_query,\
        mock.patch.object(
            database_to_staging, "storeFile",
            return_value=synapseclient.Entity()) as patch_storefile,\
        mock.patch.object(
            syn, "get",
            return_value=synapseclient.Entity(
                path="/foo/bar/PANEL1.txt")) as patch_syn_get,\
        mock.patch.object(
            os, "rename") as patch_os_rename:

        database_to_staging.store_gene_panel_files(
            syn,
            fileviewSynId,
            genieVersion,
            data_gene_panel,
            consortiumReleaseSynId,
            current_release_staging)

        patch_syn_table_query.assert_called_once_with(
            "select id from %s where cBioFileFormat = 'genePanel' "
            "and fileStage = 'staging'" % fileviewSynId)

        patch_storefile.assert_called_once_with(
            syn,
            os.path.join(
                database_to_staging.GENIE_RELEASE_DIR,
                "PANEL1_vTEST.txt"),
            parent=consortiumReleaseSynId,
            genieVersion=genieVersion,
            name="PANEL1.txt",
            cBioFileFormat="genePanel",
            staging=current_release_staging)

        patch_syn_get.assert_called_once_with('syn3333')
        patch_os_rename.assert_called_once_with(
            "/foo/bar/PANEL1.txt",
            os.path.join(
                database_to_staging.GENIE_RELEASE_DIR,
                "PANEL1_vTEST.txt"))
