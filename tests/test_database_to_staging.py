"""Tests database to staging functions"""

import os
import subprocess
from unittest import mock
from unittest.mock import patch

import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal
import pytest
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


def get_run_genie_filters_test_cases():
    return [
        {
            "name": "staging_true",
            "genie_version": "STAGING",
            "variant_filtering_synId": "synZZZZZ",
            "clinicaldf": pd.DataFrame(dict(SAMPLE_ID=["SAGE-1"])),
            "beddf": pd.DataFrame(dict(SAMPLE_ID=["TEST-1"])),
            "center_mappingdf": pd.DataFrame(dict(vcf2maf=["synZZZZZ"])),
            "processing_date": "TEST_DATE",
            "skip_mutationsincis": True,
            "consortium_release_cutoff": "CUTOFF_DATE",
            "test": False,
            "current_release_staging": True,
        },
        {
            "name": "testing_true",
            "genie_version": "TEST",
            "variant_filtering_synId": "synZZZZZ",
            "clinicaldf": pd.DataFrame(dict(SAMPLE_ID=["SAGE-1"])),
            "beddf": pd.DataFrame(dict(SAMPLE_ID=["TEST-1"])),
            "center_mappingdf": pd.DataFrame(dict(vcf2maf=["synZZZZZ"])),
            "processing_date": "TEST_DATE",
            "skip_mutationsincis": True,
            "consortium_release_cutoff": "CUTOFF_DATE",
            "test": False,
            "current_release_staging": True,
        },
        {
            "name": "production_true",
            "genie_version": "PROD",
            "variant_filtering_synId": "synZZZZZ",
            "clinicaldf": pd.DataFrame(dict(SAMPLE_ID=["SAGE-1"])),
            "beddf": pd.DataFrame(dict(SAMPLE_ID=["TEST-1"])),
            "center_mappingdf": pd.DataFrame(dict(vcf2maf=["synZZZZZ"])),
            "processing_date": "TEST_DATE",
            "skip_mutationsincis": False,
            "consortium_release_cutoff": "CUTOFF_DATE",
            "test": False,
            "current_release_staging": False,
        },
    ]


@pytest.mark.parametrize(
    "test_cases", get_run_genie_filters_test_cases(), ids=lambda x: x["name"]
)
def test_that_run_genie_filters_has_expected_calls(syn, test_cases):
    with patch.object(
        database_to_staging,
        "runMAFinBED",
        return_value=pd.DataFrame(dict(VARIANTS_TO_REMOVE=["GENIE-1", "GENIE-2"])),
    ) as patch_run_mafinbed, patch.object(
        database_to_staging,
        "mutation_in_cis_filter",
        return_value=(set(["GENIE-SAMPLE-1"]), set(["GENIE-SAMPLE-4"])),
    ) as patch_mut_in_cis, patch.object(
        database_to_staging, "no_genepanel_filter", return_value=set(["GENIE-SAMPLE-2"])
    ) as patch_no_genepanel_filter, patch.object(
        database_to_staging, "seq_date_filter", return_value=set(["GENIE-SAMPLE-3"])
    ) as patch_seq_date_filter:
        filters_results = database_to_staging.run_genie_filters(
            syn,
            test_cases["genie_version"],
            test_cases["variant_filtering_synId"],
            test_cases["clinicaldf"],
            test_cases["beddf"],
            test_cases["center_mappingdf"],
            test_cases["processing_date"],
            test_cases["skip_mutationsincis"],
            test_cases["consortium_release_cutoff"],
            test_cases["test"],
            test_cases["current_release_staging"],
        )

        patch_run_mafinbed.assert_called_once_with(
            syn,
            test_cases["center_mappingdf"],
            test=test_cases["test"],
            staging=test_cases["current_release_staging"],
            genieVersion=test_cases["genie_version"],
        )

        patch_mut_in_cis.assert_called_once_with(
            syn,
            test_cases["skip_mutationsincis"],
            test_cases["variant_filtering_synId"],
            test_cases["center_mappingdf"],
            genieVersion=test_cases["genie_version"],
            test=test_cases["test"],
            staging=test_cases["current_release_staging"],
        )

        patch_no_genepanel_filter.assert_called_once_with(
            test_cases["clinicaldf"],
            test_cases["beddf"],
        )

        patch_seq_date_filter.assert_called_once_with(
            test_cases["clinicaldf"],
            test_cases["processing_date"],
            test_cases["consortium_release_cutoff"],
        )

        assert_frame_equal(
            filters_results[0],
            pd.DataFrame(dict(VARIANTS_TO_REMOVE=["GENIE-1", "GENIE-2"])),
        )
        assert filters_results[1] == set(
            ["GENIE-SAMPLE-2", "GENIE-SAMPLE-3", "GENIE-SAMPLE-1"]
        )
        assert filters_results[2] == set(["GENIE-SAMPLE-1", "GENIE-SAMPLE-2"])
        assert filters_results[3] == set(["GENIE-SAMPLE-4"])


@pytest.mark.parametrize(
    "test, staging",
    [(True, False), (False, True), (False, False), (True, True)],
    ids=["testing_mode", "staging_mode", "prod_mode", "both_test_and_staging_mode"],
)
def test_that_runMAFinBED_calls_expected_calls(syn, test, staging):
    with patch.object(os.path, "dirname", return_value="test_file_dir/"), patch.object(
        database_to_staging, "get_run_maf_in_bed_script_cmd", return_value="test_cmd"
    ) as patch_get_cmd, patch.object(
        subprocess, "check_call"
    ) as patch_check_call, patch.object(
        database_to_staging,
        "store_filtered_variants",
        return_value=pd.DataFrame(
            dict(
                Chromosome=["1", "X", "Y"],
                removeVariants=[
                    "1 1278471 127818 A G GENIE-1",
                    "X 1278471 127818 A G GENIE-1",
                    "Y 1278471 127818 A G GENIE-1",
                ],
            )
        ),
    ) as patch_store_filtered_variants:

        # setting some testing vars
        test_script_dir = "test_file_dir/"
        test_notinbed_file_path = os.path.join(test_script_dir, "../R/notinbed.csv")
        test_center_mapping_df = pd.DataFrame(dict(vcf2maf=["synZZZZZ"]))
        test_genie_ver = "GENIE_VERSION"

        result = database_to_staging.runMAFinBED(
            syn,
            center_mappingdf=test_center_mapping_df,
            test=test,
            staging=staging,
            genieVersion="GENIE_VERSION",
        )
        patch_get_cmd.assert_called_once_with(
            notinbed_file = test_notinbed_file_path,
            script_dir = test_script_dir,
            test=test,
            staging=staging,
        )
        patch_check_call.assert_called_once_with("test_cmd")
        patch_store_filtered_variants.assert_called_once_with(
            syn=syn,
            notinbed_file=test_notinbed_file_path,
            center_mapping_df=test_center_mapping_df,
            genie_version=test_genie_ver,
        )
        assert_series_equal(
            result,
            pd.Series(
                [
                    "1 1278471 127818 A G GENIE-1",
                    "X 1278471 127818 A G GENIE-1",
                    "Y 1278471 127818 A G GENIE-1",
                ],
                name="removeVariants",
            ),
        )


@pytest.mark.parametrize(
    "test, staging, expected_cmd",
    [
        (
            True,
            False,
            [
                "Rscript",
                "/home/test_dir/Genie/../R/MAFinBED.R",
                "/home/test_dir/Genie/../R/notinbed.csv",
                "--testing",
            ],
        ),
        (
            False,
            True,
            [
                "Rscript",
                "/home/test_dir/Genie/../R/MAFinBED.R",
                "/home/test_dir/Genie/../R/notinbed.csv",
                "--staging",
            ],
        ),
        (
            False,
            False,
            [
                "Rscript",
                "/home/test_dir/Genie/../R/MAFinBED.R",
                "/home/test_dir/Genie/../R/notinbed.csv",
            ],
        ),
        (
            True,
            True,
            [
                "Rscript",
                "/home/test_dir/Genie/../R/MAFinBED.R",
                "/home/test_dir/Genie/../R/notinbed.csv",
                "--testing"
            ],
        ),
    ],
    ids=["testing_mode", "staging_mode", "prod_mode", "both_testing_and_staging_mode"],
)
def test_that_get_run_maf_in_bed_script_cmd_returns_correct_cmd(
    test, staging, expected_cmd
):
    script_dir = "/home/test_dir/Genie/"
    cmd = database_to_staging.get_run_maf_in_bed_script_cmd(
        notinbed_file=os.path.join(script_dir, "../R/notinbed.csv"),
        script_dir=script_dir,
        test=test,
        staging=staging,
    )
    assert cmd == expected_cmd
