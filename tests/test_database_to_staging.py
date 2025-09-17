"""Tests database to staging functions"""

import logging
import os
import subprocess
from unittest import mock
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import synapseclient
from genie import database_to_staging, extract, load, process_functions
from pandas.testing import assert_frame_equal, assert_series_equal

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

    with (
        mock.patch.object(
            syn, "tableQuery", return_value=Tablequerydf(gene_paneldf)
        ) as patch_syn_table_query,
        mock.patch.object(
            load, "store_file", return_value=synapseclient.Entity()
        ) as patch_storefile,
        mock.patch.object(
            syn,
            "get",
            return_value=synapseclient.Entity(
                path="/foo/bar/PANEL1.txt", versionNumber=2
            ),
        ) as patch_syn_get,
        mock.patch.object(os, "rename") as patch_os_rename,
    ):
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
    with (
        patch.object(
            syn, "create_snapshot_version", return_value=2
        ) as patch_create_version,
        patch.object(
            extract, "get_syntabledf", return_value=assay_infodf
        ) as patch_table_query,
        patch.object(
            load, "store_file", return_value=synapseclient.Entity()
        ) as patch_storefile,
    ):
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


@pytest.mark.parametrize(
    "input_data, filter_col, expected_result",
    [
        (
            pd.DataFrame(
                dict(
                    SV_STATUS=["GERMLINE", "GERMLINE"], Sample_ID=["GENIE-1", "GENIE-2"]
                )
            ),
            "SV_STATUS",
            pd.DataFrame(columns=["SV_STATUS", "Sample_ID"]),
        ),
        (
            pd.DataFrame(
                dict(
                    SV_STATUS=["GERMLINE", "SOMATIC"], Sample_ID=["GENIE-1", "GENIE-2"]
                )
            ),
            "SV_STATUS",
            pd.DataFrame(dict(SV_STATUS=["SOMATIC"], Sample_ID=["GENIE-2"])),
        ),
        (
            pd.DataFrame(
                dict(SV_STATUS=["SOMATIC", "SOMATIC"], Sample_ID=["GENIE-1", "GENIE-2"])
            ),
            "SV_STATUS",
            pd.DataFrame(
                dict(SV_STATUS=["SOMATIC", "SOMATIC"], Sample_ID=["GENIE-1", "GENIE-2"])
            ),
        ),
        (
            pd.DataFrame(
                dict(
                    SV_Status=["GERMLINE", "SOMATIC"], Sample_ID=["GENIE-1", "GENIE-2"]
                )
            ),
            "SV_STATUS",
            pd.DataFrame(dict(SV_Status=["SOMATIC"], Sample_ID=["GENIE-2"])),
        ),
    ],
    ids=["all_germline", "some_germline", "no_germline", "diff_status_col_case"],
)
def test_that_filter_out_germline_variants_returns_expected(
    input_data, filter_col, expected_result
):
    result = database_to_staging.filter_out_germline_variants(input_data, filter_col)
    assert_frame_equal(result, expected_result, check_index_type=False)


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
            "test": True,
            "current_release_staging": False,
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
    with (
        patch.object(
            database_to_staging,
            "runMAFinBED",
            return_value=pd.DataFrame(dict(VARIANTS_TO_REMOVE=["GENIE-1", "GENIE-2"])),
        ) as patch_run_mafinbed,
        patch.object(
            database_to_staging,
            "mutation_in_cis_filter",
            return_value=(set(["GENIE-SAMPLE-1"]), set(["GENIE-SAMPLE-4"])),
        ) as patch_mut_in_cis,
        patch.object(
            database_to_staging,
            "no_genepanel_filter",
            return_value=set(["GENIE-SAMPLE-2"]),
        ) as patch_no_genepanel_filter,
        patch.object(
            database_to_staging, "seq_date_filter", return_value=set(["GENIE-SAMPLE-3"])
        ) as patch_seq_date_filter,
    ):
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
    with (
        patch.object(os.path, "dirname", return_value="test_file_dir/"),
        patch.object(
            database_to_staging,
            "get_run_maf_in_bed_script_cmd",
            return_value="test_cmd",
        ) as patch_get_cmd,
        patch.object(subprocess, "check_call") as patch_check_call,
        patch.object(
            database_to_staging,
            "store_maf_in_bed_filtered_variants",
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
        ) as patch_store_maf_in_bed_filtered_variants,
    ):
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
            notinbed_file=test_notinbed_file_path,
            script_dir=test_script_dir,
            test=test,
            staging=staging,
        )
        patch_check_call.assert_called_once_with("test_cmd")
        patch_store_maf_in_bed_filtered_variants.assert_called_once_with(
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
                "--testing",
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


def test_that_store_maf_in_bed_filtered_variants_has_expected_calls(syn):
    mock_removed_variants_df = pd.DataFrame(
        dict(
            Chromosome=[2, 3],
            Start_Position=[214142, 214142],
            End_Position=[214142, 214142],
            Reference_Allele=["G", "A"],
            Tumor_Seq_Allele2=["T", "A"],
            Tumor_Sample_Barcode=["GENIE-SAGE-1", "GENIE-TEST-1"],
            Center=["SAGE", "TEST"],
        )
    )
    mock_center_mapping_df = pd.DataFrame(
        dict(center=["SAGE", "TEST"], stagingSynId=["syn001", "syn002"])
    )
    with (
        patch.object(
            pd, "read_csv", return_value=mock_removed_variants_df
        ) as patch_read_csv,
        patch.object(load, "store_file") as patch_store_file,
        patch.object(os, "unlink") as patch_unlink,
        patch.object(pd.DataFrame, "to_csv") as patch_to_csv,
    ):
        test_file_path = "mafinbed_filtered_variants.csv"
        # Call the function
        result = database_to_staging.store_maf_in_bed_filtered_variants(
            syn=syn,
            notinbed_file="test_bed.csv",
            center_mapping_df=mock_center_mapping_df,
            genie_version="TESTING",
        )

        # Assertions
        patch_read_csv.assert_called_once_with("test_bed.csv")
        patch_to_csv.assert_has_calls(
            [
                mock.call(test_file_path, index=False),
                mock.call(test_file_path, index=False),
            ]
        )

        # Check the store_file call arguments
        patch_store_file.assert_has_calls(
            [
                mock.call(
                    syn=syn,
                    filepath=test_file_path,
                    parentid="syn001",
                    version_comment="TESTING",
                ),
                mock.call(
                    syn=syn,
                    filepath=test_file_path,
                    parentid="syn002",
                    version_comment="TESTING",
                ),
            ]
        )
        patch_unlink.assert_has_calls(
            [
                mock.call(test_file_path),
                mock.call(test_file_path),
            ]
        )

        # make sure the removed_variantsdf matches the expected output
        assert_frame_equal(
            result,
            pd.DataFrame(
                dict(
                    Chromosome=[2, 3],
                    Start_Position=[214142, 214142],
                    End_Position=[214142, 214142],
                    Reference_Allele=["G", "A"],
                    Tumor_Seq_Allele2=["T", "A"],
                    Tumor_Sample_Barcode=["GENIE-SAGE-1", "GENIE-TEST-1"],
                    Center=["SAGE", "TEST"],
                    removeVariants=[
                        "2 214142 214142 G T GENIE-SAGE-1",
                        "3 214142 214142 A A GENIE-TEST-1",
                    ],
                )
            ),
        )


@pytest.mark.parametrize(
    "test, staging, skip_mutations_in_cis",
    [
        (True, False, False),
        (False, False, False),
        (False, True, False),
    ],
    ids=["testing_mode", "prod_mode", "staging_mode"],
)
def test_that_mutation_in_cis_filter_has_expected_calls_when_mutations_in_cis_is_not_skipped(
    syn, test, staging, skip_mutations_in_cis
):
    with (
        patch.object(os.path, "dirname", return_value="test_file_dir/"),
        patch.object(
            database_to_staging,
            "get_mutation_in_cis_filter_script_cmd",
            return_value="test_cmd",
        ) as patch_get_cmd,
        patch.object(subprocess, "check_call") as patch_check_call,
        patch.object(
            database_to_staging, "store_mutation_in_cis_files_to_staging"
        ) as patch_store_mutation_in_cis_files,
        patch.object(
            database_to_staging,
            "get_mutation_in_cis_filtered_samples",
            return_value=pd.Series(
                [
                    "GENIE-SAGE-1",
                    "GENIE-SAGE-2",
                ],
                name="Tumor_Sample_Barcode",
            ),
        ) as patch_get_filtered_samples,
        patch.object(
            database_to_staging,
            "get_mutation_in_cis_flagged_variants",
            return_value=pd.Series(
                [
                    "1 1278471 127818 A G GENIE-1",
                    "X 1278471 127818 A G GENIE-1",
                    "Y 1278471 127818 A G GENIE-1",
                ],
                name="flaggedVariants",
            ),
        ) as patch_get_flagged_variants,
    ):
        # setting some testing vars
        test_center_mapping_df = pd.DataFrame(
            dict(Center=["SAGE"], stagingSynId=["synZZZZZ"])
        )
        test_variant_filtering_synId = "synZZZZZ"
        test_genie_ver = "GENIE_VERSION"

        result = database_to_staging.mutation_in_cis_filter(
            syn=syn,
            skipMutationsInCis=skip_mutations_in_cis,
            center_mappingdf=test_center_mapping_df,
            variant_filtering_synId=test_variant_filtering_synId,
            genieVersion=test_genie_ver,
            test=test,
            staging=staging,
        )
        # check expected calls with expected params
        patch_get_cmd.assert_called_once_with(test=test, staging=staging)
        patch_check_call.assert_called_once_with("test_cmd")
        patch_store_mutation_in_cis_files.assert_called_once_with(
            syn=syn,
            center_mappingdf=test_center_mapping_df,
            variant_filtering_synId=test_variant_filtering_synId,
            genieVersion=test_genie_ver,
        )

        patch_get_filtered_samples.assert_called_once_with(
            syn=syn,
            variant_filtering_synId=test_variant_filtering_synId,
        )
        patch_get_flagged_variants.assert_called_once_with(
            syn=syn,
            variant_filtering_synId=test_variant_filtering_synId,
        )
        # check expected remove_samples gets propagated through
        assert_series_equal(
            result[0],
            pd.Series(
                [
                    "GENIE-SAGE-1",
                    "GENIE-SAGE-2",
                ],
                name="Tumor_Sample_Barcode",
            ),
        )

        # check expected flagged variants gets propagated through
        assert_series_equal(
            result[1],
            pd.Series(
                [
                    "1 1278471 127818 A G GENIE-1",
                    "X 1278471 127818 A G GENIE-1",
                    "Y 1278471 127818 A G GENIE-1",
                ],
                name="flaggedVariants",
            ),
        )


@pytest.mark.parametrize(
    "test, staging, skip_mutations_in_cis",
    [
        (True, False, True),
        (False, False, True),
        (False, True, True),
    ],
    ids=[
        "testing_mode_skip_mutations_in_cis",
        "prod_mode_skip_mutations_in_cis",
        "staging_mode_skip_mutations_in_cis",
    ],
)
def test_that_mutation_in_cis_filter_has_expected_calls_when_mutations_in_cis_is_skipped(
    syn, test, staging, skip_mutations_in_cis
):
    with (
        patch.object(os.path, "dirname", return_value="test_file_dir/"),
        patch.object(
            database_to_staging,
            "get_mutation_in_cis_filter_script_cmd",
            return_value="test_cmd",
        ) as patch_get_cmd,
        patch.object(subprocess, "check_call") as patch_check_call,
        patch.object(
            database_to_staging, "store_mutation_in_cis_files_to_staging"
        ) as patch_store_mutation_in_cis_files,
        patch.object(
            database_to_staging,
            "get_mutation_in_cis_filtered_samples",
            return_value=pd.Series(
                [
                    "GENIE-SAGE-1",
                    "GENIE-SAGE-2",
                ],
                name="Tumor_Sample_Barcode",
            ),
        ) as patch_get_filtered_samples,
        patch.object(
            database_to_staging,
            "get_mutation_in_cis_flagged_variants",
            return_value=pd.Series(
                [
                    "1 1278471 127818 A G GENIE-1",
                    "X 1278471 127818 A G GENIE-1",
                    "Y 1278471 127818 A G GENIE-1",
                ],
                name="flaggedVariants",
            ),
        ) as patch_get_flagged_variants,
    ):
        # setting some testing vars
        test_center_mapping_df = pd.DataFrame(
            dict(Center=["SAGE"], stagingSynId=["synZZZZZ"])
        )
        test_variant_filtering_synId = "synZZZZZ"
        test_genie_ver = "GENIE_VERSION"

        result = database_to_staging.mutation_in_cis_filter(
            syn=syn,
            skipMutationsInCis=skip_mutations_in_cis,
            center_mappingdf=test_center_mapping_df,
            variant_filtering_synId=test_variant_filtering_synId,
            genieVersion=test_genie_ver,
            test=test,
            staging=staging,
        )
        # check expected calls with expected params when
        # mutation in cis is skipped or in staging mode
        patch_get_cmd.assert_not_called()
        patch_check_call.assert_not_called()
        patch_store_mutation_in_cis_files.assert_not_called()

        patch_get_filtered_samples.assert_called_once_with(
            syn=syn,
            variant_filtering_synId=test_variant_filtering_synId,
        )
        patch_get_flagged_variants.assert_called_once_with(
            syn=syn,
            variant_filtering_synId=test_variant_filtering_synId,
        )
        # check expected remove_samples gets propagated through
        assert_series_equal(
            result[0],
            pd.Series(
                [
                    "GENIE-SAGE-1",
                    "GENIE-SAGE-2",
                ],
                name="Tumor_Sample_Barcode",
            ),
        )

        # check expected flagged variants gets propagated through
        assert_series_equal(
            result[1],
            pd.Series(
                [
                    "1 1278471 127818 A G GENIE-1",
                    "X 1278471 127818 A G GENIE-1",
                    "Y 1278471 127818 A G GENIE-1",
                ],
                name="flaggedVariants",
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
                "/home/test_dir/Genie/../R/mergeCheck.R",
                "--testing",
            ],
        ),
        (
            False,
            False,
            [
                "Rscript",
                "/home/test_dir/Genie/../R/mergeCheck.R",
            ],
        ),
        (
            False,
            True,
            ["Rscript", "/home/test_dir/Genie/../R/mergeCheck.R", "--staging"],
        ),
    ],
    ids=["testing_mode", "prod_mode", "staging_mode"],
)
def test_that_get_mutation_in_cis_filter_cmd_returns_correct_cmd(
    test, staging, expected_cmd
):
    with patch.object(os.path, "dirname", return_value="/home/test_dir/Genie/"):
        cmd = database_to_staging.get_mutation_in_cis_filter_script_cmd(
            test=test, staging=staging
        )
        assert cmd == expected_cmd


def test_that_get_mutation_in_cis_filter_cmd_raises_value_error():
    with patch.object(os.path, "dirname", return_value="/home/test_dir/Genie/"):
        with pytest.raises(
            ValueError,
            match="Mutation in cis only available in staging or testing mode not both",
        ):
            database_to_staging.get_mutation_in_cis_filter_script_cmd(
                test=True, staging=True
            )


@pytest.mark.parametrize(
    "flagged_df, expected_flagged_variants",
    [
        (
            pd.DataFrame(
                dict(
                    Chromosome=[2],
                    Start_Position=[214142],
                    HGVSp_Short=[""],
                    Reference_Allele=["G"],
                    Tumor_Seq_Allele2=["T"],
                    Tumor_Sample_Barcode=["GENIE-SAGE-1"],
                )
            ),
            pd.Series(["2 214142  G T GENIE-SAGE-1"], name="flaggedVariants"),
        ),
    ],
    ids=["nonempty_query"],
)
def test_that_get_mutation_in_cis_flagged_variants_returns_expected_flagged_variants(
    syn, flagged_df, expected_flagged_variants
):
    with patch.object(
        extract, "get_syntabledf", return_value=flagged_df
    ) as patch_get_syntabledf:
        result = database_to_staging.get_mutation_in_cis_flagged_variants(
            syn, variant_filtering_synId="synZZZZZ"
        )
        patch_get_syntabledf.assert_called_once_with(
            syn=syn,
            query_string="SELECT * FROM synZZZZZ where "
            "Flag = 'Flag' and Tumor_Sample_Barcode is not null",
        )
        assert_series_equal(result, expected_flagged_variants)


@pytest.mark.parametrize(
    "filtered_df, expected_removed_samples",
    [
        (
            pd.DataFrame(dict(Tumor_Sample_Barcode=["GENIE-SAGE-1"])),
            pd.Series(["GENIE-SAGE-1"], name="Tumor_Sample_Barcode"),
        ),
        (
            pd.DataFrame(
                dict(
                    Tumor_Sample_Barcode=[
                        "GENIE-SAGE-1",
                        "GENIE-SAGE-1",
                        "GENIE-SAGE-2",
                    ]
                )
            ),
            pd.Series(["GENIE-SAGE-1", "GENIE-SAGE-2"], name="Tumor_Sample_Barcode"),
        ),
    ],
    ids=["nonempty_query", "query_with_dups"],
)
def test_that_get_mutation_in_cis_filtered_samples_returns_expected_variants(
    syn, filtered_df, expected_removed_samples
):
    with patch.object(
        extract, "get_syntabledf", return_value=filtered_df
    ) as patch_get_syntabledf:
        result = database_to_staging.get_mutation_in_cis_filtered_samples(
            syn, variant_filtering_synId="synZZZZZ"
        )
        patch_get_syntabledf.assert_called_once_with(
            syn=syn,
            query_string="SELECT Tumor_Sample_Barcode FROM synZZZZZ where "
            "Flag = 'TOSS' and Tumor_Sample_Barcode is not null",
        )
        assert_series_equal(result, expected_removed_samples, check_index=False)


@pytest.fixture
def mock_center_mappingdf():
    yield pd.DataFrame(
        {"center": ["CenterA", "CenterB"], "stagingSynId": ["syn001", "syn002"]}
    )


@pytest.fixture
def mock_merge_check_df():
    yield pd.DataFrame({"Center": ["CenterA", "CenterB", None], "test_col": [1, 2, 3]})


def test_store_mutation_in_cis_files_to_staging(
    syn,
    mock_center_mappingdf,
    mock_merge_check_df,
):
    with (
        patch.object(
            extract, "get_syntabledf", return_value=mock_merge_check_df
        ) as patch_get_syntabledf,
        patch.object(load, "store_file") as patch_store_file,
        patch.object(os, "unlink") as patch_unlink,
        patch.object(pd.DataFrame, "to_csv") as patch_to_csv,
    ):
        test_file_path = "mutationsInCis_filtered_samples.csv"
        # Call the function
        database_to_staging.store_mutation_in_cis_files_to_staging(
            syn, mock_center_mappingdf, "syn123", "TESTING"
        )

        # Assertions
        patch_get_syntabledf.assert_called_once_with(
            syn=syn,
            query_string="select * from syn123 where Center in ('CenterA','CenterB')",
        )
        patch_to_csv.assert_has_calls(
            [
                mock.call(test_file_path, index=False),
                mock.call(test_file_path, index=False),
            ]
        )

        # Check the store_file call arguments
        patch_store_file.assert_has_calls(
            [
                mock.call(
                    syn=syn,
                    filepath=test_file_path,
                    parentid="syn001",
                    version_comment="TESTING",
                ),
                mock.call(
                    syn=syn,
                    filepath=test_file_path,
                    parentid="syn002",
                    version_comment="TESTING",
                ),
            ]
        )
        patch_unlink.assert_has_calls(
            [
                mock.call(test_file_path),
                mock.call(test_file_path),
            ]
        )


def test_store_mutation_in_cis_files_to_staging_no_centers(
    syn,
    mock_center_mappingdf,
):
    with (
        patch.object(extract, "get_syntabledf") as patch_get_syntabledf,
        patch.object(load, "store_file") as patch_store_file,
        patch.object(os, "unlink") as patch_unlink,
    ):
        # Set up mocks
        mock_mergeCheckDf = pd.DataFrame({"Center": [], "test_col": []})
        patch_get_syntabledf.return_value = mock_mergeCheckDf

        # Call the function
        database_to_staging.store_mutation_in_cis_files_to_staging(
            syn, mock_center_mappingdf, "syn123", "TESTING"
        )

        # Assertions
        patch_get_syntabledf.assert_called_once_with(
            syn=syn,
            query_string="select * from syn123 where Center in ('CenterA','CenterB')",
        )
        patch_store_file.assert_not_called()
        patch_unlink.assert_not_called()


@pytest.mark.parametrize(
    "current_release_staging",
    [
        (True),
        (False),
    ],
    ids=["current_release_staging", "other_release_staging"],
)
def test_store_sv_files(syn, current_release_staging):
    svdf = pd.DataFrame(
        dict(
            SV_Status=["SOMATIC", "GERMLINE", "SOMATIC"],
            SAMPLE_ID=["GENIE-1", "GENIE-2", "GENIE-3"],
            CENTER=["C1", "C2", "C1"],
        )
    )

    center_mappingdf = pd.DataFrame(
        {"center": ["C1", "C2"], "stagingSynId": ["syn123", "syn456"]}
    )
    database_to_staging.GENIE_RELEASE_DIR = "./"
    database_to_staging.SV_CENTER_PATH = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR, "data_sv_%s.txt"
    )
    with (
        patch("builtins.open") as mock_open,
        patch.object(
            extract, "get_syntabledf", return_value=svdf
        ) as patch_get_syntabledf,
        patch.object(syn, "create_snapshot_version", return_value=1),
        patch.object(pd.DataFrame, "to_csv") as patch_to_csv,
        patch.object(load, "store_file") as patch_store_file,
        patch.object(
            process_functions, "removePandasDfFloat", return_value="test_sv_text"
        ),
    ):
        # call the function
        sv_sample = database_to_staging.store_sv_files(
            syn,
            release_synid="syn123",
            genie_version="TESTING",
            synid="syn5678",
            keep_for_center_consortium_samples=["GENIE-1", "GENIE-2"],
            keep_for_merged_consortium_samples=["GENIE-1", "GENIE-3"],
            current_release_staging=current_release_staging,
            center_mappingdf=center_mappingdf,
        )

        # validate the calls
        patch_get_syntabledf.assert_called_once_with(syn, "select * from syn5678")
        if not current_release_staging:
            patch_to_csv.assert_has_calls(
                [
                    mock.call(
                        database_to_staging.SV_CENTER_PATH % "C1", sep="\t", index=False
                    ),
                    mock.call(
                        database_to_staging.SV_CENTER_PATH % "C2", sep="\t", index=False
                    ),
                ]
            )

            # Check the store_file call arguments
            patch_store_file.assert_has_calls(
                [
                    mock.call(
                        syn=syn,
                        filepath=database_to_staging.SV_CENTER_PATH % "C1",
                        version_comment="TESTING",
                        parentid="syn123",
                    ),
                    mock.call(
                        syn=syn,
                        filepath=database_to_staging.SV_CENTER_PATH % "C2",
                        version_comment="TESTING",
                        parentid="syn456",
                    ),
                ]
            )

        mock_open.assert_called_once_with(
            os.path.join(database_to_staging.GENIE_RELEASE_DIR, "data_sv.txt"), "w"
        )
        file_handle = mock_open.return_value.__enter__.return_value
        file_handle.write.assert_called_with("test_sv_text")
        patch_store_file.assert_called_with(
            syn=syn,
            filepath=os.path.join(database_to_staging.GENIE_RELEASE_DIR, "data_sv.txt"),
            parentid="syn123",
            version_comment="TESTING",
            name="data_sv.txt",
            used="syn5678.1",
        )
        assert sv_sample == ["GENIE-1", "GENIE-3"]


@pytest.mark.parametrize(
    "wes_seqassayids, expected_output",
    [
        (
            ["ID1"],
            pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-2"],
                    "mutations": ["ID2"],
                    "cna": ["NA"],
                    "sv": ["ID2"],
                }
            ),
        ),
        (
            ["ID3"],
            pd.DataFrame(
                {
                    "SAMPLE_ID": ["GENIE-1", "GENIE-2"],
                    "mutations": ["ID1", "ID2"],
                    "cna": ["ID1", "NA"],
                    "sv": ["ID1", "ID2"],
                }
            ),
        ),
    ],
    ids=["filter_wes_seqassayids", "no_filter_wes_seqassayids"],
)
def test_store_data_gene_matrix(syn, wes_seqassayids, expected_output):
    database_to_staging.GENIE_RELEASE_DIR = "./"
    clinicaldf = pd.DataFrame(
        {"SAMPLE_ID": ["GENIE-1", "GENIE-2"], "SEQ_ASSAY_ID": ["ID1", "ID2"]}
    )
    with (
        patch.object(pd.DataFrame, "to_csv") as patch_to_csv,
        patch.object(load, "store_file") as patch_store_file,
    ):
        # call the function
        data_gene_matrix = database_to_staging.store_data_gene_matrix(
            syn,
            genie_version="TESTING",
            clinicaldf=clinicaldf,
            cna_samples=["GENIE-1"],
            release_synid="syn123",
            wes_seqassayids=wes_seqassayids,
            sv_samples=["GENIE-1", "GENIE-2"],
        )

        # validate the calls
        patch_to_csv.assert_called_once_with(
            "./data_gene_matrix.txt", sep="\t", index=False
        )
        patch_store_file.assert_called_once_with(
            syn=syn,
            filepath=os.path.join(
                database_to_staging.GENIE_RELEASE_DIR, "data_gene_matrix.txt"
            ),
            parentid="syn123",
            version_comment="TESTING",
            name="data_gene_matrix.txt",
        )
        pd.testing.assert_frame_equal(
            data_gene_matrix.reset_index(drop=True), expected_output
        )


@pytest.mark.parametrize(
    "input_col,expected_to_redact_vector,expected_to_redact_pediatric_vector",
    [
        (
            pd.Series([4380, 23725, 32120, 33215]),  # in years: 12, 65, 88, 91,
            pd.Series([False, False, False, True]),
            pd.Series([False, False, False, False]),
        ),
        (
            pd.Series([6570, 32485]),  # in years: 18, 89
            pd.Series([False, False]),
            pd.Series([False, False]),
        ),
        (
            pd.Series(
                [
                    ">32485",
                    "<6570",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                ]  # in years: ">89", "<18"
            ),
            pd.Series([True, False, False, False, False]),
            pd.Series([False, True, False, False, False]),
        ),
        (
            pd.Series(
                [
                    ">32485",
                    "<6570",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                    np.nan,
                ]  # in years: ">89", "<18"
            ),
            pd.Series([True, False, False, False, False, False]),
            pd.Series([False, True, False, False, False, False]),
        ),
    ],
    ids=[
        "redact_89+_not_ped",
        "no_redaction_for_cutoff_value",
        "redact_range_values",
        "no_redaction_for_NAs",
    ],
)
def test_to_redact_interval(
    input_col, expected_to_redact_vector, expected_to_redact_pediatric_vector
):
    # call the function
    to_redact, to_redact_pediatric = database_to_staging._to_redact_interval(input_col)

    # validate the calls
    assert to_redact.equals(expected_to_redact_vector)
    assert to_redact_pediatric.equals(expected_to_redact_pediatric_vector)


@pytest.mark.parametrize(
    "input_col,expected_col",
    [
        (
            pd.Series(
                [
                    4380,
                    23725,
                    32120,
                    33215,
                    32485,
                    6570,
                ]  # in years: 12, 65, 88, 91, 89, 18
            ),
            pd.Series([4380, 23725, 32120, 33215, 32485, 6570]),
        ),
        (
            pd.Series(
                [
                    ">32485",
                    "<6570",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                ]  # in years: ">89", "<18"
            ),
            pd.Series(
                [
                    "cannotReleaseHIPAA",
                    "withheld",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                ]
            ),
        ),
        (
            pd.Series(
                [
                    ">32485",
                    "<6570",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                    np.nan,
                ]  # in years: ">89", "<18"
            ),
            pd.Series(
                [
                    "cannotReleaseHIPAA",
                    "withheld",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                    np.nan,
                ]
            ),
        ),
    ],
    ids=[
        "no_redaction_for_numeric_values",
        "redact_range_values",
        "no_redaction_for_NAs",
    ],
)
def test__redact_year(input_col, expected_col):
    # call the function
    output = database_to_staging._redact_year(input_col)

    # validate the calls
    assert output.equals(expected_col)


@pytest.mark.parametrize(
    "input_col,expected_col",
    [
        (
            pd.Series(
                [
                    4380,
                    23725,
                    32120,
                    33215,
                    32485,
                    6570,
                ]  # in years: 12, 65, 88, 91, 89, 18
            ),
            pd.Series([4380, 23725, 32120, 33215, 32485, 6570]),
        ),
        (
            pd.Series(
                [
                    ">32485",
                    "<6570",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                ]  # in years: ">89", "<18"
            ),
            pd.Series(
                [
                    ">32485",
                    "withheld",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                ]
            ),
        ),
        (
            pd.Series(
                [
                    ">32485",
                    "<6570",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                    np.nan,
                ]  # in years: ">89", "<18"
            ),
            pd.Series(
                [
                    ">32485",
                    "withheld",
                    "Not Collected",
                    "Unknown",
                    "Not Applicable",
                    np.nan,
                ]
            ),
        ),
    ],
    ids=[
        "no_redaction_for_numeric_values",
        "redact_range_values",
        "no_redaction_for_NAs",
    ],
)
def test_redact_ped_year(input_col, expected_col):
    # call the function
    output = database_to_staging._redact_ped_year(input_col)

    # validate the calls
    assert output.equals(expected_col)


@pytest.mark.parametrize(
    "df_col_year1,df_col_year2, expected_to_redact",
    [
        (
            pd.Series([1892, 1993, 1993]),
            pd.Series([2033, 2033, 2009]),  # in years: 141, 40, 16
            pd.Series([True, False, False]),
        ),
        (
            pd.Series(["Not Collected", "Unknown", ">89", "<18"]),
            pd.Series([1992, 1993, 1992, 1993]),
            pd.Series([False, False, False, False]),
        ),
        (
            pd.Series(["Not Collected", np.nan]),
            pd.Series([1992, 1993]),
            pd.Series([False, False]),
        ),
    ],
    ids=[
        "redact_89+_diff",
        "no_redaction_cuz_no_diff_calculated_for_string",
        "no_redaction_cuz_no_diff_calculated_for_NAs",
    ],
)
def test_to_redact_difference(df_col_year1, df_col_year2, expected_to_redact):
    # call the function
    to_redact = database_to_staging._to_redact_difference(df_col_year1, df_col_year2)

    # validate the calls
    assert to_redact.equals(expected_to_redact)


def get_redact_phi_test_cases():
    return [
        {
            "name": "no_redaction",
            "input_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "a", "c", "d"],
                    "SAMPLE_ID": ["a", "b", "c", "d"],
                    "BIRTH_YEAR": [1992, 1993, 1992, 1993],
                    "YEAR_CONTACT": [
                        2033,
                        2030,
                        2010,
                        2082,
                    ],  # difference to BIRTH_YEAR: 41, 37, 18, 89
                    "YEAR_DEATH": [
                        2044,
                        2047,
                        2012,
                        2082,
                    ],  # difference to BIRTH_YEAR: 52, 54, 20, 89
                    "AGE_AT_SEQ_REPORT": [
                        13870,
                        17885,
                        6935,
                        32485,
                    ],  # in years 38*365, 49*365, 19*365, 89*365
                    "INT_CONTACT": [
                        14966,
                        13506,
                        6571,
                        32485,
                    ],  # in years: 41*365+1, 37*365+1, 18*365+1, 89*365
                    "INT_DOD": [
                        18981,
                        19712,
                        7300,
                        32485,
                    ],  # in years: 52*365+1, 54*365+2, 20*365, 89*365
                    "Other_col": ["a", "b", "c", "d"],
                }
            ),
            "expected_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "a", "c", "d"],
                    "SAMPLE_ID": ["a", "b", "c", "d"],
                    "BIRTH_YEAR": [1992, 1993, 1992, 1993],
                    "YEAR_CONTACT": [
                        2033,
                        2030,
                        2010,
                        2082,
                    ],  # difference to BIRTH_YEAR: 41, 37, 18, 89
                    "YEAR_DEATH": [
                        2044,
                        2047,
                        2012,
                        2082,
                    ],  # difference to BIRTH_YEAR: 52, 54, 20, 89
                    "AGE_AT_SEQ_REPORT": [
                        13870,
                        17885,
                        6935,
                        32485,
                    ],  # in years 38*365, 49*365, 19*365, 89*365
                    "INT_CONTACT": [
                        14966,
                        13506,
                        6571,
                        32485,
                    ],  # in years: 41*365+1, 37*365+1, 18*365+1, 89*365
                    "INT_DOD": [
                        18981,
                        19712,
                        7300,
                        32485,
                    ],  # in years: 52*365+1, 54*365+2, 20*365, 89*365
                    "Other_col": ["a", "b", "c", "d"],
                }
            ),
        },
        {
            "name": "redact_89+_not_ped",
            "input_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b", "c", "d", "e", "e"],  # redact b, d, e
                    "SAMPLE_ID": ["a", "b", "c", "d", "e", "f"],
                    "BIRTH_YEAR": [1992, 1992, 1992, 1992, 1992, 1992],
                    "YEAR_CONTACT": [
                        2008,
                        2080,
                        2008,
                        2082,
                        2080,
                        2081,
                    ],  # difference to BIRTH_YEAR: 16, 88, 16, 90, 88, 89
                    "YEAR_DEATH": [
                        2011,
                        2082,
                        2008,
                        2082,
                        2081,
                        2081,
                    ],  # difference to BIRTH_YEAR: 19, 90, 16, 90, 89, 89
                    "AGE_AT_SEQ_REPORT": [
                        6570,
                        32485,
                        5840,
                        32850,
                        32120,
                        32485,
                    ],  # in years: 18*365, 89*365, 16*365, 90*365, 88*365, 89*365
                    "INT_CONTACT": [
                        5841,
                        32121,
                        5841,
                        32851,
                        32121,
                        32485,
                    ],  # in years: 16*365+1, 88*365+1, 16*365+1, 90*365+1, 88*365+1, 89*365
                    "INT_DOD": [
                        6935,
                        32850,
                        5840,
                        32850,
                        32485,
                        32486,
                    ],  # in years: 19*365, 90*365, 16*365, 90*365, 89*365, 89*365+1
                    "Other_col": ["a", "b", "c", "d", "e", "f"],
                }
            ),
            "expected_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b", "c", "d", "e", "e"],
                    "SAMPLE_ID": ["a", "b", "c", "d", "e", "f"],
                    "BIRTH_YEAR": [
                        1992,
                        "cannotReleaseHIPAA",
                        1992,
                        "cannotReleaseHIPAA",
                        1992,
                        "cannotReleaseHIPAA",
                    ],
                    "YEAR_CONTACT": [
                        2008,
                        2080,
                        2008,
                        2082,
                        2080,
                        2081,
                    ],  # difference to BIRTH_YEAR: 16, 88, 16, 90, 88, 89
                    "YEAR_DEATH": [
                        2011,
                        2082,
                        2008,
                        2082,
                        2081,
                        2081,
                    ],  # difference to BIRTH_YEAR: 19, 90, 16, 90, 89, 89
                    "AGE_AT_SEQ_REPORT": [
                        6570,
                        32485,
                        5840,
                        ">32485",
                        32120,
                        32485,
                    ],  # in years: 18*365, 89*365, 16*365, 90*365, 88*365, 89*365
                    "INT_CONTACT": [
                        5841,
                        32121,
                        5841,
                        ">32485",
                        32121,
                        32485,
                    ],  # in years: 16*365+1, 88*365+1, 16*365+1, 90*365+1, 88*365+1, 89*365
                    "INT_DOD": [
                        6935,
                        ">32485",
                        5840,
                        ">32485",
                        32485,
                        ">32485",
                    ],  # in years: 19*365, 90*365, 16*365, 90*365, 89*365, 89*365+1
                    "Other_col": ["a", "b", "c", "d", "e", "f"],
                }
            ),
        },
        {
            "name": "redact_range_values_for_ped",
            "input_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b", "c", "d"],
                    "SAMPLE_ID": ["a", "b", "c", "d"],
                    "BIRTH_YEAR": [1992, 1992, 1992, 1992],
                    "YEAR_CONTACT": [
                        "<18",
                        2080,
                        2008,
                        2082,
                    ],  # difference to BIRTH_YEAR: 16, 88, 16, 90
                    "YEAR_DEATH": [
                        2011,
                        ">89",
                        2008,
                        2082,
                    ],  # difference to BIRTH_YEAR: 19, 90, 16, 90
                    "AGE_AT_SEQ_REPORT": [
                        6570,
                        32485,
                        "<6570",
                        ">32485",
                    ],  # in years: 18*365, 89*365, 16*365, 90*365
                    "INT_CONTACT": [
                        "<6570",
                        32121,
                        5841,
                        32851,
                    ],  # in years: 16*365+1, 88*365+1, 16*365+1, 90*365+1
                    "INT_DOD": [
                        6935,
                        ">32485",
                        "<18",
                        32850,
                    ],  # in years: 19*365, 90*365, "<18", 90*365
                    "Other_col": ["a", "b", "c", "d"],
                }
            ),
            "expected_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b", "c", "d"],
                    "SAMPLE_ID": ["a", "b", "c", "d"],
                    "BIRTH_YEAR": [
                        1992,
                        "cannotReleaseHIPAA",
                        1992,
                        "cannotReleaseHIPAA",
                    ],
                    "YEAR_CONTACT": [
                        "withheld",
                        2080,
                        2008,
                        2082,
                    ],  # difference to BIRTH_YEAR: 16, 88, 16, 90
                    "YEAR_DEATH": [
                        2011,
                        ">89",
                        2008,
                        2082,
                    ],  # difference to BIRTH_YEAR: 19, 90, 16, 90
                    "AGE_AT_SEQ_REPORT": [
                        6570,
                        32485,
                        "withheld",
                        ">32485",
                    ],  # in years: 18*365, 89*365, 16*365, 90*365
                    "INT_CONTACT": [
                        "withheld",
                        32121,
                        5841,
                        ">32485",
                    ],  # in years: 16*365+1, 88*365+1, 16*365+1, 90*365+1
                    "INT_DOD": [
                        6935,
                        ">32485",
                        "withheld",
                        ">32485",
                    ],  # in years: 19*365, 90*365, "<18", 90*365
                    "Other_col": ["a", "b", "c", "d"],
                }
            ),
        },
        {
            "name": "no_redaction_for_NAs",
            "input_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b", "c", "d", "e", "f"],
                    "SAMPLE_ID": ["a", "b", "c", "d", "e", "f"],
                    "BIRTH_YEAR": [1992, 1992, 1992, 1992, 1993, 1994],
                    "YEAR_CONTACT": [
                        "<18",
                        2080,
                        2008,
                        2082,
                        2023,
                        np.nan,
                    ],  # difference to BIRTH_YEAR: 16, 88, 16, 90, 30, 31
                    "YEAR_DEATH": [
                        2011,
                        ">89",
                        2008,
                        2082,
                        2025,
                        2025,
                    ],  # difference to BIRTH_YEAR: 19, 90, 16, 90, 32, 31
                    "AGE_AT_SEQ_REPORT": [
                        6570,
                        32485,
                        "<6570",
                        ">32485",
                        11315,
                        np.nan,
                    ],  # in years: 18*365, 89*365, 16*365, 90*365, 31*365
                    "INT_CONTACT": [
                        "<6570",
                        32121,
                        5841,
                        32851,
                        10951,
                        np.nan,
                    ],  # in years: 16*365+1, 88*365+1, 16*365+1, 90*365+1, 30*365+1
                    "INT_DOD": [
                        6935,
                        ">32485",
                        5840,
                        32850,
                        11680,
                        11315,
                    ],  # in years: 19*365, 90*365, 16*365, 90*365, 32*365, 31*365
                    "Other_col": ["a", "b", "c", "d", "e", "f"],
                }
            ),
            "expected_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b", "c", "d", "e", "f"],
                    "SAMPLE_ID": ["a", "b", "c", "d", "e", "f"],
                    "BIRTH_YEAR": [
                        1992,
                        "cannotReleaseHIPAA",
                        1992,
                        "cannotReleaseHIPAA",
                        1993,
                        1994,
                    ],
                    "YEAR_CONTACT": [
                        "withheld",
                        2080,
                        2008,
                        2082,
                        2023,
                        np.nan,
                    ],  # difference to BIRTH_YEAR: 16, 88, 16, 90, 30
                    "YEAR_DEATH": [
                        2011,
                        ">89",
                        2008,
                        2082,
                        2025,
                        2025,
                    ],  # difference to BIRTH_YEAR: 19, 90, 16, 90, 32, 31
                    "AGE_AT_SEQ_REPORT": [
                        6570,
                        32485,
                        "withheld",
                        ">32485",
                        11315,
                        np.nan,
                    ],  # in years: 18*365, 89*365, 16*365, 90*365, 31*365
                    "INT_CONTACT": [
                        "withheld",
                        32121,
                        5841,
                        ">32485",
                        10951,
                        np.nan,
                    ],  # in years: 16*365+1, 88*365+1, 16*365+1, 90*365+1, 30*365+1
                    "INT_DOD": [
                        6935,
                        ">32485",
                        5840,
                        ">32485",
                        11680,
                        11315,
                    ],  # in years: 19*365, 90*365, 16*365, 90*365, 32*365, 31*365
                    "Other_col": ["a", "b", "c", "d", "e", "f"],
                }
            ),
        },
        {
            "name": "redact_range_birth_year",
            "input_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b"],
                    "SAMPLE_ID": ["a", "b"],
                    "BIRTH_YEAR": ["<18", ">89"],
                    "YEAR_CONTACT": [
                        2011,
                        2010,
                    ],
                    "YEAR_DEATH": [
                        2013,
                        2013,
                    ],
                    "AGE_AT_SEQ_REPORT": [
                        "<18",
                        np.nan,
                    ],
                    "INT_CONTACT": [
                        np.nan,
                        np.nan,
                    ],
                    "INT_DOD": [
                        np.nan,
                        ">89",
                    ],
                    "Other_col": ["a", "b"],
                }
            ),
            "expected_df": pd.DataFrame(
                {
                    "PATIENT_ID": ["a", "b"],
                    "SAMPLE_ID": ["a", "b"],
                    "BIRTH_YEAR": [
                        "withheld",
                        "cannotReleaseHIPAA",
                    ],
                    "YEAR_CONTACT": [
                        2011,
                        2010,
                    ],
                    "YEAR_DEATH": [
                        2013,
                        2013,
                    ],
                    "AGE_AT_SEQ_REPORT": [
                        "withheld",
                        np.nan,
                    ],
                    "INT_CONTACT": [
                        np.nan,
                        np.nan,
                    ],
                    "INT_DOD": [
                        np.nan,
                        ">32485",
                    ],
                    "Other_col": ["a", "b"],
                }
            ),
        },
    ]


@pytest.mark.parametrize(
    "test_cases", get_redact_phi_test_cases(), ids=lambda x: x["name"]
)
def test_redact_phi(test_cases):
    # call the function
    output = database_to_staging.redact_phi(test_cases["input_df"])

    # validate_calss
    assert_frame_equal(output, test_cases["expected_df"], check_dtype=False)
