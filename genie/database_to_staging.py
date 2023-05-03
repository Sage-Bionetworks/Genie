#!/usr/local/bin/ python3
"""
Functions for releasing GENIE consortium releases
"""
import logging
import math
import os
import re
import subprocess
from typing import List

import pandas as pd
import pyranges
import synapseclient

from . import extract, load, process_functions, transform

logger = logging.getLogger(__name__)

# GENIE CONSTANTS
# GENIE_RELEASE_DIR = os.path.join(
#    os.path.dirname(os.path.abspath(__file__)),"GENIE_release")
GENIE_RELEASE_DIR = os.path.join(os.path.expanduser("~/.synapseCache"), "GENIE_release")
CASE_LIST_PATH = os.path.join(GENIE_RELEASE_DIR, "case_lists")
CNA_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR, "data_CNA_%s.txt")
SAMPLE_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR, "data_clinical_supp_sample_%s.txt")
PATIENT_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, "data_clinical_supp_patient_%s.txt"
)
MUTATIONS_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, "data_mutations_extended_%s.txt"
)
FUSIONS_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR, "data_fusions_%s.txt")
SEG_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR, "data_cna_hg19_%s.seg")
SV_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR, "data_sv_%s.txt")
BED_DIFFS_SEQASSAY_PATH = os.path.join(GENIE_RELEASE_DIR, "diff_%s.csv")


# TODO: Add to transform.py
def _to_redact_interval(df_col):
    """
    Determines interval values that are <18 and >89 that need to be redacted
    Returns bool because BIRTH_YEAR needs to be redacted as well based
    on the results

    Args:
        df_col: Dataframe column/pandas.Series of an interval column

    Returns:
        tuple: pandas.Series: to redact boolean vector
               pandas.Series: to redact pediatric boolean vector

    """
    phi_cutoff = 365 * 89
    pediatric_cutoff = 365 * 18
    # Some centers pre-redact their values by adding < or >. These
    # must be redacted
    contain_greaterthan = df_col.astype(str).str.contains(">", na=False)
    contain_lessthan = df_col.astype(str).str.contains("<", na=False)
    # Add in errors='coerce' to turn strings into NaN
    col_int = pd.to_numeric(df_col, errors="coerce")
    to_redact = (col_int > phi_cutoff) | contain_greaterthan
    to_redact_pediatric = (col_int < pediatric_cutoff) | contain_lessthan
    return to_redact, to_redact_pediatric


# TODO: Add to transform.py
def _redact_year(df_col):
    """Redacts year values that have < or >

    Args:
        df_col: Dataframe column/pandas.Series of a year column

    Returns:
        pandas.Series: Redacted series

    """
    year = df_col.astype(str)
    contain_greaterthan = year.str.contains(">", na=False)
    contain_lessthan = year.str.contains("<", na=False)
    df_col[contain_greaterthan] = "cannotReleaseHIPAA"
    df_col[contain_lessthan] = "withheld"
    return df_col


# TODO: Add to transform.py
def _to_redact_difference(df_col_year1, df_col_year2):
    """Determine if difference between year2 and year1 is > 89

    Args:
        df_col_year1: Dataframe column/pandas.Series of a year column
        df_col_year2: Dataframe column/pandas.Series of a year column

    Returns:
        pandas.Series: to redact boolean vector

    """
    # Add in errors='coerce' to turn strings into NaN
    year1 = pd.to_numeric(df_col_year1, errors="coerce")
    year2 = pd.to_numeric(df_col_year2, errors="coerce")
    to_redact = year2 - year1 > 89
    return to_redact


# TODO: Add to transform.py
def redact_phi(
    clinicaldf, interval_cols_to_redact=["AGE_AT_SEQ_REPORT", "INT_CONTACT", "INT_DOD"]
):
    """Redacts the PHI by re-annotating the clinical file

    Args:
        clinicaldf: merged clinical dataframe
        interval_cols_to_redact: List of interval columns to redact.
                                 Defaults to ['AGE_AT_SEQ_REPORT',
                                              'INT_CONTACT',
                                              'INT_DOD']

    Returns:
        pandas.DataFrame: Redacted clinical dataframe

    """
    # Moved to cannotReleaseHIPAA and withheld because the HIPAA
    # years would change every single year. (e.g. <1926, >1998 would be
    # inaccurate every year)
    for col in interval_cols_to_redact:
        to_redact, to_redactpeds = _to_redact_interval(clinicaldf[col])
        clinicaldf.loc[to_redact, "BIRTH_YEAR"] = "cannotReleaseHIPAA"
        clinicaldf.loc[to_redact, col] = ">32485"
        clinicaldf.loc[to_redactpeds, "BIRTH_YEAR"] = "withheld"
        clinicaldf.loc[to_redactpeds, col] = "<6570"
    # Redact BIRTH_YEAR values that have < or >
    # Birth year has to be done separately because it is not an interval
    clinicaldf["BIRTH_YEAR"] = _redact_year(clinicaldf["BIRTH_YEAR"])
    to_redact = _to_redact_difference(
        clinicaldf["BIRTH_YEAR"], clinicaldf["YEAR_CONTACT"]
    )
    clinicaldf.loc[to_redact, "BIRTH_YEAR"] = "cannotReleaseHIPAA"
    to_redact = _to_redact_difference(
        clinicaldf["BIRTH_YEAR"], clinicaldf["YEAR_DEATH"]
    )
    clinicaldf.loc[to_redact, "BIRTH_YEAR"] = "cannotReleaseHIPAA"

    return clinicaldf


# TODO: Add to transform.py
def remove_maf_samples(mafdf: pd.DataFrame, keep_samples: list) -> pd.DataFrame:
    """Remove samples from maf file

    Args:
        mafdf: Maf dataframe
        keep_samples: Samples to keep

    Returns:
        Filtered maf dataframe

    """
    keep_maf = mafdf["Tumor_Sample_Barcode"].isin(keep_samples)
    mafdf = mafdf.loc[keep_maf,]
    return mafdf


def get_whitelist_variants_idx(mafdf):
    """Get boolean vector for variants that are known somatic sites
    This is to override the germline filter
    """
    columns = ["Chromosome", "Start", "End", "Symbol"]
    whitelist = pd.read_csv(
        "https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.19/data/known_somatic_sites.bed",
        sep="\t",
        comment="#",
        header=None,
        names=columns,
    )
    rangesdf = mafdf[
        ["Chromosome", "Start_Position", "End_Position", "Hugo_Symbol", "HGVSp_Short"]
    ]
    rangesdf = rangesdf.rename(
        columns={"Start_Position": "Start", "End_Position": "End"}
    )
    maf_ranges = pyranges.PyRanges(rangesdf)
    whitelist_ranges = pyranges.PyRanges(whitelist)
    whitelisted_variants = maf_ranges.intersect(whitelist_ranges, how="containment")
    whitelist_variantsdf = whitelisted_variants.as_df()
    if not whitelist_variantsdf.empty:
        variants = (
            whitelist_variantsdf["Hugo_Symbol"]
            + " "
            + whitelist_variantsdf["HGVSp_Short"]
        )
    else:
        variants = []
    maf_variants = mafdf["Hugo_Symbol"] + " " + mafdf["HGVSp_Short"]
    # For some reason intersect and overlap doesn't work when
    # Start and End are the same. Here is an example that won't be
    # matched by the intersect function
    # variant: chr9-10-10
    # Bed: chr9-9-10
    match_start_end = mafdf["Start_Position"].isin(whitelist["Start"]) | mafdf[
        "End_Position"
    ].isin(whitelist["End"])
    return maf_variants.isin(variants) | match_start_end


# TODO: Add to transform.py
def configure_maf(mafdf, remove_variants, flagged_variants):
    """Configures each maf dataframe, does germline filtering

    Args:
        mafdf: Maf dataframe
        remove_variants: Variants to remove
        flagged_variants: Variants to flag

    Returns:
        configured maf row
    """
    variant = mafdf[
        [
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode",
        ]
    ].apply(lambda x: " ".join(x.map(str)), axis=1)
    mergecheck_variant = mafdf[
        [
            "Chromosome",
            "Start_Position",
            "HGVSp_Short",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode",
        ]
    ].apply(lambda x: " ".join(x.map(str)), axis=1)

    # Flag mutation in cis variants
    mafdf["mutationInCis_Flag"] = mergecheck_variant.isin(flagged_variants)
    # Remove common variants
    # na=False to resolve this linked error
    # https://stackoverflow.com/questions/52297740
    # common_variants = mafdf['FILTER'].astype(str).str.contains(
    #     "common_variant", na=False
    # )
    # Germline Filter
    gnomad_cols = [
        "gnomAD_AFR_AF",
        "gnomAD_AMR_AF",
        "gnomAD_ASJ_AF",
        "gnomAD_EAS_AF",
        "gnomAD_FIN_AF",
        "gnomAD_NFE_AF",
        "gnomAD_OTH_AF",
        "gnomAD_SAS_AF",
    ]
    # location of germline variants
    common_variants_idx = mafdf.loc[:, gnomad_cols].max(axis=1, skipna=True) > 0.0005

    # Remove specific variants
    to_remove_variants = variant.isin(remove_variants)
    # Genome Nexus successfully annotated (vcf2maf does not have this column)
    if mafdf.get("Annotation_Status") is None:
        mafdf["Annotation_Status"] = "SUCCESS"
    # Make sure to only get variants that were successfully annotated
    success = mafdf["Annotation_Status"] == "SUCCESS"
    # Get whitelisted variants
    whitelist_variants_idx = get_whitelist_variants_idx(mafdf)
    mafdf = mafdf.loc[
        (
            (~common_variants_idx | whitelist_variants_idx)
            & ~to_remove_variants
            & success
        ),
    ]
    # May not need to do this because these columns are always
    # returned as numerical values now
    # fillnas = ['t_depth', 't_ref_count', 't_alt_count',
    #            'n_depth', 'n_ref_count', 'n_alt_count']
    # for col in fillnas:
    #     mafdf[col][mafdf[col].astype(str) == "."] = float('nan')
    n_depth_ind = mafdf["n_depth"].astype(str).isin(["NA", "0.0", "0", "nan"])
    mafdf.loc[n_depth_ind, "Match_Norm_Seq_Allele2"] = ""
    mafdf.loc[n_depth_ind, "Match_Norm_Seq_Allele1"] = ""
    # Calculate missing t_depth, t_ref_count, t_alt_count
    t_counts = calculate_missing_variant_counts(
        depth=mafdf["t_depth"],
        alt_count=mafdf["t_alt_count"],
        ref_count=mafdf["t_ref_count"],
    )
    mafdf["t_depth"] = t_counts["depth"]
    mafdf["t_ref_count"] = t_counts["ref_count"]
    mafdf["t_alt_count"] = t_counts["alt_count"]
    # Calculate missing n_depth, n_ref_count, n_alt_count
    n_counts = calculate_missing_variant_counts(
        depth=mafdf["n_depth"],
        alt_count=mafdf["n_alt_count"],
        ref_count=mafdf["n_ref_count"],
    )
    mafdf["n_depth"] = n_counts["depth"]
    mafdf["n_ref_count"] = n_counts["ref_count"]
    mafdf["n_alt_count"] = n_counts["alt_count"]

    return mafdf


# TODO: Add to transform.py
def calculate_missing_variant_counts(
    depth: pd.Series, alt_count: pd.Series, ref_count: pd.Series
) -> dict:
    """Calculate missing counts. t_depth = t_alt_count + t_ref_count

    Args:
        depth: Allele Depth
        alt_count: Allele alt counts
        ref_count: Allele ref counts

    Returns:
        filled in depth, alt_count and ref_count values

    """
    # Avoid SettingWithCopyWarning
    depth = depth.copy()
    alt_count = alt_count.copy()
    ref_count = ref_count.copy()
    # t_depth = t_ref_count + t_alt_count
    null_depth = depth.isnull()
    # The notation null_depth_ref means all the reference values for which
    # depth is NA
    null_depth_ref = ref_count[null_depth]
    null_depth_alt = alt_count[null_depth]
    depth.loc[null_depth] = null_depth_ref + null_depth_alt
    # t_ref_count = t_depth - t_alt_count
    null_ref = ref_count.isnull()
    null_ref_depth = depth[null_ref]
    null_ref_alt = alt_count[null_ref]
    ref_count[null_ref] = null_ref_depth - null_ref_alt
    # t_alt_count = t_depth - t_ref_count
    null_alt = alt_count.isnull()
    null_alt_depth = depth[null_alt]
    null_alt_ref = ref_count[null_alt]
    alt_count[null_alt] = null_alt_depth - null_alt_ref
    return {"depth": depth, "ref_count": ref_count, "alt_count": alt_count}


# TODO: Add to transform.py
def runMAFinBED(
    syn,
    center_mappingdf,
    test=False,
    genieVersion="test",
    genie_user=None,
    genie_pass=None,
):
    """
    Run MAF in BED script, filter data and update MAFinBED database

    Args:
        syn: Synapse object
        center_mappingdf: center mapping dataframe
        test: Testing parameter. Default is False.
        genieVersion: GENIE version. Default is test.
        genie_user: Synapse username. Default is None.
        genie_pass: Synapse password.  Default is None.

    Returns:
        pd.Series: Variants to remove
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    mafinbed_script = os.path.join(script_dir, "../R/MAFinBED.R")
    # TODO: Use tempfile
    notinbed_file = os.path.join(script_dir, "../R/notinbed.csv")
    # The MAFinBED script filters out the centers that aren't being processed
    command = ["Rscript", mafinbed_script, notinbed_file]
    if genie_user is not None and genie_pass is not None:
        command.extend(["--syn_user", genie_user, "--syn_pass", genie_pass])
    if test:
        command.append("--testing")
    subprocess.check_call(command)

    # mutationSynId = databaseSynIdMappingDf['Id'][
    #     databaseSynIdMappingDf['Database'] == "vcf2maf"][0]
    # removedVariants = syn.tableQuery(
    #     "select Chromosome, Start_Position, End_Position, Reference_Allele, "
    #     "Tumor_Seq_Allele2, Tumor_Sample_Barcode, Center from {} where inBED"
    #     " is False and Center in ('{}')".format(
    #         mutationSynId, "','".join(center_mappingdf.center)))
    # removedVariantsDf = removedVariants.asDataFrame()
    removed_variantsdf = pd.read_csv(notinbed_file)
    removed_variantsdf["removeVariants"] = (
        removed_variantsdf["Chromosome"].astype(str)
        + " "
        + removed_variantsdf["Start_Position"].astype(str)
        + " "
        + removed_variantsdf["End_Position"].astype(str)
        + " "
        + removed_variantsdf["Reference_Allele"].astype(str)
        + " "
        + removed_variantsdf["Tumor_Seq_Allele2"].astype(str)
        + " "
        + removed_variantsdf["Tumor_Sample_Barcode"].astype(str)
    )
    # Store filtered variants
    for center in removed_variantsdf["Center"].unique():
        center_mutation = removed_variantsdf[removed_variantsdf["Center"] == center]
        # mafText = process.removePandasDfFloat(center_mutation)
        center_mutation.to_csv("mafinbed_filtered_variants.csv", index=False)
        load.store_file(
            syn=syn,
            filepath="mafinbed_filtered_variants.csv",
            parentid=center_mappingdf["stagingSynId"][
                center_mappingdf["center"] == center
            ][0],
            version_comment=genieVersion,
        )
        os.unlink("mafinbed_filtered_variants.csv")
    return removed_variantsdf["removeVariants"]


# TODO: Add to transform.py
def seq_date_filter(clinicalDf, processingDate, consortiumReleaseCutOff):
    """
    Filter samples by seq date

    Args:
        clinicalDf: Clinical dataframe
        processingDate: Processing date in form of Apr-XXXX
        consortiumReleaseCutOff: Release cut off days

    Returns:
        list: Samples to remove
    """
    removeSeqDateSamples = process_functions.seqDateFilter(
        clinicalDf, processingDate, consortiumReleaseCutOff
    )
    return removeSeqDateSamples


def sample_class_filter(clinical_df: pd.DataFrame) -> list:
    """Filter samples by SAMPLE_CLASS

    Args:
        clinical_df (pd.DataFrame): Clinical dataframe

    Returns:
        list: List of samples to filter out
    """
    if clinical_df.get("SAMPLE_CLASS") is not None:
        remove_samples = clinical_df["SAMPLE_ID"][
            clinical_df["SAMPLE_CLASS"] == "cfDNA"
        ].tolist()
    else:
        remove_samples = []
    return remove_samples


# TODO: Add to transform.py
def mutation_in_cis_filter(
    syn,
    skipMutationsInCis,
    variant_filtering_synId,
    center_mappingdf,
    genieVersion,
    test=False,
    genie_user=None,
    genie_pass=None,
):
    """
    Run mutation in cis filter, look up samples to remove

    Args:
        syn: Synapse object
        skipMutationsInCis: Skip this filter
        variant_filtering_synId: mergeCheck database dataframe
        center_mappingdf: center mapping dataframe
        genieVersion: GENIE version. Default is test.
        test: Testing parameter. Default is False.
        genie_user: Synapse username. Default is None.
        genie_pass: Synapse password.  Default is None.

    Returns:
        pd.Series: Samples to remove
    """
    if not skipMutationsInCis:
        mergeCheck_script = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "../R/mergeCheck.R"
        )
        command = ["Rscript", mergeCheck_script]
        # TODO: deprecate genie password soon
        if genie_user is not None and genie_pass is not None:
            command.extend(["--syn_user", genie_user, "--syn_pass", genie_pass])
        if test:
            command.append("--testing")
        # TODO: use subprocess.run instead
        subprocess.check_call(command)
        # Store each centers mutations in cis to their staging folder
        center_str = "','".join(center_mappingdf.center)
        query_str = (
            f"select * from {variant_filtering_synId} where Center in ('{center_str}')"
        )
        mergeCheckDf = extract.get_syntabledf(syn=syn, query_string=query_str)
        for center in mergeCheckDf.Center.unique():
            if not pd.isnull(center):
                stagingSynId = center_mappingdf.stagingSynId[
                    center_mappingdf["center"] == center
                ]
                mergeCheckDf[mergeCheckDf["Center"] == center].to_csv(
                    "mutationsInCis_filtered_samples.csv", index=False
                )
                load.store_file(
                    syn=syn,
                    filepath="mutationsInCis_filtered_samples.csv",
                    parentid=stagingSynId[0],
                    version_comment=genieVersion,
                )
                os.unlink("mutationsInCis_filtered_samples.csv")
    query_str = (
        f"SELECT Tumor_Sample_Barcode FROM {variant_filtering_synId} where "
        "Flag = 'TOSS' and Tumor_Sample_Barcode is not null"
    )
    filtered_samplesdf = extract.get_syntabledf(syn=syn, query_string=query_str)
    # #Alex script #1 removed patients
    remove_samples = filtered_samplesdf["Tumor_Sample_Barcode"].drop_duplicates()

    query_str = (
        f"SELECT * FROM {variant_filtering_synId} where "
        "Flag = 'Flag' and Tumor_Sample_Barcode is not null"
    )
    flag_variantsdf = extract.get_syntabledf(syn=syn, query_string=query_str)

    flag_variantsdf["flaggedVariants"] = (
        flag_variantsdf["Chromosome"].astype(str)
        + " "
        + flag_variantsdf["Start_Position"].astype(str)
        + " "
        + flag_variantsdf["HGVSp_Short"].astype(str)
        + " "
        + flag_variantsdf["Reference_Allele"].astype(str)
        + " "
        + flag_variantsdf["Tumor_Seq_Allele2"].astype(str)
        + " "
        + flag_variantsdf["Tumor_Sample_Barcode"].astype(str)
    )
    return (remove_samples, flag_variantsdf["flaggedVariants"])


# TODO: Add to transform.py
def seq_assay_id_filter(clinicaldf):
    """
    (Deprecated)
    Remove samples that are part of SEQ_ASSAY_IDs with less
    than 50 samples

    Args:
        clinicalDf: Sample clinical dataframe

    Returns:
        pd.Series: samples to remove
    """
    remove_seqassayid = clinicaldf["SEQ_ASSAY_ID"].value_counts()[
        clinicaldf["SEQ_ASSAY_ID"].value_counts() < 50
    ]
    clinicaldf = clinicaldf[
        clinicaldf["SEQ_ASSAY_ID"].isin(remove_seqassayid.keys().values)
    ]
    return clinicaldf["SAMPLE_ID"]


def no_genepanel_filter(clinicaldf, beddf):
    """
    Remove samples that don't have bed files associated with
    them

    Args:
        clinicaldf:  Clinical data
        beddf: bed data

    Returns:
        pd.Series: samples to remove
    """

    logger.info("NO GENE PANEL FILTER")
    has_seq_assay = clinicaldf["SEQ_ASSAY_ID"].isin(beddf["SEQ_ASSAY_ID"])
    remove_samples = clinicaldf["SAMPLE_ID"][~has_seq_assay]
    logger.info(
        "Removing samples with no bed file: {}".format(",".join(remove_samples))
    )
    return remove_samples


# TODO: add to load.py
def store_gene_panel_files(
    syn,
    fileviewSynId,
    genieVersion,
    data_gene_panel,
    consortiumReleaseSynId,
    wes_seqassayids,
):
    # Only need to upload these files once
    logger.info("STORING GENE PANELS FILES")
    wes_genepanel_filenames = [
        "data_gene_panel_{}.txt".format(seqassayid) for seqassayid in wes_seqassayids
    ]
    # Format string for tableQuery statement
    wes_genepanel_str = "','".join(wes_genepanel_filenames)
    # Only need to upload these files once
    logger.info("STORING GENE PANELS FILES")
    # This line of code is required to make sure any new files are
    # pulled into the file view.
    syn.tableQuery(f"select * from {fileviewSynId} limit 1")
    genePanelDf = extract.get_syntabledf(
        syn,
        f"select id from {fileviewSynId} where "
        "cBioFileFormat = 'genePanel' and "
        "fileStage = 'staging' and "
        f"name not in ('{wes_genepanel_str}')",
    )
    genePanelEntities = []
    panelNames = set(data_gene_panel["mutations"])
    print(f"EXISTING GENE PANELS: {','.join(panelNames)}")
    for synId in genePanelDf["id"]:
        genePanel = syn.get(synId)
        genePanelName = os.path.basename(genePanel.path)
        newGenePanelPath = os.path.join(GENIE_RELEASE_DIR, genePanelName)
        gene_panel = genePanelName.replace(".txt", "").replace("data_gene_panel_", "")
        print(gene_panel)
        if gene_panel in panelNames:
            os.rename(genePanel.path, newGenePanelPath)
            annotations = {"cBioFileFormat": "genePanel"}
            genePanelEntities.append(
                load.store_file(
                    syn=syn,
                    filepath=newGenePanelPath,
                    parentid=consortiumReleaseSynId,
                    version_comment=genieVersion,
                    name=genePanelName,
                    annotations=annotations,
                    used=f"{synId}.{genePanel.versionNumber}",
                )
            )
    return genePanelEntities


# TODO: add to load.py
def store_sv_files(
    syn: synapseclient.Synapse,
    release_synid: str,
    genie_version: str,
    synid: str,
    keep_for_center_consortium_samples: List[str],
    keep_for_merged_consortium_samples: List[str],
    current_release_staging: str,
    center_mappingdf: pd.DataFrame,
):
    """
    Create, filter, configure, and store structural variant file

    Args:
        syn: Synapse object
        release_synid: Synapse id to store release file
        genie_version: GENIE version (ie. v6.1-consortium)
        synid: SV database synid
        keep_for_center_consortium_samples: Samples to keep for center files
        keep_for_merged_consortium_samples: Samples to keep for merged file
        current_release_staging: Staging flag
        center_mappingdf: Center mapping dataframe
    """
    logger.info("MERING, FILTERING, STORING SV FILES")
    sv_df = extract.get_syntabledf(
        syn,
        f"select * from {synid}",
    )
    version = syn.create_snapshot_version(synid, comment=genie_version)

    # sv_df["ENTREZ_GENE_ID"].mask(
    #     sv_df["ENTREZ_GENE_ID"] == 0, float("nan"), inplace=True
    # )

    if not current_release_staging:
        sv_staging_df = sv_df[
            sv_df["SAMPLE_ID"].isin(keep_for_center_consortium_samples)
        ]
        for center in center_mappingdf.center:
            center_fusion = sv_staging_df[sv_staging_df["CENTER"] == center]
            if not center_fusion.empty:
                center_fusion.to_csv(SV_CENTER_PATH % center, sep="\t", index=False)
                load.store_file(
                    syn=syn,
                    filepath=SV_CENTER_PATH % center,
                    version_comment=genie_version,
                    parentid=center_mappingdf["stagingSynId"][
                        center_mappingdf["center"] == center
                    ][0],
                )

    sv_df = sv_df[sv_df["SAMPLE_ID"].isin(keep_for_merged_consortium_samples)]
    sv_df.rename(columns=transform._col_name_to_titlecase, inplace=True)
    sv_text = process_functions.removePandasDfFloat(sv_df)
    sv_path = os.path.join(GENIE_RELEASE_DIR, "data_sv.txt")
    with open(sv_path, "w") as sv_file:
        sv_file.write(sv_text)
    load.store_file(
        syn=syn,
        filepath=sv_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="data_sv.txt",
        used=f"{synid}.{version}",
    )


# TODO: Add to transform.py
def append_or_create_release_maf(dataframe: pd.DataFrame, filepath: str):
    """Creates a file with the dataframe or appends to a existing file.

    Args:
        df: pandas.dataframe to write out
        filepath: Filepath to append or create

    """
    if not os.path.exists(filepath) or os.stat(filepath).st_size == 0:
        data = process_functions.removePandasDfFloat(dataframe)
        with open(filepath, "w") as f_out:
            f_out.write(data)
    else:
        data = process_functions.removePandasDfFloat(dataframe, header=False)
        with open(filepath, "a") as f_out:
            f_out.write(data)


# TODO: Add to load.py
def store_maf_files(
    syn,
    genie_version,
    flatfiles_view_synid,
    release_synid,
    clinicaldf,
    center_mappingdf,
    keep_for_merged_consortium_samples,
    keep_for_center_consortium_samples,
    remove_mafinbed_variants,
    flagged_mutationInCis_variants,
    current_release_staging,
):
    """
    Create, filter, configure, and store maf file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        flatfiles_view_synid: Synapse id of fileview with all the flat files
        release_synid: Synapse id to store release file
        clinicaldf: Clinical dataframe with SAMPLE_ID and CENTER
        center_mappingdf: Center mapping dataframe
        keep_for_merged_consortium_samples: Samples to keep for merged file
        keep_for_center_consortium_samples: Samples to keep for center files
        remove_mafinbed_variants: Variants to remove
        flagged_mutationInCis_variants: Variants to flag
        current_release_staging: Staging flag
    """

    logger.info("FILTERING, STORING MUTATION FILES")
    centerMafSynIdsDf = extract.get_syntabledf(
        syn=syn,
        query_string=f"select id from {flatfiles_view_synid} where name like '%mutation%'",
    )
    mutations_path = os.path.join(GENIE_RELEASE_DIR, "data_mutations_extended.txt")
    with open(mutations_path, "w"):
        pass
    # Create maf file per center for their staging directory
    for center in clinicaldf["CENTER"].unique():
        with open(MUTATIONS_CENTER_PATH % center, "w"):
            pass
    used_entities = []
    # Must get the headers (because can't assume headers are the same order)
    maf_ent = syn.get(centerMafSynIdsDf.id[0])
    headerdf = pd.read_csv(maf_ent.path, sep="\t", comment="#", nrows=0)
    column_order = headerdf.columns
    for _, mafSynId in enumerate(centerMafSynIdsDf.id):
        maf_ent = syn.get(mafSynId)
        logger.info(maf_ent.path)
        # Extract center name
        center = maf_ent.path.split("_")[3].replace(".txt", "")
        if center in center_mappingdf.center.tolist():
            used_entities.append(f"{maf_ent.id}.{maf_ent.versionNumber}")
            mafchunks = pd.read_csv(
                maf_ent.path, sep="\t", comment="#", chunksize=100000
            )

            for mafchunk in mafchunks:
                # Reorder column headers
                mafchunk = mafchunk[column_order]
                # Get center for center staging maf
                # Configure maf
                configured_mafdf = configure_maf(
                    mafchunk, remove_mafinbed_variants, flagged_mutationInCis_variants
                )
                # Create maf for release
                merged_mafdf = remove_maf_samples(
                    configured_mafdf, keep_for_merged_consortium_samples
                )
                append_or_create_release_maf(merged_mafdf, mutations_path)
                # Create maf for center staging
                center_mafdf = remove_maf_samples(
                    configured_mafdf, keep_for_center_consortium_samples
                )
                append_or_create_release_maf(
                    center_mafdf, MUTATIONS_CENTER_PATH % center
                )

    load.store_file(
        syn=syn,
        filepath=mutations_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="data_mutations_extended.txt",
        used=used_entities,
    )

    if not current_release_staging:
        for center in clinicaldf["CENTER"].unique():
            staging_synid = center_mappingdf["stagingSynId"][
                center_mappingdf["center"] == center
            ][0]
            load.store_file(
                syn=syn,
                filepath=MUTATIONS_CENTER_PATH % center,
                version_comment=genie_version,
                parentid=staging_synid,
            )


# TODO: Add to transform.py
def run_genie_filters(
    syn,
    genie_user,
    genie_pass,
    genie_version,
    variant_filtering_synId,
    clinicaldf,
    beddf,
    center_mappingdf,
    processing_date,
    skip_mutationsincis,
    consortium_release_cutoff,
    test,
):
    """
    Run GENIE filters and returns variants and samples to remove

    Args:
        syn: Synapse object
        genie_user: Synapse username
        genie_pass: Synapse password
        genie_version: GENIE version (ie. v6.1-consortium)
        variant_filtering_synId: Synapse id of mutationInCis table
        clinicaldf: Clinical dataframe with SAMPLE_ID and SEQ_ASSAY_ID
        beddf: Bed dataframe
        center_mappingdf: Center mapping dataframe
        processing_date: Processing date
        skip_mutationsincis: Skip mutation in cis filter
        consortium_release_cutoff: Release cutoff in days
        test: Test flag

    Returns:
        pandas.Series: Variants to remove
        set: samples to remove for release files
        set: samples to remove for center files
        pandas.Series: Variants to flag
    """

    # ADD CHECKS TO CODE BEFORE UPLOAD.
    # Throw error if things don't go through
    logger.info("RUN GENIE FILTERS")
    # STORING CLINICAL FILES INTO CBIOPORTAL
    # FILTERING
    logger.info("MAF IN BED FILTER")
    remove_mafinbed_variants = runMAFinBED(
        syn,
        center_mappingdf,
        test=test,
        genieVersion=genie_version,
        genie_user=genie_user,
        genie_pass=genie_pass,
    )

    logger.info("MUTATION IN CIS FILTER")
    (
        remove_mutationincis_samples,
        flagged_mutationincis_variants,
    ) = mutation_in_cis_filter(
        syn,
        skip_mutationsincis,
        variant_filtering_synId,
        center_mappingdf,
        genieVersion=genie_version,
        test=test,
        genie_user=genie_user,
        genie_pass=genie_pass,
    )
    remove_no_genepanel_samples = no_genepanel_filter(clinicaldf, beddf)

    logger.info("SEQ DATE FILTER")
    remove_seqdate_samples = seq_date_filter(
        clinicaldf, processing_date, consortium_release_cutoff
    )

    # Only certain samples are removed for the files that go into
    # staging center folder
    remove_center_consortium_samples = set(remove_mutationincis_samples).union(
        set(remove_no_genepanel_samples)
    )
    # Most filteres are applied for the files that go into the merged
    # consortium release
    remove_merged_consortium_samples = set(remove_seqdate_samples)

    remove_merged_consortium_samples = remove_merged_consortium_samples.union(
        remove_center_consortium_samples
    )

    return (
        remove_mafinbed_variants,
        remove_merged_consortium_samples,
        remove_center_consortium_samples,
        flagged_mutationincis_variants,
    )


# TODO: Add to etl.py
def store_assay_info_files(
    syn, genie_version, assay_info_synid, clinicaldf, release_synid
):
    """Creates, stores assay information and gets WES panel list

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        assay_info_synid: Assay information database synid
        clinicaldf: Clinical dataframe with SAMPLE_ID and SEQ_ASSAY_ID
        release_synid: Synapse id to store release file

    Returns:
        List of whole exome sequencing SEQ_ASSAY_IDs
    """
    logger.info("Creates assay information file")
    assay_info_path = os.path.join(GENIE_RELEASE_DIR, "assay_information.txt")
    seq_assay_str = "','".join(clinicaldf["SEQ_ASSAY_ID"].unique())
    version = syn.create_snapshot_version(assay_info_synid, comment=genie_version)
    assay_infodf = extract.get_syntabledf(
        syn,
        f"select * from {assay_info_synid} where SEQ_ASSAY_ID "
        f"in ('{seq_assay_str}')",
    )
    assay_infodf.to_csv(assay_info_path, sep="\t", index=False)
    load.store_file(
        syn=syn,
        filepath=assay_info_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="assay_information.txt",
        used=f"{assay_info_synid}.{version}",
    )
    wes_index = assay_infodf["library_strategy"] == "WXS"
    wes_panels = assay_infodf["SEQ_ASSAY_ID"][wes_index]
    return wes_panels.tolist()


# TODO: Add to etl.py
def store_clinical_files(
    syn,
    genie_version,
    clinicaldf,
    oncotree_url,
    sample_cols,
    patient_cols,
    remove_center_consortium_samples,
    remove_merged_consortium_samples,
    release_synid,
    current_release_staging,
    center_mappingdf,
    used=None,
):
    """
    Create, filter, configure, and store clinical file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        clinicaldf: Clinical dataframe with SAMPLE_ID and SEQ_ASSAY_ID
        oncotree_url: Oncotree URL
        sample_cols: Clinical sample columns
        patient_cols: Clinical patient columns
        remove_center_consortium_samples: Samples to remove for center files
        remove_merged_consortium_samples: Samples to remove for merged file
        release_synid: Synapse id to store release file
        current_release_staging: Staging flag
        center_mappingdf: Center mapping dataframe

    Returns:
        pandas.DataFrame: configured clinical dataframe
        pandas.Series: samples to keep for center files
        pandas.Series: samples to keep for release files
    """

    logger.info("CONFIGURING CLINICAL FILES")
    logger.info("REMOVING PHI")
    # clinicaldf = redact_phi(clinicaldf)
    logger.info("ADD CANCER TYPES")
    # This removes support for both oncotree urls (only support json)
    oncotree_dict = process_functions.get_oncotree_code_mappings(oncotree_url)
    # Add in unknown key which maps to UNKNOWN everything
    oncotree_dict["UNKNOWN"] = {
        "CANCER_TYPE": "UNKNOWN",
        "CANCER_TYPE_DETAILED": "UNKNOWN",
        "ONCOTREE_PRIMARY_NODE": "UNKNOWN",
        "ONCOTREE_SECONDARY_NODE": "UNKNOWN",
    }

    clinicaldf["CANCER_TYPE"] = [
        oncotree_dict[code.upper()]["CANCER_TYPE"]
        if code.upper() in oncotree_dict.keys()
        else float("nan")
        for code in clinicaldf["ONCOTREE_CODE"]
    ]

    clinicaldf["CANCER_TYPE_DETAILED"] = [
        oncotree_dict[code.upper()]["CANCER_TYPE_DETAILED"]
        if code.upper() in oncotree_dict.keys()
        else float("nan")
        for code in clinicaldf["ONCOTREE_CODE"]
    ]

    clinicaldf["ONCOTREE_PRIMARY_NODE"] = [
        oncotree_dict[code.upper()]["ONCOTREE_PRIMARY_NODE"]
        if code.upper() in oncotree_dict.keys()
        else float("nan")
        for code in clinicaldf["ONCOTREE_CODE"]
    ]

    clinicaldf["ONCOTREE_SECONDARY_NODE"] = [
        oncotree_dict[code.upper()]["ONCOTREE_SECONDARY_NODE"]
        if code.upper() in oncotree_dict.keys()
        else float("nan")
        for code in clinicaldf["ONCOTREE_CODE"]
    ]

    # All cancer types that are null contain deprecated oncotree codes
    # And should be removed
    clinicaldf = clinicaldf[~clinicaldf["CANCER_TYPE"].isnull()]
    # Suggest using AGE_AT_SEQ_REPORT_DAYS instead so that the
    # descriptions can match
    clinicaldf["AGE_AT_SEQ_REPORT_DAYS"] = clinicaldf["AGE_AT_SEQ_REPORT"]
    clinicaldf["AGE_AT_SEQ_REPORT"] = [
        int(math.floor(int(float(age)) / 365.25))
        if process_functions.checkInt(age)
        else age
        for age in clinicaldf["AGE_AT_SEQ_REPORT"]
    ]
    clinicaldf["AGE_AT_SEQ_REPORT"][clinicaldf["AGE_AT_SEQ_REPORT"] == ">32485"] = ">89"
    clinicaldf["AGE_AT_SEQ_REPORT"][clinicaldf["AGE_AT_SEQ_REPORT"] == "<6570"] = "<18"

    ############################################################
    # CENTER SPECIFIC CODE FOR RIGHT NOW (REMOVE UHN-555-V1)
    ############################################################
    # clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'] != "UHN-555-V1"]
    # clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'] != "PHS-TRISEQ-V1"]

    # clinicalDf = clinicalDf[clinicalDf['CENTER'] != "WAKE"]
    # clinicalDf = clinicalDf[clinicalDf['CENTER'] != "CRUK"]
    ############################################################
    ############################################################

    clinicaldf.drop_duplicates("SAMPLE_ID", inplace=True)

    logger.info("STORING CLINICAL FILES")
    # samples must be removed after reading in the clinical file again
    staging_clinicaldf = clinicaldf[
        ~clinicaldf["SAMPLE_ID"].isin(remove_center_consortium_samples)
    ]
    if not current_release_staging:
        for center in center_mappingdf.center:
            center_clinical = staging_clinicaldf[staging_clinicaldf["CENTER"] == center]
            center_sample = center_clinical[sample_cols].drop_duplicates("SAMPLE_ID")
            center_patient = center_clinical[patient_cols].drop_duplicates("PATIENT_ID")
            center_sample.to_csv(SAMPLE_CENTER_PATH % center, sep="\t", index=False)
            center_patient.to_csv(PATIENT_CENTER_PATH % center, sep="\t", index=False)
            load.store_file(
                syn=syn,
                filepath=SAMPLE_CENTER_PATH % center,
                version_comment=genie_version,
                parentid=center_mappingdf["stagingSynId"][
                    center_mappingdf["center"] == center
                ][0],
            )
            load.store_file(
                syn=syn,
                filepath=PATIENT_CENTER_PATH % center,
                version_comment=genie_version,
                parentid=center_mappingdf["stagingSynId"][
                    center_mappingdf["center"] == center
                ][0],
            )

    clinicaldf = clinicaldf[
        ~clinicaldf["SAMPLE_ID"].isin(remove_merged_consortium_samples)
    ]

    keep_center_consortium_samples = staging_clinicaldf.SAMPLE_ID
    keep_merged_consortium_samples = clinicaldf.SAMPLE_ID
    # This mapping table is the GENIE clinical code to description
    # mapping to generate the headers of the clinical file
    mapping = extract.get_syntabledf(syn=syn, query_string="SELECT * FROM syn9621600")
    clinical_path = os.path.join(GENIE_RELEASE_DIR, "data_clinical.txt")
    clinical_sample_path = os.path.join(GENIE_RELEASE_DIR, "data_clinical_sample.txt")
    clinical_patient_path = os.path.join(GENIE_RELEASE_DIR, "data_clinical_patient.txt")
    process_functions.addClinicalHeaders(
        clinicaldf,
        mapping,
        patient_cols,
        sample_cols,
        clinical_sample_path,
        clinical_patient_path,
    )
    load.store_file(
        syn=syn,
        filepath=clinical_sample_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="data_clinical_sample.txt",
        used=used,
    )

    load.store_file(
        syn=syn,
        filepath=clinical_patient_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="data_clinical_patient.txt",
        used=used,
    )

    clinicaldf.to_csv(clinical_path, sep="\t", index=False)
    load.store_file(
        syn=syn,
        filepath=clinical_path,
        parentid=release_synid,
        name="data_clinical.txt",
        used=used,
        version_comment="database",
    )

    return (clinicaldf, keep_center_consortium_samples, keep_merged_consortium_samples)


# TODO: Add to etl.py
def store_cna_files(
    syn,
    flatfiles_view_synid,
    keep_for_center_consortium_samples,
    keep_for_merged_consortium_samples,
    center_mappingdf,
    genie_version,
    release_synid,
    current_release_staging,
):
    """
    Create, filter and store cna file

    Args:
        syn: Synapse object
        flatfiles_view_synid: Synapse id of fileview with all the flat files
        keep_for_center_consortium_samples: Samples to keep for center files
        keep_for_merged_consortium_samples: Samples to keep for merged file
        center_mappingdf: Center mapping dataframe
        genie_version: GENIE version (ie. v6.1-consortium)
        release_synid: Synapse id to store release file
    Returns:
        list: CNA samples
    """
    logger.info("MERING, FILTERING, STORING CNA FILES")
    cna_path = os.path.join(GENIE_RELEASE_DIR, "data_CNA.txt")
    query_str = ("select id from {} " "where name like 'data_CNA%'").format(
        flatfiles_view_synid
    )
    center_cna_synidsdf = extract.get_syntabledf(syn, query_str)
    # Grab all unique symbols and form cna_template
    all_symbols = set()
    for cna_synid in center_cna_synidsdf["id"]:
        cna_ent = syn.get(cna_synid)
        with open(cna_ent.path, "r") as cna_file:
            # Read first line first
            cna_file.readline()
            # Get all hugo symbols
            all_symbols = all_symbols.union(
                set(line.split("\t")[0] for line in cna_file)
            )
    cna_template = pd.DataFrame({"Hugo_Symbol": list(all_symbols)})
    cna_template.sort_values("Hugo_Symbol", inplace=True)
    cna_template.to_csv(cna_path, sep="\t", index=False)
    # Loop through to create finalized CNA file
    with_center_hugo_symbol = pd.Series("Hugo_Symbol")
    with_center_hugo_symbol = pd.concat(
        [with_center_hugo_symbol, pd.Series(keep_for_center_consortium_samples)]
    )

    with_merged_hugo_symbol = pd.Series("Hugo_Symbol")
    with_merged_hugo_symbol = pd.concat(
        [with_merged_hugo_symbol, pd.Series(keep_for_merged_consortium_samples)]
    )

    cna_samples = []
    used_entities = []
    for cna_synId in center_cna_synidsdf["id"]:
        cna_ent = syn.get(cna_synId)
        center = cna_ent.name.replace("data_CNA_", "").replace(".txt", "")
        logger.info(cna_ent.path)
        if center in center_mappingdf.center.tolist():
            used_entities.append(f"{cna_synId}.{cna_ent.versionNumber}")
            center_cna = pd.read_csv(cna_ent.path, sep="\t")
            merged_cna = cna_template.merge(center_cna, on="Hugo_Symbol", how="outer")
            merged_cna.sort_values("Hugo_Symbol", inplace=True)

            if not current_release_staging:
                merged_cna = merged_cna[
                    merged_cna.columns[merged_cna.columns.isin(with_center_hugo_symbol)]
                ]

                cna_text = process_functions.removePandasDfFloat(merged_cna)
                # Replace blank with NA's
                cna_text = cna_text.replace("\t\t", "\tNA\t")
                cna_text = cna_text.replace("\t\t", "\tNA\t")
                cna_text = cna_text.replace("\t\n", "\tNA\n")

                # Store center CNA file in staging dir
                with open(CNA_CENTER_PATH % center, "w") as cna_file:
                    cna_file.write(cna_text)
                load.store_file(
                    syn=syn,
                    filepath=CNA_CENTER_PATH % center,
                    version_comment=genie_version,
                    parentid=center_mappingdf["stagingSynId"][
                        center_mappingdf["center"] == center
                    ][0],
                )
            # This is to remove more samples for the final cna file
            merged_cna = merged_cna[
                merged_cna.columns[merged_cna.columns.isin(with_merged_hugo_symbol)]
            ]

            cna_text = process_functions.removePandasDfFloat(merged_cna)
            cna_text = cna_text.replace("\t\t", "\tNA\t")
            cna_text = cna_text.replace("\t\t", "\tNA\t")
            cna_text = cna_text.replace("\t\n", "\tNA\n")

            with open(CNA_CENTER_PATH % center, "w") as cna_file:
                cna_file.write(cna_text)
            # Join CNA file
            cna_samples.extend(merged_cna.columns[1:].tolist())
            linux_join_command = ["join", cna_path, CNA_CENTER_PATH % center]
            output = subprocess.check_output(linux_join_command)
            with open(cna_path, "w") as cna_file:
                cna_file.write(output.decode("utf-8").replace(" ", "\t"))

    load.store_file(
        syn=syn,
        filepath=cna_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="data_CNA.txt",
        used=used_entities,
    )

    return cna_samples


# TODO: Add to etl.py
def store_seg_files(
    syn,
    genie_version,
    seg_synid,
    release_synid,
    keep_for_center_consortium_samples,
    keep_for_merged_consortium_samples,
    center_mappingdf,
    current_release_staging,
):
    """
    Create, filter and store seg file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        seg_synid: Seg database synid
        release_synid: Synapse id to store release file
        keep_for_center_consortium_samples: Samples to keep for center files
        keep_for_merged_consortium_samples: Samples to keep for merged file
        center_mappingdf: Center mapping dataframe
        current_release_staging: Staging flag
    """
    logger.info("MERING, FILTERING, STORING SEG FILES")
    seg_path = os.path.join(GENIE_RELEASE_DIR, "data_cna_hg19.seg")
    version = syn.create_snapshot_version(seg_synid, comment=genie_version)

    segdf = extract.get_syntabledf(
        syn=syn,
        query_string=f"SELECT ID,CHROM,LOCSTART,LOCEND,NUMMARK,SEGMEAN,CENTER FROM {seg_synid}",
    )
    segdf = segdf.rename(
        columns={
            "CHROM": "chrom",
            "LOCSTART": "loc.start",
            "LOCEND": "loc.end",
            "SEGMEAN": "seg.mean",
            "NUMMARK": "num.mark",
        }
    )
    if not current_release_staging:
        staging_segdf = segdf[segdf["ID"].isin(keep_for_center_consortium_samples)]
        for center in center_mappingdf.center:
            center_seg = staging_segdf[staging_segdf["CENTER"] == center]
            if not center_seg.empty:
                del center_seg["CENTER"]
                segtext = process_functions.removePandasDfFloat(center_seg)
                with open(SEG_CENTER_PATH % center, "w") as seg_file:
                    seg_file.write(segtext)
                load.store_file(
                    syn=syn,
                    filepath=SEG_CENTER_PATH % center,
                    version_comment=genie_version,
                    parentid=center_mappingdf["stagingSynId"][
                        center_mappingdf["center"] == center
                    ][0],
                )
    del segdf["CENTER"]
    segdf = segdf[segdf["ID"].isin(keep_for_merged_consortium_samples)]
    segtext = process_functions.removePandasDfFloat(segdf)
    with open(seg_path, "w") as seg_file:
        seg_file.write(segtext)
    load.store_file(
        syn=syn,
        filepath=seg_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="data_cna_hg19.seg",
        used=f"{seg_synid}.{version}",
    )


# TODO: Add to etl.py
def store_data_gene_matrix(
    syn,
    genie_version,
    clinicaldf,
    cna_samples,
    release_synid,
    wes_seqassayids,
    used=None,
):
    """
    Create and store data gene matrix file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        clinicaldf: Clinical dataframe with SAMPLE_ID and SEQ_ASSAY_ID
        cna_samples: Samples with CNA
        release_synid: Synapse id to store release file

    Returns:
        pandas.DataFrame: data gene matrix dataframe
    """
    logger.info("STORING DATA GENE MATRIX FILE")
    data_gene_matrix_path = os.path.join(GENIE_RELEASE_DIR, "data_gene_matrix.txt")
    # Samples have already been removed
    data_gene_matrix = pd.DataFrame(columns=["SAMPLE_ID", "SEQ_ASSAY_ID"])
    data_gene_matrix = pd.concat(
        [data_gene_matrix, clinicaldf[["SAMPLE_ID", "SEQ_ASSAY_ID"]]]
    )
    data_gene_matrix = data_gene_matrix.rename(columns={"SEQ_ASSAY_ID": "mutations"})
    data_gene_matrix = data_gene_matrix[data_gene_matrix["SAMPLE_ID"] != ""]
    data_gene_matrix.drop_duplicates("SAMPLE_ID", inplace=True)
    # Gene panel file is written below CNA, because of the "cna" column
    # Add in CNA column into gene panel file
    cna_seqids = data_gene_matrix["mutations"][
        data_gene_matrix["SAMPLE_ID"].isin(cna_samples)
    ].unique()
    data_gene_matrix["cna"] = data_gene_matrix["mutations"]
    data_gene_matrix["cna"][~data_gene_matrix["cna"].isin(cna_seqids)] = "NA"
    wes_panel_mut = data_gene_matrix["mutations"].isin(wes_seqassayids)
    data_gene_matrix = data_gene_matrix[~wes_panel_mut]
    wes_panel_cna = data_gene_matrix["cna"].isin(wes_seqassayids)
    data_gene_matrix = data_gene_matrix[~wes_panel_cna]

    data_gene_matrix.to_csv(data_gene_matrix_path, sep="\t", index=False)

    load.store_file(
        syn=syn,
        filepath=data_gene_matrix_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="data_gene_matrix.txt",
    )
    return data_gene_matrix


# TODO: Add to etl.py
def store_bed_files(
    syn,
    genie_version,
    beddf,
    seq_assay_ids,
    center_mappingdf,
    current_release_staging,
    release_synid,
    used=None,
):
    """
    Store bed files, store the bed regions that had symbols remapped
    Filters bed file by clinical dataframe seq assays

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        beddf: Bed dataframe
        seq_assay_ids: All SEQ_ASSAY_IDs in the clinical file
        center_mappingdf: Center mapping dataframe
        current_release_staging: Staging flag
        release_synid: Synapse id to store release file
    """
    logger.info("STORING COMBINED BED FILE")
    combined_bed_path = os.path.join(GENIE_RELEASE_DIR, "genomic_information.txt")
    if not current_release_staging:
        for seq_assay in beddf["SEQ_ASSAY_ID"].unique():
            bed_seq_df = beddf[beddf["SEQ_ASSAY_ID"] == seq_assay]
            center = seq_assay.split("-")[0]
            bed_seq_df = bed_seq_df[bed_seq_df["Hugo_Symbol"] != bed_seq_df["ID"]]
            # There should always be a match here, because there should never
            # be a SEQ_ASSAY_ID that starts without the center name
            # If there is, check the bed db for SEQ_ASSAY_ID
            center_ind = center_mappingdf["center"] == center
            if not bed_seq_df.empty:
                bed_seq_df.to_csv(BED_DIFFS_SEQASSAY_PATH % seq_assay, index=False)
                load.store_file(
                    syn=syn,
                    filepath=BED_DIFFS_SEQASSAY_PATH % seq_assay,
                    version_comment=genie_version,
                    parentid=center_mappingdf["stagingSynId"][center_ind][0],
                )
    # This clinicalDf is already filtered through most of the filters
    beddf = beddf[beddf["SEQ_ASSAY_ID"].isin(seq_assay_ids)]
    beddf.to_csv(combined_bed_path, sep="\t", index=False)
    load.store_file(
        syn=syn,
        filepath=combined_bed_path,
        parentid=release_synid,
        version_comment=genie_version,
        name="genomic_information.txt",
        used=used,
    )


# TODO: Add to etl.py
def stagingToCbio(
    syn,
    processingDate,
    genieVersion,
    CENTER_MAPPING_DF,
    databaseSynIdMappingDf,
    oncotree_url=None,
    consortiumReleaseCutOff=183,
    current_release_staging=False,
    skipMutationsInCis=False,
    test=False,
    genie_user=None,
    genie_pass=None,
):
    """
    Main function that takes the GENIE database and creates release files

    Args:
        syn: Synapse object
        processingDate: Processing date in form of Apr-XXXX
        genieVersion: GENIE version. Default is test.
        CENTER_MAPPING_DF: center mapping dataframe
        databaseSynIdMappingDf: Database to Synapse Id mapping
        oncotree_url: Oncotree link
        consortiumReleaseCutOff: Release cut off days
        current_release_staging: Is it staging. Default is False.
        skipMutationsInCis: Skip mutation in cis filter. Default is False.
        test: Testing parameter. Default is False.
        genie_user: Synapse username. Default is None.
        genie_pass: Synapse password.  Default is None.

    Returns:
        list: Gene panel entities
    """
    if not os.path.exists(GENIE_RELEASE_DIR):
        os.mkdir(GENIE_RELEASE_DIR)
    consortiumReleaseSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "consortium"
    ][0]
    centerMafFileViewSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "centerMafView"
    ][0]
    patientSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "patient"
    ][0]
    sampleSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "sample"
    ][0]
    bedSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "bed"
    ][0]
    fileviewSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "fileview"
    ][0]
    segSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "seg"
    ][0]
    variant_filtering_synId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "mutationsInCis"
    ][0]
    fusionSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "fusions"
    ][0]
    sv_synid = databaseSynIdMappingDf["Id"][databaseSynIdMappingDf["Database"] == "sv"][
        0
    ]
    # Grab assay information
    assay_info_ind = databaseSynIdMappingDf["Database"] == "assayinfo"
    assay_info_synid = databaseSynIdMappingDf["Id"][assay_info_ind][0]

    # Using center mapping df to gate centers in release fileStage
    center_query_str = "','".join(CENTER_MAPPING_DF.center)
    patient_snapshot = syn.create_snapshot_version(patientSynId, comment=genieVersion)
    patient_used = f"{patientSynId}.{patient_snapshot}"
    patientDf = extract.get_syntabledf(
        syn, f"SELECT * FROM {patientSynId} where CENTER in ('{center_query_str}')"
    )
    sample_snapshot = syn.create_snapshot_version(sampleSynId, comment=genieVersion)
    sample_used = f"{sampleSynId}.{sample_snapshot}"
    sampleDf = extract.get_syntabledf(
        syn, f"SELECT * FROM {sampleSynId} where CENTER in ('{center_query_str}')"
    )
    bed_snapshot = syn.create_snapshot_version(bedSynId, comment=genieVersion)
    bed_used = f"{bedSynId}.{bed_snapshot}"
    bedDf = extract.get_syntabledf(
        syn,
        "SELECT Chromosome,Start_Position,End_Position,Hugo_Symbol,ID,"
        "SEQ_ASSAY_ID,Feature_Type,includeInPanel,clinicalReported FROM"
        f" {bedSynId} where CENTER in ('{center_query_str}')",
    )

    # Clinical release scope filter
    # If private -> Don't release to public
    clinicalReleaseScopeDf = extract.get_syntabledf(
        syn, "SELECT * FROM syn8545211 where releaseScope <> 'private'"
    )

    patientCols = clinicalReleaseScopeDf["fieldName"][
        clinicalReleaseScopeDf["level"] == "patient"
    ].tolist()
    sampleCols = clinicalReleaseScopeDf["fieldName"][
        clinicalReleaseScopeDf["level"] == "sample"
    ].tolist()

    # Remove this when these columns are removed from both databases
    if sampleDf.get("AGE_AT_SEQ_REPORT_NUMERICAL") is not None:
        del sampleDf["AGE_AT_SEQ_REPORT_NUMERICAL"]
    del sampleDf["CENTER"]
    # Remove this when these columns are removed from both databases
    if patientDf.get("BIRTH_YEAR_NUMERICAL") is not None:
        del patientDf["BIRTH_YEAR_NUMERICAL"]
    # del patientDf['BIRTH_YEAR_NUMERICAL']

    totalSample = ["PATIENT_ID"]
    totalSample.extend(sampleCols)
    sampleCols = totalSample
    # Make sure to only grab samples that have patient information
    sampleDf = sampleDf[sampleDf["PATIENT_ID"].isin(patientDf["PATIENT_ID"])]
    clinicalDf = sampleDf.merge(patientDf, on="PATIENT_ID", how="outer")
    # Remove patients without any sample or patient ids
    clinicalDf = clinicalDf[~clinicalDf["SAMPLE_ID"].isnull()]
    clinicalDf = clinicalDf[~clinicalDf["PATIENT_ID"].isnull()]

    (
        remove_mafInBed_variants,
        removeForMergedConsortiumSamples,
        removeForCenterConsortiumSamples,
        flagged_mutationInCis_variants,
    ) = run_genie_filters(
        syn,
        genie_user,
        genie_pass,
        genieVersion,
        variant_filtering_synId,
        clinicalDf,
        bedDf,
        CENTER_MAPPING_DF,
        processingDate,
        skipMutationsInCis,
        consortiumReleaseCutOff,
        test,
    )

    (
        clinicalDf,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
    ) = store_clinical_files(
        syn,
        genieVersion,
        clinicalDf,
        oncotree_url,
        sampleCols,
        patientCols,
        removeForCenterConsortiumSamples,
        removeForMergedConsortiumSamples,
        consortiumReleaseSynId,
        current_release_staging,
        CENTER_MAPPING_DF,
        used=[sample_used, patient_used],
    )

    assert not clinicalDf["SAMPLE_ID"].duplicated().any()

    store_maf_files(
        syn,
        genieVersion,
        centerMafFileViewSynId,
        consortiumReleaseSynId,
        clinicalDf[["SAMPLE_ID", "CENTER"]],
        CENTER_MAPPING_DF,
        keepForMergedConsortiumSamples,
        keepForCenterConsortiumSamples,
        remove_mafInBed_variants,
        flagged_mutationInCis_variants,
        current_release_staging,
    )

    cnaSamples = store_cna_files(
        syn,
        centerMafFileViewSynId,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        CENTER_MAPPING_DF,
        genieVersion,
        consortiumReleaseSynId,
        current_release_staging,
    )

    wes_panelids = store_assay_info_files(
        syn, genieVersion, assay_info_synid, clinicalDf, consortiumReleaseSynId
    )

    data_gene_matrix = store_data_gene_matrix(
        syn, genieVersion, clinicalDf, cnaSamples, consortiumReleaseSynId, wes_panelids
    )

    genePanelEntities = store_gene_panel_files(
        syn,
        fileviewSynId,
        genieVersion,
        data_gene_matrix,
        consortiumReleaseSynId,
        wes_panelids,
    )

    store_sv_files(
        syn,
        consortiumReleaseSynId,
        genieVersion,
        sv_synid,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        current_release_staging,
        CENTER_MAPPING_DF,
    )

    store_seg_files(
        syn,
        genieVersion,
        segSynId,
        consortiumReleaseSynId,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        CENTER_MAPPING_DF,
        current_release_staging,
    )

    store_bed_files(
        syn,
        genieVersion,
        bedDf,
        clinicalDf["SEQ_ASSAY_ID"].unique(),
        CENTER_MAPPING_DF,
        current_release_staging,
        consortiumReleaseSynId,
        used=bed_used,
    )

    return genePanelEntities


# TODO: Add to transform.py
def revise_metadata_files(syn, consortiumid, genie_version=None):
    """
    Rewrite metadata files with the correct GENIE version

    Args:
        syn: Synapse object
        consortiumid: Synapse id of consortium release folder
        genie_version: GENIE version, Default to None
    """
    release_files = syn.getChildren(consortiumid)
    meta_file_ents = [
        syn.get(
            i["id"], downloadLocation=GENIE_RELEASE_DIR, ifcollision="overwrite.local"
        )
        for i in release_files
        if "meta" in i["name"] and i["name"] != "meta_fusions.txt"
    ]

    for meta_ent in meta_file_ents:
        with open(meta_ent.path, "r+") as meta:
            meta_text = meta.read()
            if "meta_study" not in meta_ent.path:
                version = ""
            else:
                version = re.search(".+GENIE.+v(.+)", meta_text).group(1)
            # Fix this line
            genie_version = version if genie_version is None else genie_version

            if version != genie_version:
                meta_text = meta_text.replace(
                    "GENIE Cohort v{}".format(version),
                    "GENIE Cohort v{}".format(genie_version),
                )

                meta_text = meta_text.replace(
                    "GENIE v{}".format(version), "GENIE v{}".format(genie_version)
                )

                meta.seek(0)
                meta.write(meta_text)
                meta.truncate()
        load.store_file(
            syn=syn,
            filepath=meta_ent.path,
            parentid=consortiumid,
            version_comment=genie_version,
        )


# TODO: Add to load.py
def search_or_create_folder(
    syn: synapseclient.Synapse, parentid: str, folder_name: str
) -> str:
    """
    Searches for an existing Synapse Folder given a parent id
    and creates the Synapse folder if it doesn't exist

    Args:
        syn (synapseclient.Synapse): Synapse connection
        parentid (str): Synapse Id of a project or folder
        folder_name (str): Folde rname

    Returns:
        str: Synapse Folder id
    """
    folder_id = syn.findEntityId(name=folder_name, parent=parentid)
    # if case_lists doesn't exist
    if folder_id is None:
        folder_ent = synapseclient.Folder(name=folder_name, parent=parentid)
        folder_id = syn.store(folder_ent).id
    return folder_id


# TODO: Add to load.py
def create_link_version(
    syn,
    genie_version,
    case_list_entities,
    gene_panel_entities,
    database_synid_mappingdf,
    release_type="consortium",
):
    """
    Create release links from the actual entity and version

    TODO: Refactor to use fileviews

    Args:
        syn: Synapse object
        genie_version: GENIE version number
        case_list_entities: Case list entities
        gene_panel_entities: Gene panel entities
        database_synid_mappingdf: dataframe containing database to
                                  synapse id mapping
        release_type: 'consortium' or 'public' release
    """
    # Grab major release numbers (ie 1,2,3 ...)
    major_release = genie_version.split(".")[0]
    all_releases_synid = database_synid_mappingdf["Id"][
        database_synid_mappingdf["Database"] == "release"
    ].values[0]
    # Create major release folder
    major_release_folder_synid = search_or_create_folder(
        syn, all_releases_synid, f"Release {major_release}"
    )
    # If the major release folder didn't exist, go ahead and create the
    # release folder
    release_folder_synid = search_or_create_folder(
        syn, major_release_folder_synid, genie_version
    )
    # Search or create case lists folder
    caselist_folder_synid = search_or_create_folder(
        syn, release_folder_synid, "case_lists"
    )

    # caselistId = findCaseListId(syn, release_folder_synid)
    consortium_synid = database_synid_mappingdf["Id"][
        database_synid_mappingdf["Database"] == release_type
    ].values[0]
    consortium_release_files = syn.getChildren(consortium_synid)

    for release_file in consortium_release_files:
        not_folder = release_file["type"] != "org.sagebionetworks.repo.model.Folder"
        # data_clinical.txt MUST be pulled in when doing consortium release
        not_public = (
            release_file["name"] != "data_clinical.txt" or release_type == "consortium"
        )
        is_gene_panel = release_file["name"].startswith("data_gene_panel")
        is_depreciated_file = release_file["name"] in ["data_fusions.txt"]

        if not_folder and not_public and not is_gene_panel and not is_depreciated_file:
            syn.store(
                synapseclient.Link(
                    release_file["id"],
                    parent=release_folder_synid,
                    targetVersion=release_file["versionNumber"],
                )
            )

    release_files = syn.getChildren(release_folder_synid)
    clinical_ent = [
        ents["id"] for ents in release_files if ents["name"] == "data_clinical.txt"
    ]
    if clinical_ent:
        # Set private permission for the data_clinical.txt link
        syn.setPermissions(clinical_ent[0], principalId=3346558, accessType=[])
        syn.setPermissions(clinical_ent[0], principalId=3326313, accessType=[])

    for ents in case_list_entities:
        syn.store(
            synapseclient.Link(
                ents.id, parent=caselist_folder_synid, targetVersion=ents.versionNumber
            )
        )

    # Store gene panels
    for ents in gene_panel_entities:
        syn.store(
            synapseclient.Link(
                ents.id, parent=release_folder_synid, targetVersion=ents.versionNumber
            )
        )

    return {
        "release_folder": release_folder_synid,
        "caselist_folder": caselist_folder_synid,
    }
