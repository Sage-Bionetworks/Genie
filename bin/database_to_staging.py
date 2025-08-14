"""
Workflow to trigger a consortium release

``` mermaid
flowchart TD
    A["Start consortium release"] --> B["Remove old files in GENIE release dir"]
    %% B --> C["Login to Synapse"]
    %% C --> D{"Test / Staging / Prod Environment?"}
    %% D -->|Test| E["Set test databaseSynIdMappingId"]
    %% D -->|Staging| F["Set staging databaseSynIdMappingId, skip mutations in cis"]
    %% D -->|Production| G["Set production databaseSynIdMappingId"]
    B --> H["Prepare for processing by getting databaseSynIdMapping table, oncotree link, cbioportal folder, database synapse ids, etc"]
    %% H --> I{"Oncotree link provided?"}
    %% I -->|No| J["Extract Oncotree URL from database"]
    %% I -->|Yes| K["Use provided Oncotree URL"]
    %% J & K --> L["Check Oncotree URL accessibility"]
    %% L --> M["Validate cBioPortal path exists"]
    %% M --> N["Get Synapse IDs for consortium, process tracker, etc."]
    %% N --> O["Create or retrieve case_lists folder"]
    %% O --> P{"Staging?"}
    %% P -- No --> Q["Start process tracking"]
    H --> R["Query for all the centers for which data should be released"]

    %% Expanded stagingToCbio logic
    %% R --> S2["Create GENIE release directory if missing"]
    %% S1 --> S2["Extract Synapse Table IDs (patient, sample, maf, bed, seg, etc.)"]
    R --> S3["Create snapshots & pull patient, sample, bed data"]
    S3 --> S4["Merge patient + sample tables into clinicalDf"]
    S4 --> S5["Run GENIE filters"]
    S5 --> SD1["Variant Filters"]
    SD1 --> SD2["Germline filter"]
    SD1 --> SD3["MAF in BED"]
    SD2 & SD3 --> SD5["List of variants to remove"]
    S5 --> SE1["Sample Filters"]
    SE1 --> SE2["SEQ_DATE"]
    SE1 --> SE3["No Bed file"]
    SE1 --> SE4["Oncotree"]
    SE1 --> SE5["Mutation In Cis"]
    SE2 & SE3 & SE4 & SE5 --> SE6["List of samples to remove"]
    SE6 --> SS3["Merge/Filter/store clinical file"]
    %% SS3 --> S6["Merge/Filter/store clinical file"]
    SD5 & SS3 --> S7["Merge/Filter/store MAF file"]
    SS3 --> S8["Merge/Filter/store CNA file"]
    SS3 --> S9["Merge/Filter/store Assay information file"]
    SS3 --> S10["Merge/Filter/store SV file"]
    S8 & S7 & S9 & S10 --> S11["Merge/Filter/store Data Gene Matrix"]
    S11 --> S12["Download and upload gene panel files"]
    SS3 --> S13["Merge/Filter/store SEG file"]
    SS3 --> S14["Merge/Filter/store BED files"]
    %% SS3 --> S15["Return list of gene panel entities"]

    SS3 & S12 & S13 & S14 --> T["Remove old case list files"]
    T --> U["Generate new case lists and upload to Synapse"]
    U --> V["Revise metadata files with correct GENIE version"]
    V --> W["Run cBioPortal validation script"]
    %% W --> X{"Production?"}
    %% X -->|Yes| Y["Upload validation logs to Synapse"]
    W --> Z["Create release folder with links to files"]
    %% Z --> AA{"Production?"}
    %% AA -->|Yes| AB["End process tracking"]
    Z --> AC["Run dashboard updater"]
    AC --> AD["Generate dashboard HTML"]
    AD --> AF["End"]
```

"""
import argparse
import datetime
import logging
import os
import subprocess

import synapseclient
from genie import (
    create_case_lists,
    dashboard_table_updater,
    database_to_staging,
    load,
    process_functions,
)

# import time


logger = logging.getLogger(__name__)

PWD = os.path.dirname(os.path.abspath(__file__))


def generate_dashboard_html(
    genie_version: str, staging: bool = False, testing: bool = False
):
    """Generates dashboard html writeout that gets uploaded to the
    release folder

    Args:
        genie_version: GENIE release
        staging: Use staging files. Default is False
        testing: Use testing files. Default is False

    """
    markdown_render_cmd = [
        "Rscript",
        os.path.join(PWD, "../R/dashboard_markdown_generator.R"),
        genie_version,
        "--template_path",
        os.path.join(PWD, "../templates/dashboardTemplate.Rmd"),
    ]

    if staging:
        markdown_render_cmd.append("--staging")
    if testing:
        markdown_render_cmd.append("--testing")
    subprocess.check_call(markdown_render_cmd)


def generate_data_guide(genie_version, oncotree_version=None, database_mapping=None):
    """Generates the GENIE data guide"""

    template_path = os.path.join(PWD, "../templates/data_guide_template.Rnw")
    with open(template_path, "r") as template_file:
        template_str = template_file.read()

    replacements = {
        "{{release}}": genie_version,
        "{{database_synid}}": database_mapping,
        "{{oncotree}}": oncotree_version.replace("_", "\\_"),
        "{{genie_banner}}": os.path.join(PWD, "../genie_banner.png"),
    }

    for search in replacements:
        replacement = replacements[search]
        # If no replacement value is passed in, don't replace
        if replacement is not None:
            template_str = template_str.replace(search, replacement)

    with open(os.path.join(PWD, "data_guide.Rnw"), "w") as data_guide_file:
        data_guide_file.write(template_str)

    subprocess.check_call(
        ["R", "CMD", "Sweave", "--pdf", os.path.join(PWD, "data_guide.Rnw")]
    )
    return "data_guide.pdf"


def main(
    genie_version,
    processing_date,
    cbioportal_path,
    oncotree_link=None,
    consortium_release_cutoff=184,
    test=False,
    staging=False,
    debug=False,
    skip_mutationsincis=False,
):
    """
    - Does parameter checks
    - Updates process tracking start
    - initiates database to staging
    - create case lists
    - revise meta files
    - run cBioPortal validation
    - create link versions
    - update process tracking end
    - Create dashboard tables and plots

    Args:
        genie_version: GENIE version,
        processing_date: processing date
        cbioportal_path: Path to cbioportal validator
        oncotree_link: Link to oncotree codes
        consortium_release_cutoff: release cut off value in days
        test: Test flag, uses test databases
        staging: Staging flag, uses staging databases
        debug:  Synapse debug flag
        skip_mutationsincis: Skip mutation in cis filter
    """
    # HACK: Delete all existing files first
    process_functions.rmFiles(database_to_staging.GENIE_RELEASE_DIR)

    syn = process_functions.synapse_login(debug=debug)
    # HACK: Use project id instead of this...
    if test:
        databaseSynIdMappingId = "syn11600968"
        genie_version = "TESTING"
    elif staging:
        databaseSynIdMappingId = "syn12094210"
    else:
        databaseSynIdMappingId = "syn10967259"
    # Database/folder syn id mapping
    databaseSynIdMapping = syn.tableQuery(
        "select * from {}".format(databaseSynIdMappingId)
    )
    databaseSynIdMappingDf = databaseSynIdMapping.asDataFrame()
    # databaseSynIdMappingDf.index = databaseSynIdMappingDf.Database
    # del databaseSynIdMappingDf['Database']
    # databaseSynIdMappingDf.to_dict()

    if oncotree_link is None:
        oncoLink = databaseSynIdMappingDf["Id"][
            databaseSynIdMappingDf["Database"] == "oncotreeLink"
        ].values[0]
        oncoLinkEnt = syn.get(oncoLink)
        oncotree_link = oncoLinkEnt.externalURL

    # Check if you can connect to oncotree link,
    # if not then don't run validation / processing
    process_functions.checkUrl(oncotree_link)

    cbioValidatorPath = os.path.join(
        cbioportal_path, "core/src/main/scripts/importer/validateData.py"
    )
    assert os.path.exists(cbioValidatorPath), "Please specify correct cbioportalPath"
    syn.table_query_timeout = 50000

    consortiumSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "consortium"
    ].values[0]
    processTrackerSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "processTracker"
    ].values[0]
    # get syn id of case list folder in consortium release
    # caseListSynId = findCaseListId(syn, consortiumSynId)
    caseListSynId = database_to_staging.search_or_create_folder(
        syn, consortiumSynId, "case_lists"
    )

    if not staging:
        load.update_process_trackingdf(
            syn=syn,
            process_trackerdb_synid=processTrackerSynId,
            center="SAGE",
            process_type="dbToStage",
            start=True,
        )

    centerMappingSynId = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "centerMapping"
    ].values[0]
    # Only release files where release is true
    center_mapping = syn.tableQuery(
        "SELECT * FROM {} where release is true".format(centerMappingSynId)
    )
    center_mappingdf = center_mapping.asDataFrame()
    processingDate = datetime.datetime.strptime(processing_date, "%b-%Y")

    logger.info("STAGING TO CONSORTIUM")
    genePanelEntities = database_to_staging.stagingToCbio(
        syn,
        processingDate,
        genie_version,
        center_mappingdf,
        databaseSynIdMappingDf,
        oncotree_url=oncotree_link,
        consortiumReleaseCutOff=consortium_release_cutoff,
        current_release_staging=staging,
        skipMutationsInCis=skip_mutationsincis,
        test=test,
    )

    # Create case lists files
    logger.info("CREATE CASE LIST FILES")
    # Remove old caselists first
    if not os.path.exists(database_to_staging.CASE_LIST_PATH):
        os.mkdir(database_to_staging.CASE_LIST_PATH)
    caselists = os.listdir(database_to_staging.CASE_LIST_PATH)
    for caselist in caselists:
        os.remove(os.path.join(database_to_staging.CASE_LIST_PATH, caselist))
    clinical_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "data_clinical.txt",
    )
    assay_information_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "assay_information.txt",
    )
    create_case_lists.main(
        clinical_path,
        assay_information_path,
        database_to_staging.CASE_LIST_PATH,
        "genie_private",
    )
    caseListFiles = os.listdir(database_to_staging.CASE_LIST_PATH)
    caseListEntities = []
    for casePath in caseListFiles:
        casePath = os.path.join(database_to_staging.CASE_LIST_PATH, casePath)
        caseListEntities.append(
            load.store_file(
                syn=syn,
                filepath=casePath,
                parentid=caseListSynId,
                version_comment=genie_version,
            )
        )

    logger.info("REMOVING UNNECESSARY FILES")
    # genie_files = os.listdir(database_to_staging.GENIE_RELEASE_DIR)
    # for genie_file in genie_files:
    #     if (
    #         genie_version not in genie_file
    #         and "meta" not in genie_file
    #         and "case_lists" not in genie_file
    #     ):
    #         os.remove(os.path.join(database_to_staging.GENIE_RELEASE_DIR, genie_file))
    # os.remove(clinical_path)

    logger.info("REVISE METADATA FILES")
    database_to_staging.revise_metadata_files(syn, consortiumSynId, genie_version)

    logger.info("CBIO VALIDATION")

    # Must be exit 0 because the validator sometimes fails,
    # but we still want to capture the output

    command = [
        cbioValidatorPath,
        "-s",
        database_to_staging.GENIE_RELEASE_DIR,
        "-n",
        "; exit 0",
    ]
    cbioOutput = subprocess.check_output(" ".join(command), shell=True)
    logger.info(cbioOutput.decode("utf-8"))

    cbio_validator_log = f"cbioValidatorLogsConsortium_{genie_version}.txt"
    if not test and not staging:
        log_folder_synid = databaseSynIdMappingDf["Id"][
            databaseSynIdMappingDf["Database"] == "logs"
        ].values[0]
        with open(cbio_validator_log, "w") as cbio_log:
            cbio_log.write(cbioOutput.decode("utf-8"))
        syn.store(synapseclient.File(cbio_validator_log, parentId=log_folder_synid))
        os.remove(cbio_validator_log)
    # HACK: Instead of doing this, files should be written to a tempdir...
    # logger.info("REMOVING OLD FILES")

    # process_functions.rmFiles(database_to_staging.CASE_LIST_PATH)
    # private_cna_meta_path = os.path.join(
    #     database_to_staging.GENIE_RELEASE_DIR, "genie_private_meta_cna_hg19_seg.txt"
    # )
    # if os.path.exists(private_cna_meta_path):
    #     os.unlink(private_cna_meta_path)

    logger.info("CREATING LINK VERSION")
    # Returns release and case list folder
    _ = database_to_staging.create_link_version(
        syn, genie_version, caseListEntities, genePanelEntities, databaseSynIdMappingDf
    )

    if not staging:
        load.update_process_trackingdf(
            syn=syn,
            process_trackerdb_synid=processTrackerSynId,
            center="SAGE",
            process_type="dbToStage",
            start=False,
        )

    logger.info("DASHBOARD UPDATE")
    # Only run dashboard update if not testing or staging
    if not args.test and not args.staging:
        dashboard_table_updater.run_dashboard(
            syn, databaseSynIdMappingDf, genie_version, staging=staging
        )
    generate_dashboard_html(genie_version, staging=staging, testing=test)
    logger.info("DASHBOARD UPDATE COMPLETE")
    logger.info("AUTO GENERATE DATA GUIDE")

    # TODO: remove data guide code
    # oncotree_version = oncotree_link.split("=")[1]
    # data_guide_pdf = generate_data_guide(
    #    genie_version,
    #    oncotree_version=oncotree_version,
    #    database_mapping=databaseSynIdMappingId,
    # )
    # load.store_file(
    #    syn=syn,
    #    filepath=data_guide_pdf,
    #    version_comment=genie_version,
    #    parentid=folders["release_folder"],
    # )
    logger.info("COMPLETED DATABASE TO STAGING")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Release GENIE consortium files")

    parser.add_argument(
        "processingDate",
        type=str,
        metavar="Jan-2017",
        help="The processing date of GENIE in Month-Year format" " (ie. Apr-2017)",
    )

    parser.add_argument(
        "cbioportalPath",
        type=str,
        metavar="/path/to/cbioportal",
        help="Make sure you clone the cbioportal github: "
        "git clone https://github.com/cBioPortal/cbioportal.git",
    )

    parser.add_argument(
        "genieVersion", type=str, metavar="1.0.1", help="GENIE release version"
    )

    parser.add_argument("--oncotree_link", type=str, help="Link to oncotree code")

    parser.add_argument(
        "--consortiumReleaseCutOff",
        type=int,
        metavar="184",
        default=184,
        help="Consortium release cut off time in days",
    )

    parser.add_argument(
        "--skipMutationsInCis",
        action="store_true",
        help="Skip running mutation in cis script",
    )

    parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")

    parser.add_argument("--debug", action="store_true", help="Synapse debug feature")

    test_group = parser.add_mutually_exclusive_group()
    test_group.add_argument("--test", action="store_true", help="Run test")

    test_group.add_argument(
        "--staging", action="store_true", help="Store into staging folder"
    )

    args = parser.parse_args()

    main(
        genie_version=args.genieVersion,
        processing_date=args.processingDate,
        cbioportal_path=args.cbioportalPath,
        oncotree_link=args.oncotree_link,
        consortium_release_cutoff=args.consortiumReleaseCutOff,
        test=args.test,
        staging=args.staging,
        debug=args.debug,
        skip_mutationsincis=args.skipMutationsInCis,
    )
