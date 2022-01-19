import argparse
import datetime
import logging
import os
import subprocess
import synapseclient

# import time

from genie import (
    create_case_lists,
    dashboard_table_updater,
    database_to_staging,
    process_functions,
)

logger = logging.getLogger(__name__)

PWD = os.path.dirname(os.path.abspath(__file__))


def generate_dashboard_html(
    genie_version, staging=False, genie_user=None, genie_pass=None
):
    """Generates dashboard html writeout that gets uploaded to the
    release folder

    Args:
        syn: Synapse connection
        genie_version: GENIE release
        staging: Use staging files. Default is False
        genie_user: GENIE synapse username
        genie_pass: GENIE synapse password

    """
    markdown_render_cmd = [
        "Rscript",
        os.path.join(PWD, "../R/dashboard_markdown_generator.R"),
        genie_version,
        "--template_path",
        os.path.join(PWD, "../templates/dashboardTemplate.Rmd"),
    ]

    if genie_user is not None and genie_pass is not None:
        markdown_render_cmd.extend(["--syn_user", genie_user, "--syn_pass", genie_pass])
    if staging:
        markdown_render_cmd.append("--staging")
    subprocess.check_call(markdown_render_cmd)


def generate_data_guide(
    genie_version,
    oncotree_version=None,
    database_mapping=None,
    genie_user=None,
    genie_pass=None,
):
    """Generates the GENIE data guide"""

    template_path = os.path.join(PWD, "../templates/data_guide_template.Rnw")
    with open(template_path, "r") as template_file:
        template_str = template_file.read()

    replacements = {
        "{{release}}": genie_version,
        "{{database_synid}}": database_mapping,
        "{{oncotree}}": oncotree_version.replace("_", "\\_"),
        "{{username}}": genie_user,
        "{{password}}": genie_pass,
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
    pemfile=None,
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
        pemfile: Path to private key file
        test: Test flag, uses test databases
        staging: Staging flag, uses staging databases
        debug:  Synapse debug flag
        skip_mutationsincis: Skip mutation in cis filter
    """
    syn = process_functions.synLogin(pemfile, debug=debug)
    genie_user = os.environ.get("GENIE_USER")
    if pemfile is not None:
        genie_pass = process_functions.get_password(pemfile)
    else:
        genie_pass = None

    if test:
        databaseSynIdMappingId = "syn11600968"
        genie_version = "TESTING"
    elif staging:
        skip_mutationsincis = True
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
    caseListSynId, _ = database_to_staging.search_and_create_folder(
        syn, consortiumSynId, "case_lists"
    )

    if not staging:
        database_to_staging.update_process_trackingdf(
            syn, processTrackerSynId, "SAGE", "dbToStage", start=True
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
        genie_user=genie_user,
        genie_pass=genie_pass,
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
        "data_clinical_{}.txt".format(genie_version),
    )
    assay_information_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "assay_information_{}.txt".format(genie_version),
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
            database_to_staging.store_file(
                syn, casePath, parent=caseListSynId, genieVersion=genie_version
            )
        )

    logger.info("REMOVING UNNECESSARY FILES")
    genie_files = os.listdir(database_to_staging.GENIE_RELEASE_DIR)
    for genie_file in genie_files:
        if (
            genie_version not in genie_file
            and "meta" not in genie_file
            and "case_lists" not in genie_file
        ):
            os.remove(os.path.join(database_to_staging.GENIE_RELEASE_DIR, genie_file))
    os.remove(clinical_path)

    logger.info("REVISE METADATA FILES")
    database_to_staging.revise_metadata_files(
        syn, staging, consortiumSynId, genie_version
    )

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
    logger.info("REMOVING OLD FILES")

    process_functions.rmFiles(database_to_staging.CASE_LIST_PATH)
    private_cna_meta_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR, "genie_private_meta_cna_hg19_seg.txt"
    )
    if os.path.exists(private_cna_meta_path):
        os.unlink(private_cna_meta_path)

    logger.info("CREATING LINK VERSION")
    # Returns release and case list folder
    folders = database_to_staging.create_link_version(
        syn, genie_version, caseListEntities, genePanelEntities, databaseSynIdMappingDf
    )

    if not staging:
        database_to_staging.update_process_trackingdf(
            syn, processTrackerSynId, "SAGE", "dbToStage", start=False
        )

    if not test:
        logger.info("DASHBOARD UPDATE")
        dashboard_table_updater.run_dashboard(
            syn, databaseSynIdMappingDf, genie_version, staging=staging
        )
        generate_dashboard_html(
            genie_version, staging=staging, genie_user=genie_user, genie_pass=genie_pass
        )
        logger.info("DASHBOARD UPDATE COMPLETE")
        logger.info("AUTO GENERATE DATA GUIDE")

    oncotree_version = oncotree_link.split("=")[1]
    data_guide_pdf = generate_data_guide(
        genie_version,
        oncotree_version=oncotree_version,
        database_mapping=databaseSynIdMappingId,
        genie_user=genie_user,
        genie_pass=genie_pass,
    )
    database_to_staging.store_file(
        syn,
        data_guide_pdf,
        genieVersion=genie_version,
        parent=folders["release_folder"],
    )
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
        metavar=184,
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
        pemfile=args.pemFile,
        test=args.test,
        staging=args.staging,
        debug=args.debug,
        skip_mutationsincis=args.skipMutationsInCis,
    )
