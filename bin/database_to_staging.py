import argparse
import datetime
import logging
import os
import subprocess
import synapseclient
# import time

from genie import database_to_staging
from genie import process_functions
from genie import create_case_lists
from genie import dashboard_table_updater

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main(genie_version,
         processing_date,
         cbioportal_path,
         oncotree_link=None,
         consortium_release_cutoff=184,
         pemfile=None,
         test=False,
         staging=False,
         debug=False,
         skip_mutationsincis=False):
    '''
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
    '''
    syn = process_functions.synLogin(pemfile, debug=debug)
    genie_user = os.environ.get('GENIE_USER')
    if pemfile is not None:
        genie_pass = process_functions.get_password(pemfile)
    else:
        genie_pass = None

    if test:
        databaseSynIdMappingId = 'syn11600968'
        genie_version = "TESTING"
    elif staging:
        skip_mutationsincis = True
        databaseSynIdMappingId = 'syn12094210'
    else:
        databaseSynIdMappingId = 'syn10967259'
    # Database/folder syn id mapping
    databaseSynIdMapping = syn.tableQuery(
        'select * from {}'.format(databaseSynIdMappingId))
    databaseSynIdMappingDf = databaseSynIdMapping.asDataFrame()
    # databaseSynIdMappingDf.index = databaseSynIdMappingDf.Database
    # del databaseSynIdMappingDf['Database']
    # databaseSynIdMappingDf.to_dict()

    if oncotree_link is None:
        oncoLink = databaseSynIdMappingDf['Id'][
            databaseSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
        oncoLinkEnt = syn.get(oncoLink)
        oncotree_link = oncoLinkEnt.externalURL

    # Check if you can connect to oncotree link,
    # if not then don't run validation / processing
    process_functions.checkUrl(oncotree_link)

    cbioValidatorPath = os.path.join(
        cbioportal_path, "core/src/main/scripts/importer/validateData.py")
    assert os.path.exists(cbioValidatorPath),\
        "Please specify correct cbioportalPath"
    syn.table_query_timeout = 50000

    consortiumSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'consortium'].values[0]
    processTrackerSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'processTracker'].values[0]
    # get syn id of case list folder in consortium release
    # caseListSynId = findCaseListId(syn, consortiumSynId)
    caseListSynId, already_exists = \
        database_to_staging.search_and_create_folder(
            syn, consortiumSynId, "case_lists")

    if not staging:
        database_to_staging.update_process_trackingdf(
            syn, processTrackerSynId, 'SAGE', 'dbToStage', start=True)

    centerMappingSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'centerMapping'].values[0]
    # Only release files where release is true
    center_mapping = syn.tableQuery(
        'SELECT * FROM {} where release is true'.format(centerMappingSynId))
    center_mappingdf = center_mapping.asDataFrame()
    processingDate = datetime.datetime.strptime(processing_date, '%b-%Y')

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
        genie_pass=genie_pass)

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
        'data_clinical_{}.txt'.format(genie_version))
    gene_matrix_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "data_gene_matrix_{}.txt".format(genie_version))
    create_case_lists.main(
        clinical_path,
        gene_matrix_path,
        database_to_staging.CASE_LIST_PATH,
        "genie_private")
    caseListFiles = os.listdir(database_to_staging.CASE_LIST_PATH)
    caseListEntities = []
    for casePath in caseListFiles:
        casePath = os.path.join(database_to_staging.CASE_LIST_PATH, casePath)
        caseListEntities.append(database_to_staging.storeFile(
            syn,
            casePath,
            parent=caseListSynId,
            staging=staging,
            caseLists=True,
            genieVersion=genie_version))

    logger.info("REMOVING UNNECESSARY FILES")
    genie_files = os.listdir(database_to_staging.GENIE_RELEASE_DIR)
    for genieFile in genie_files:
        if genie_version not in genieFile and \
             "meta" not in genieFile and "case_lists" not in genieFile:
            os.remove(os.path.join(database_to_staging.GENIE_RELEASE_DIR,
                                   genieFile))
    os.remove(clinical_path)

    logger.info("REVISE METADATA FILES")
    database_to_staging.revise_metadata_files(syn, staging,
                                              consortiumSynId,
                                              genie_version)

    logger.info("CBIO VALIDATION")
    '''
    Must be exit 0 because the validator sometimes fails,
    but we still want to capture the output
    '''
    command = [cbioValidatorPath, '-s', database_to_staging.GENIE_RELEASE_DIR,
               '-n', '; exit 0']
    cbioOutput = subprocess.check_output(" ".join(command), shell=True)
    logger.info(cbioOutput.decode("utf-8"))

    cbio_validator_log = \
        "cbioValidatorLogsConsortium_{}.txt".format(genie_version)
    if not test and not staging:
        log_folder_synid = databaseSynIdMappingDf['Id'][
            databaseSynIdMappingDf['Database'] == 'logs'].values[0]
        with open(cbio_validator_log, "w") as cbioLog:
            cbioLog.write(cbioOutput.decode("utf-8"))
        syn.store(synapseclient.File(
            cbio_validator_log, parentId=log_folder_synid))
        os.remove(cbio_validator_log)
    logger.info("REMOVING OLD FILES")

    process_functions.rmFiles(database_to_staging.CASE_LIST_PATH)
    private_cna_meta_path = \
        '%s/genie_private_meta_cna_hg19_seg.txt' % database_to_staging.GENIE_RELEASE_DIR
    if os.path.exists(private_cna_meta_path):
        os.unlink(private_cna_meta_path)

    logger.info("CREATING LINK VERSION")
    database_to_staging.create_link_version(
        syn, genie_version, caseListEntities,
        genePanelEntities, databaseSynIdMappingDf)

    if not staging:
        database_to_staging.update_process_trackingdf(
            syn, processTrackerSynId, 'SAGE', 'dbToStage', start=False)

    logger.info("COMPLETED DATABASE TO STAGING")

    if not test:
        logger.info("DASHBOARD UPDATE")
        dashboard_table_updater.run_dashboard(
            syn,
            databaseSynIdMappingDf,
            genie_version,
            staging=staging)
        dashboard_markdown_html_commands = [
            'Rscript',
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'dashboard_markdown_generator.R'),
            genie_version]

        if genie_user is not None and genie_pass is not None:
            dashboard_markdown_html_commands.extend(
                ['--syn_user', genie_user, '--syn_pass', genie_pass])
        if staging:
            dashboard_markdown_html_commands.append('--staging')
        subprocess.check_call(dashboard_markdown_html_commands)
        logger.info("DASHBOARD UPDATE COMPLETE")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Release GENIE consortium files')

    parser.add_argument(
        "processingDate",
        type=str,
        metavar="Jan-2017",
        help="The processing date of GENIE in Month-Year format"
             " (ie. Apr-2017)")

    parser.add_argument(
        "cbioportalPath",
        type=str,
        metavar="/path/to/cbioportal",
        help="Make sure you clone the cbioportal github: "
             "git clone https://github.com/cBioPortal/cbioportal.git")

    parser.add_argument(
        "genieVersion",
        type=str,
        metavar="1.0.1",
        help="GENIE release version")

    parser.add_argument(
        "--oncotreeLink",
        type=str,
        help="Link to oncotree code")

    parser.add_argument(
        "--consortiumReleaseCutOff",
        type=int,
        metavar=184,
        default=184,
        help="Consortium release cut off time in days")

    parser.add_argument(
        "--skipMutationsInCis",
        action='store_true',
        help="Skip running mutation in cis script")

    parser.add_argument(
        "--pemFile",
        type=str,
        help="Path to PEM file (genie.pem)")

    parser.add_argument(
        "--debug",
        action='store_true',
        help="Synapse debug feature")

    test_group = parser.add_mutually_exclusive_group()
    test_group.add_argument(
        "--test",
        action='store_true',
        help="Run test")

    test_group.add_argument(
        "--staging",
        action='store_true',
        help="Store into staging folder")

    args = parser.parse_args()

    main(genie_version=args.genieVersion,
         processing_date=args.processingDate,
         cbioportal_path=args.cbioportalPath,
         oncotree_link=args.oncotreeLink,
         consortium_release_cutoff=args.consortiumReleaseCutOff,
         pemfile=args.pemFile,
         test=args.test,
         staging=args.staging,
         debug=args.debug,
         skip_mutationsincis=args.skipMutationsInCis)