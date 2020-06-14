import argparse
import datetime
import logging
import os
import synapseclient
import subprocess
import time

from genie import dashboard_table_updater
from genie import process_functions
from genie import consortium_to_public
from genie import database_to_staging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):
    cbioValidatorPath = os.path.join(
        args.cbioportalPath,
        "core/src/main/scripts/importer/validateData.py")
    assert os.path.exists(cbioValidatorPath), \
        "Please specify correct cbioportalPath"
    assert not (args.test and args.staging), \
        "You can only specify --test or --staging, not both"
    try:
        processingDate = datetime.datetime.strptime(
            args.processingDate, '%b-%Y')
    except ValueError:
        raise ValueError(
            "Process date must be in the format "
            "abbreviated_month-YEAR ie. Oct-2017")

    syn = process_functions.synLogin(args.pemFile, debug=args.debug)
    # Get all the possible public releases
    # Get configuration
    if args.test:
        databaseSynIdMappingId = 'syn11600968'
        args.genieVersion = "TESTpublic"
    elif args.staging:
        databaseSynIdMappingId = 'syn12094210'
    else:
        databaseSynIdMappingId = 'syn10967259'
    databaseSynIdMapping = syn.tableQuery(
        'select * from %s' % databaseSynIdMappingId)
    databaseSynIdMappingDf = databaseSynIdMapping.asDataFrame()
    public_synid = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'public'].values[0]

    releaseSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'release'].values[0]

    officialPublic = \
        consortium_to_public.get_public_to_consortium_synid_mapping(
            syn, releaseSynId, test=args.test)

    assert args.genieVersion in officialPublic.keys(), \
        "genieVersion must be one of these: {}.".format(
            ", ".join(officialPublic.keys()))

    args.releaseId = officialPublic[args.genieVersion]
    if not args.test and not args.staging:
        processTrackerSynId = databaseSynIdMappingDf['Id'][
            databaseSynIdMappingDf['Database'] == 'processTracker'].values[0]
        processTracker = syn.tableQuery(
            "SELECT timeStartProcessing FROM %s where center = 'SAGE' "
            "and processingType = 'public'" % processTrackerSynId)
        processTrackerDf = processTracker.asDataFrame()
        processTrackerDf['timeStartProcessing'][0] = str(int(time.time()*1000))
        syn.store(synapseclient.Table(processTrackerSynId, processTrackerDf))

    caseListEntities, genePanelEntities = \
        consortium_to_public.consortiumToPublic(
            syn, processingDate, args.genieVersion,
            args.releaseId, databaseSynIdMappingDf,
            publicReleaseCutOff=args.publicReleaseCutOff)

    database_to_staging.revise_metadata_files(syn,
                                              args.staging,
                                              public_synid,
                                              args.genieVersion)

    logger.info("CBIO VALIDATION")
    # Must be exit 0 because the validator sometimes fails,
    # but we still want to capture the output
    command = ['python', cbioValidatorPath, '-s',
               database_to_staging.GENIE_RELEASE_DIR, '-n', '; exit 0']
    cbio_output = subprocess.check_output(" ".join(command), shell=True)
    cbio_decoded_output = cbio_output.decode("utf-8")
    logger.info(cbio_decoded_output)
    if not args.test and not args.staging:
        log_folder_synid = databaseSynIdMappingDf['Id'][
            databaseSynIdMappingDf['Database'] == 'logs'].values[0]
        # Use tempfiles
        cbio_log_file = "cbioValidatorLogsPublic_{}.txt".format(
            args.genieVersion)
        with open(cbio_log_file, "w") as cbioLog:
            cbioLog.write(cbio_decoded_output)
        syn.store(synapseclient.File(cbio_log_file, parentId=log_folder_synid))
        os.remove(cbio_log_file)
    logger.info("REMOVING OLD FILES")
    process_functions.rmFiles(database_to_staging.CASE_LIST_PATH)
    seg_meta_file = '{}/genie_public_meta_cna_hg19_seg.txt'.format(
        database_to_staging.GENIE_RELEASE_DIR)
    if os.path.exists(seg_meta_file):
        os.unlink(seg_meta_file)

    logger.info("CREATING LINK VERSION")
    consortium_to_public.createLinkVersion(
        syn, args.genieVersion, caseListEntities,
        genePanelEntities, databaseSynIdMappingDf)
    # Don't update process tracker is testing or staging
    if not args.test and not args.staging:
        processTracker = syn.tableQuery(
            "SELECT timeEndProcessing FROM %s where center = 'SAGE' and "
            "processingType = 'public'" % processTrackerSynId)
        processTrackerDf = processTracker.asDataFrame()
        processTrackerDf['timeEndProcessing'][0] = str(int(time.time()*1000))
        syn.store(synapseclient.Table(processTrackerSynId, processTrackerDf))

    if not args.test:
        logger.info("DASHBOARD UPDATE")
        dashboard_table_updater.run_dashboard(
            syn, databaseSynIdMappingDf,
            args.genieVersion, staging=args.staging, public=True)
        dashboard_markdown_gen_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'dashboard_markdown_generator.R')
        dashboard_markdown_html_commands = ['Rscript',
                                            dashboard_markdown_gen_path,
                                            args.genieVersion]
        if args.staging:
            dashboard_markdown_html_commands.append('--staging')
        subprocess.check_call(dashboard_markdown_html_commands)
        logger.info("DASHBOARD UPDATE COMPLETE")

    logger.info("COMPLETED CONSORTIUM TO PUBLIC")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("processingDate",
                        type=str,
                        metavar="Jan-2017",
                        help="The process date of GENIE in Month-Year format "
                             "(ie. Apr-2017)")

    parser.add_argument("cbioportalPath",
                        type=str,
                        metavar="/path/to/cbioportal",
                        help="Make sure you clone the cbioportal github: "
                             "git clone https://github.com/cBioPortal/cbioportal.git")

    parser.add_argument("genieVersion",
                        type=str,
                        help="GENIE public release version")

    parser.add_argument("--publicReleaseCutOff",
                        type=int,
                        default=366,
                        help="Public release cut off time in days (Must "
                             "account for leap year, 366)")

    parser.add_argument("--staging",
                        action='store_true',
                        help="Store into staging folder")

    parser.add_argument("--test",
                        action='store_true',
                        help="Store into staging folder")

    parser.add_argument("--pemFile",
                        type=str,
                        help="Path to PEM file (genie.pem)")

    parser.add_argument("--debug",
                        action='store_true',
                        help="Synapse debug feature")
    args = parser.parse_args()
    main(args)
