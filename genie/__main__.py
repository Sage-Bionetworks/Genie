#!/usr/bin/env python3
# noqa pylint: disable=line-too-long
"""genie cli"""
import argparse
import datetime
import logging
import os
import subprocess
import time

import synapseclient

from . import (
    create_case_lists, consortium_to_public, dashboard_table_updater,
    database_to_staging, process_functions,
    validate, write_invalid_reasons
)
from .__version__ import __version__


logger = logging.getLogger(__name__)


def synapse_login(username=None, password=None):
    """
    This function logs into synapse for you if credentials are saved.
    If not saved, then user is prompted username and password.

    :returns:     Synapseclient object
    """
    try:
        syn = synapseclient.login(silent=True)
    except Exception:
        if username is None and password is None:
            raise ValueError("Please specify --syn_user, --syn_pass to specify your Synapse "
                             "login. Please view https://docs.synapse.org/articles/client_configuration.html"
                             "to learn about logging into Synapse via the Python client.")
        syn = synapseclient.login(email=username,
                                  password=password,
                                  silent=True)
    return syn


def perform_create_case_list(_, args):
    """CLI to create case lists given a clinical and assay information file"""
    create_case_lists.main(args.clinical_file_name,
                           args.assay_info_file_name,
                           args.output_dir,
                           args.study_id)


def perform_get_file_errors(syn, args):
    """CLI to get invalid reasons"""
    project = syn.get(args.project_id)
    db_mapping = syn.tableQuery(f"select * from {project.dbMapping[0]}")
    db_mappingdf = db_mapping.asDataFrame()
    error_tracker_synid = db_mappingdf['Id'][
        db_mappingdf['Database'] == "errorTracker"
    ][0]
    center_errors = write_invalid_reasons.get_center_invalid_errors(
        syn, error_tracker_synid
    )
    print(center_errors[args.center])


def perform_public_release(syn, args):
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
    genie_user = os.environ.get('GENIE_USER')
    if args.pemFile is not None:
        genie_pass = process_functions.get_password(args.pemFile)
    else:
        genie_pass = None

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
    command = [cbioValidatorPath, '-s',
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
    folders = database_to_staging.create_link_version(
        syn, args.genieVersion, caseListEntities,
        genePanelEntities, databaseSynIdMappingDf,
        release_type="public"
    )
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
        dashboard_table_updater.run_dashboard(syn, databaseSynIdMappingDf,
                                              args.genieVersion,
                                              staging=args.staging)
        database_to_staging.generate_dashboard_html(
            args.genieVersion, staging=args.staging,
            genie_user=genie_user,
            genie_pass=genie_pass
        )
        logger.info("DASHBOARD UPDATE COMPLETE")
        logger.info("AUTO GENERATE DATA GUIDE")

    onco_link = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'oncotreeLink'
    ].values[0]
    onco_link_ent = syn.get(onco_link)
    oncotree_link = onco_link_ent.externalURL
    oncotree_version = oncotree_link.split("=")[1]

    data_guide_pdf = database_to_staging.generate_data_guide(
        args.genieVersion,
        oncotree_version=oncotree_version,
        database_mapping=databaseSynIdMappingId,
        genie_user=genie_user,
        genie_pass=genie_pass
    )
    data_guide_ent = synapseclient.File(data_guide_pdf,
                                        parent=folders['release_folder'])
    syn.store(data_guide_ent)
    logger.info("COMPLETED CONSORTIUM TO PUBLIC")


def build_parser():
    parser = argparse.ArgumentParser(description='GENIE processing')

    parser.add_argument("--syn_user", type=str, help='Synapse username')

    parser.add_argument("--syn_pass", type=str, help='Synapse password')

    parser.add_argument('-v', '--version', action='version',
                        version='genie {}'.format(__version__))

    subparsers = parser.add_subparsers(title='commands',
                                       description='The following commands are available:',
                                       help='For additional help: "genie <COMMAND> -h"')

    parser_validate = subparsers.add_parser('validate', help='Validates GENIE file formats')

    parser_validate.add_argument("filepath", type=str, nargs="+",
                                 help='File(s) that you are validating.'
                                      'If you validation your clinical files and you have both sample and '
                                      'patient files, you must provide both')

    parser_validate.add_argument("center", type=str, help='Contributing Centers')

    parser_validate.add_argument("--format_registry_packages", type=str, nargs="+",
                                 default=["genie_registry"],
                                 help="Python package name(s) to get valid file formats from (default: %(default)s).")

    parser_validate.add_argument("--oncotree_link", type=str, help="Link to oncotree code")

    validate_group = parser_validate.add_mutually_exclusive_group()

    validate_group.add_argument("--filetype", type=str,
                                help='By default, the validator uses the filename to match '
                                     'the file format.  If your filename is incorrectly named, '
                                     'it will be invalid.  If you know the file format you are '
                                     'validating, you can ignore the filename validation and skip '
                                     'to file content validation. '
                                     'Note, the filetypes with SP at '
                                     'the end are for special sponsored projects.')

    validate_group.add_argument("--parentid", type=str, default=None,
                                help='Synapse id of center input folder. '
                                     'If specified, your valid files will be uploaded '
                                     'to this directory.')

    # TODO: remove this default when private genie project is ready
    parser_validate.add_argument("--project_id", type=str,
                                 default="syn3380222",
                                 help='Synapse Project ID where data is stored. (default: %(default)s).')

    parser_validate.add_argument("--nosymbol-check", action='store_true',
                                 help='Do not check hugo symbols of fusion and cna file')

    parser_validate.set_defaults(func=validate._perform_validate)

    parser_create_case = subparsers.add_parser(
        'create-case-lists', help='Creates cBioPortal case list files'
    )
    parser_create_case.add_argument("clinical_file_name",
                                    type=str,
                                    help="Clinical file path")
    parser_create_case.add_argument("assay_info_file_name",
                                    type=str,
                                    help="gene matrix file path")
    parser_create_case.add_argument("output_dir",
                                    type=str,
                                    help="Output directory")
    parser_create_case.add_argument("study_id",
                                    type=str,
                                    help="Output directory")
    parser_create_case.set_defaults(func=perform_create_case_list)

    parser_get_invalid = subparsers.add_parser(
        'get-file-errors', help='Get the file invalid reasons for a specific center'
    )
    parser_get_invalid.add_argument("center", type=str, help='Contributing Centers')
    parser_get_invalid.add_argument(
        "--project_id", type=str,
        default="syn3380222",
        help='Synapse Project ID where data is stored. (default: %(default)s).'
    )
    parser_get_invalid.set_defaults(func=perform_get_file_errors)


    parser_pr = subparsers.add_parser(
        'public-release', help='Generate public release for GENIE'
    )
    parser_pr.add_argument(
        "processingDate",
        type=str,
        metavar="Jan-2017",
        help="The process date of GENIE in Month-Year format "
             "(ie. Apr-2017)")
    parser_pr.add_argument(
        "cbioportalPath",
        type=str,
        metavar="/path/to/cbioportal",
        help="Make sure you clone the cbioportal github: "
             "git clone https://github.com/cBioPortal/cbioportal.git"
    )
    parser_pr.add_argument(
        "genieVersion",
        type=str,
        help="GENIE public release version"
    )
    parser_pr.add_argument(
        "--publicReleaseCutOff",
        type=int,
        default=366,
        help="Public release cut off time in days (Must "
             "account for leap year, 366)"
    )
    parser_pr.add_argument(
        "--staging",
        action='store_true',
        help="Store into staging folder"
    )
    parser_pr.add_argument(
        "--test",
        action='store_true',
        help="Store into staging folder"
    )
    parser_pr.add_argument(
        "--pemFile",
        type=str,
        help="Path to PEM file (genie.pem)"
    )
    parser_pr.add_argument(
        "--debug",
        action='store_true',
        help="Synapse debug feature"
    )
    parser_pr.set_defaults(func=perform_public_release)

    return parser


def main():
    """Invoke"""
    args = build_parser().parse_args()
    syn = synapse_login(args.syn_user, args.syn_pass)
    # func has to match the set_defaults
    args.func(syn, args)


if __name__ == "__main__":
    main()
