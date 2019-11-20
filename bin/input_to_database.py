#! /usr/bin/env python3
import os
import sys
import json
import argparse
import logging

from genie import config
from genie import input_to_database
from genie import write_invalid_reasons
from genie import process_functions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_config(config_obj):
    """Utility to get configuration from different sources.

    Args:
        config_obj: An object holding a configuration.
        
    Currently allows a path to a JSON file.
    """

    if os.path.exists(config_obj) and config_obj.endswith('json'):
        config = json.load(open(config_obj))
    
    return config

def get_processing_status(syn, center_mapping_id):
    """Determine if processing is in progress.
    
    """
    center_mapping_ent = syn.get(center_mapping_id, downloadFile=False)
    return center_mapping_ent.get('isProcessing', ['True'])[0] == 'True'

def set_processing_status(syn, center_mapping_id, status):
    """Set processing status.
    
    """
    center_mapping_ent = syn.get(center_mapping_id, downloadFile=False)
    center_mapping_ent.isProcessing = str(status)
    center_mapping_ent = syn.store(center_mapping_ent)

    return status

def main(process,
         project_id,
         center=None,
         pemfile=None,
         delete_old=False,
         only_validate=False,
         oncotree_link=None,
         create_new_maf_database=False,
         debug=False,
         reference=None,
         vcf2maf_path=None,
         vep_path=None,
         vep_data=None):

    syn = process_functions.synLogin(pemfile, debug=debug)
    # Must specify correct paths to vcf2maf, VEP and VEP data
    # if trying to process vcf, maf and mafSP
    if process in ['vcf', 'maf', 'mafSP'] and not only_validate:
        assert os.path.exists(vcf2maf_path), (
            "Path to vcf2maf (--vcf2mafPath) must be specified "
            "if `--process {vcf,maf,mafSP}` is used")
        assert os.path.exists(vep_path), (
            "Path to VEP (--vepPath) must be specified "
            "if `--process {vcf,maf,mafSP}` is used")
        assert os.path.exists(vep_data), (
            "Path to VEP data (--vepData) must be specified "
            "if `--process {vcf,maf,mafSP}` is used")

    # Get the Synapse Project where data is stored
    # Should have annotations to find the table lookup
    project = syn.get(project_id)
    database_to_synid_mapping_synid = project.annotations.get("dbMapping", "")[0]
    
    databaseToSynIdMapping = syn.tableQuery(
        'SELECT * FROM {}'.format(database_to_synid_mapping_synid))
    databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()

    center_mapping_id = process_functions.getDatabaseSynId(
        syn, "centerMapping",
        databaseToSynIdMappingDf=databaseToSynIdMappingDf)

    center_mapping = syn.tableQuery('SELECT * FROM %s' % center_mapping_id)
    center_mapping_df = center_mapping.asDataFrame()

    if center is not None:
        assert center in center_mapping_df.center.tolist(), (
            "Must specify one of these centers: {}".format(
                ", ".join(center_mapping_df.center)))
        centers = [center]
    else:
        center_mapping_df = \
            center_mapping_df[~center_mapping_df['inputSynId'].isnull()]
        # release is a bool column
        center_mapping_df = center_mapping_df[center_mapping_df['release']]
        centers = center_mapping_df.center

    if oncotree_link is None:
        try:
            onco_link = databaseToSynIdMappingDf['Id'][
                databaseToSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
            onco_link_ent = syn.get(onco_link)
            oncotree_link = onco_link_ent.externalURL
            # Check if you can connect to oncotree link,
            # if not then don't run validation / processing
            process_functions.checkUrl(oncotree_link)
        except:
            oncotree_link = None

    center_mapping_ent = syn.get(center_mapping_id)
    # remove this query timeout and see what happens
    syn.table_query_timeout = 50000

    # Create new maf database, should only happen once if its specified
    if create_new_maf_database:
        databaseToSynIdMappingDf = \
            input_to_database.create_and_archive_maf_database(syn, databaseToSynIdMappingDf)

    for center in centers:
        input_to_database.center_input_to_database(
            syn, project_id, center, process,
            only_validate,
            vcf2maf_path, vep_path,
            vep_data, databaseToSynIdMappingDf,
            center_mapping_df, reference=reference,
            delete_old=delete_old,
            oncotree_link=oncotree_link)

    # To ensure that this is the new entity

    error_tracker_synid = process_functions.getDatabaseSynId(
        syn, "errorTracker", databaseToSynIdMappingDf=databaseToSynIdMappingDf)
    # Only write out invalid reasons if the center
    # isnt specified and if only validate
    if center is None and only_validate:
        logger.info("WRITING INVALID REASONS TO CENTER STAGING DIRS")
        write_invalid_reasons.write_invalid_reasons(
            syn, center_mapping_df, error_tracker_synid)

if __name__ == "__main__":
    '''
    Argument parsers
    TODO: Fix case of arguments
    '''
    parser = argparse.ArgumentParser(
        description='GENIE center ')

    parser.add_argument('--config', help='JSON config file.')

    parser.add_argument(
        "process",
        choices=['vcf', 'maf', 'main', 'mafSP'],
        help='Process vcf, maf or the rest of the files')
    parser.add_argument(
        '--center',
        help='The centers')
    parser.add_argument(
        "--pemFile",
        type=str,
        help="Path to PEM file (genie.pem)")
    parser.add_argument(
        "--deleteOld",
        action='store_true',
        help="Delete all old processed and temp files")
    parser.add_argument(
        "--onlyValidate",
        action='store_true',
        help="Only validate the files, don't process")
    parser.add_argument(
        "--oncotree_link",
        type=str,
        help="Link to oncotree code")
    parser.add_argument(
        "--createNewMafDatabase",
        action='store_true',
        help="Creates a new maf database")
    parser.add_argument(
        "--project_id",
        type=str,
        help="Synapse Project ID where data is stored.")
    parser.add_argument(
        "--debug",
        action='store_true',
        help="Add debug mode to synapse")
    parser.add_argument(
        "--reference",
        type=str,
        help="Path to VCF reference file")
    parser.add_argument("--format_registry_packages", type=str, nargs="+",
                         default="genie",
                         help="Python package name(s) to get valid file formats from (default: %(default)s).")

    # DEFAULT PARAMS
    parser.add_argument(
        "--vcf2mafPath",
        type=str,
        help="Path to vcf2maf",
        default=os.path.expanduser("~/vcf2maf-1.6.14"))
    parser.add_argument(
        "--vepPath",
        type=str,
        help="Path to VEP",
        default=os.path.expanduser("~/vep"))
    parser.add_argument(
        "--vepData",
        type=str,
        help="Path to VEP data",
        default=os.path.expanduser("~/.vep"))

    args = parser.parse_args()

    main(args.process,
         project_id=args.project_id,
         center=args.center,
         pemfile=args.pemFile,
         delete_old=args.deleteOld,
         only_validate=args.onlyValidate,
         oncotree_link=args.oncotree_link,
         create_new_maf_database=args.createNewMafDatabase,
         debug=args.debug,
         reference=args.reference,
         vcf2maf_path=args.vcf2mafPath,
         vep_path=args.vepPath,
         vep_data=args.vepData)
