#! /usr/bin/env python3
"""Script to crawl Synapse folder for a center, validate, and update database tables.

"""
import argparse
from datetime import date
import logging
import os

from genie import (
    config,
    input_to_database,
    process_functions,
    write_invalid_reasons,
    validate
)

logger = logging.getLogger(__name__)


# TODO: Remove oncotree_link
def main(
    process: str,
    project_id: str,
    center=None,
    pemfile=None,
    delete_old=False,
    only_validate=False,
    oncotree_link=None,
    genie_annotation_pkg=None,
    create_new_maf_database=False,
    debug=False,
    format_registry=None,
):
    """_summary_

    Args:
        process (str): main or mutation processing
        project_id (str): Synapse project id that houses GENIE project
        center (str, optional): GENIE center. Defaults to None.
        pemfile (str, optional): Path to private key. Defaults to None.
        delete_old (bool, optional): True to delete all old input/processed files.
                                     Defaults to False.
        only_validate (bool, optional): True if only validate files. Defaults to False.
        oncotree_link (str, optional): Link to oncotree version. Defaults to None.
        genie_annotation_pkg (str, optional): vcf/maf conversion tools.
                                              Defaults to None.
        create_new_maf_database (bool, optional): To create new maf table.
                                                  Defaults to False.
        debug (bool, optional): Debug mode. Defaults to False.
        format_registry (str, optional): File format registry python package.
                                         Defaults to None.

    Raises:
        ValueError: _description_
        Exception: _description_
    """

    syn = process_functions.synLogin(pemfile, debug=debug)

    # Get the Synapse Project where data is stored
    # Should have annotations to find the table lookup
    project = syn.get(project_id)

    # Get project GENIE configurations
    database_to_synid_mapping_synid = project.annotations.get("dbMapping", "")
    genie_config = validate.get_config(
        syn=syn, synid=database_to_synid_mapping_synid[0]
    )

    # TODO: remove this once code is using genie config
    databaseToSynIdMappingDf = process_functions.get_syntabledf(
        syn=syn, query_string=f"SELECT * FROM {database_to_synid_mapping_synid[0]}"
    )

    center_mapping_id = genie_config['centerMapping']
    center_mapping_df = process_functions.get_syntabledf(
        syn=syn,
        query_string=f"SELECT * FROM {center_mapping_id} where release is true"
    )

    # Filter for specific center
    if center is not None:
        if center not in center_mapping_df.center.tolist():
            raise ValueError("Must specify one of these centers: {}".format(
                ", ".join(center_mapping_df.center)
            ))
        center_mapping_df = center_mapping_df[center_mapping_df['center'] == center]
        centers = [center]
    else:
        # exclude_sites = ['JHU', 'DFCI', 'GRCC', 'VICC', 'NKI', 'MSK',
        #                  'UHN', 'MDA', 'WAKE', 'YALE', 'UCSF', 'CRUK',
        #                  'CHOP', 'VHIO', 'SCI', 'PHS', 'COLU', 'UCHI']
        # center_mapping_df = center_mapping_df[~center_mapping_df["inputSynId"].isnull()]
        # release is a bool column
        # center_mapping_df = center_mapping_df[center_mapping_df["release"]]
        # center_mapping_df = center_mapping_df[
        #     ~center_mapping_df['center'].isin(exclude_sites)
        # ]
        centers = center_mapping_df.center
    center_mapping_df.index = center_mapping_df.center
    # TODO: remove this after code uses genie_config
    # del center_mapping_df['center']
    # Add center configurations including input/staging synapse ids
    genie_config['center_config'] = center_mapping_df.to_dict('index')

    # Get actual oncotree link
    if oncotree_link is None:
        onco_link_ent = syn.get(genie_config['oncotreeLink'])
        oncotree_link = onco_link_ent.externalURL
        genie_config['oncotreeLink'] = onco_link_ent.externalURL
    else:
        genie_config['oncotreeLink'] = oncotree_link
    # Check if you can connect to oncotree link,
    # if not then don't run validation / processing
    process_functions.checkUrl(genie_config['oncotreeLink'])

    # Add genie annotation package to config
    genie_config['genie_annotation_pkg'] = genie_annotation_pkg

    center_mapping_ent = syn.get(center_mapping_id)
    if center_mapping_ent.get("isProcessing", ["True"])[0] == "True":
        raise Exception(
            "Processing/validation is currently happening.  "
            f"Please change/add the 'isProcessing' annotation on {center_mapping_id} "
            "to False to enable processing"
        )
    else:
        center_mapping_ent.isProcessing = "True"
        center_mapping_ent = syn.store(center_mapping_ent)

    # Create new maf database, should only happen once if its specified
    if create_new_maf_database:
        today = date.today()
        table_name = f"Narrow MAF Database - {today}"
        # filetype = "vcf2maf"
        # syn7208886 is the GENIE staging project to archive maf table
        new_tables = process_functions.create_new_fileformat_table(
            syn, "vcf2maf", table_name, project_id, "syn7208886"
        )
        syn.setPermissions(new_tables["newdb_ent"].id, 3326313, [])
        databaseToSynIdMappingDf = new_tables["newdb_mappingdf"]
        genie_config['vcf2maf'] = new_tables["newdb_ent"].id

    # Get file format classes
    format_registry = config.collect_format_types(args.format_registry_packages)

    # Start GENIE processing
    for process_center in centers:
        input_to_database.center_input_to_database(
            syn=syn,
            project_id=project_id,
            center=process_center,
            process=process,
            only_validate=only_validate,
            database_to_synid_mappingdf=databaseToSynIdMappingDf,
            delete_old=delete_old,
            format_registry=format_registry,
            genie_config=genie_config
        )

    # To ensure that this is the new entity
    center_mapping_ent = syn.get(center_mapping_id)
    center_mapping_ent.isProcessing = "False"
    center_mapping_ent = syn.store(center_mapping_ent)

    error_tracker_synid = genie_config['errorTracker']
    # Only write out invalid reasons if the center
    # isnt specified and if only validate
    if center is None and only_validate:
        logger.info("WRITING INVALID REASONS TO CENTER STAGING DIRS")
        write_invalid_reasons.write(syn, center_mapping_df, error_tracker_synid)
    logger.info("INPUT TO DATABASE COMPLETE")


if __name__ == "__main__":
    # Argument parsers
    # TODO: Fix case of arguments
    parser = argparse.ArgumentParser(description="GENIE center ")
    parser.add_argument(
        "process",
        choices=["mutation", "main"],
        help="Process vcf, maf or the rest of the files",
    )
    parser.add_argument(
        "--project_id", help="Synapse Project ID where data is stored.", required=True
    )
    parser.add_argument("--center", help="The centers")
    parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")
    parser.add_argument(
        "--deleteOld",
        action="store_true",
        help="Delete all old processed and temp files",
    )
    parser.add_argument(
        "--onlyValidate",
        action="store_true",
        help="Only validate the files, don't process",
    )
    parser.add_argument("--oncotree_link", type=str, help="Link to oncotree code")
    parser.add_argument(
        "--createNewMafDatabase", action="store_true", help="Creates a new maf database"
    )
    parser.add_argument(
        "--debug", action="store_true", help="Add debug mode to synapse"
    )
    parser.add_argument("--genie_annotation_pkg", help="GENIE annotation pkg")

    # DEFAULT PARAMS
    parser.add_argument(
        "--format_registry_packages",
        type=str,
        nargs="+",
        default=["genie_registry"],
        help="Python package name(s) to get valid file formats from "
        "(default: %(default)s).",
    )
    args = parser.parse_args()

    main(
        args.process,
        project_id=args.project_id,
        center=args.center,
        pemfile=args.pemFile,
        delete_old=args.deleteOld,
        only_validate=args.onlyValidate,
        oncotree_link=args.oncotree_link,
        create_new_maf_database=args.createNewMafDatabase,
        debug=args.debug,
        genie_annotation_pkg=args.genie_annotation_pkg,
        format_registry=args.format_registry_packages,
    )
