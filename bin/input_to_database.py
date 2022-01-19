#! /usr/bin/env python3
"""Script to crawl Synapse folder for a center, validate, and update database tables.

"""
import argparse
from datetime import date
import logging
import os

from genie import input_to_database, write_invalid_reasons, process_functions, config

logger = logging.getLogger(__name__)

# TODO: Remove oncotree_link
# TODO: Remove gneie_annotation_pkg
def main(
    process,
    project_id,
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

    syn = process_functions.synLogin(pemfile, debug=debug)

    # Get the Synapse Project where data is stored
    # Should have annotations to find the table lookup
    project = syn.get(project_id)
    database_to_synid_mapping_synid = project.annotations.get("dbMapping", "")

    databaseToSynIdMapping = syn.tableQuery(
        "SELECT * FROM {}".format(database_to_synid_mapping_synid[0])
    )
    databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()

    center_mapping_id = process_functions.getDatabaseSynId(
        syn, "centerMapping", databaseToSynIdMappingDf=databaseToSynIdMappingDf
    )

    center_mapping = syn.tableQuery("SELECT * FROM %s" % center_mapping_id)
    center_mapping_df = center_mapping.asDataFrame()

    if center is not None:
        assert (
            center in center_mapping_df.center.tolist()
        ), "Must specify one of these centers: {}".format(
            ", ".join(center_mapping_df.center)
        )
        centers = [center]
    else:
        # exclude_sites = ['JHU', 'DFCI', 'GRCC', 'VICC', 'NKI', 'MSK',
        #                  'UHN', 'MDA', 'WAKE', 'YALE', 'UCSF', 'CRUK',
        #                  'CHOP', 'VHIO', 'SCI', 'PHS', 'COLU', 'UCHI']
        center_mapping_df = center_mapping_df[~center_mapping_df["inputSynId"].isnull()]
        # release is a bool column
        center_mapping_df = center_mapping_df[center_mapping_df["release"]]
        # center_mapping_df = center_mapping_df[
        #     ~center_mapping_df['center'].isin(exclude_sites)
        # ]
        centers = center_mapping_df.center

    if oncotree_link is None:
        onco_link = databaseToSynIdMappingDf["Id"][
            databaseToSynIdMappingDf["Database"] == "oncotreeLink"
        ].values[0]
        onco_link_ent = syn.get(onco_link)
        oncotree_link = onco_link_ent.externalURL
    # Check if you can connect to oncotree link,
    # if not then don't run validation / processing
    process_functions.checkUrl(oncotree_link)

    center_mapping_ent = syn.get(center_mapping_id)
    if center_mapping_ent.get("isProcessing", ["True"])[0] == "True":
        raise Exception(
            "Processing/validation is currently happening.  "
            "Please change/add the 'isProcessing' annotation on {} "
            "to False to enable processing".format(center_mapping_id)
        )
    else:
        center_mapping_ent.isProcessing = "True"
        center_mapping_ent = syn.store(center_mapping_ent)
    # remove this query timeout and see what happens
    # syn.table_query_timeout = 50000

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

    format_registry = config.collect_format_types(args.format_registry_packages)
    for process_center in centers:
        input_to_database.center_input_to_database(
            syn,
            project_id,
            process_center,
            process,
            only_validate,
            databaseToSynIdMappingDf,
            center_mapping_df,
            delete_old=delete_old,
            oncotree_link=oncotree_link,
            format_registry=format_registry,
            genie_annotation_pkg=genie_annotation_pkg,
        )

    # To ensure that this is the new entity
    center_mapping_ent = syn.get(center_mapping_id)
    center_mapping_ent.isProcessing = "False"
    center_mapping_ent = syn.store(center_mapping_ent)

    error_tracker_synid = process_functions.getDatabaseSynId(
        syn, "errorTracker", databaseToSynIdMappingDf=databaseToSynIdMappingDf
    )
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
