#! /usr/bin/env python3
"""Script to crawl Synapse folder for a center, validate, and update database tables.

"""
import argparse
from datetime import date
import logging

from genie import (
    config,
    input_to_database,
    process_functions,
    write_invalid_reasons,
)

logger = logging.getLogger(__name__)


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
    """Invoke the GENIE ETL pipeline from data input files to synapse tables

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
        ValueError: If invalid center name is specified
        Exception: If processing is already happening.
    """

    syn = process_functions.synLogin(pemfile, debug=debug)

    # Get project GENIE configurations
    genie_config = process_functions.get_genie_config(syn=syn, project_id=project_id)

    # Filter for specific center
    if center is not None:
        if center not in genie_config["center_config"].keys():
            raise ValueError(
                "Must specify one of these centers: {}".format(
                    ", ".join(genie_config["center_config"].keys())
                )
            )
        centers = [center]
    else:
        # TODO: add in logic to exclude sites from processing
        centers = list(genie_config["center_config"].keys())

    # HACK: Modify oncotree link config
    if oncotree_link is None:
        onco_link_ent = syn.get(genie_config["oncotreeLink"])
        oncotree_link = onco_link_ent.externalURL
    genie_config["oncotreeLink"] = oncotree_link
    # Check if you can connect to oncotree link,
    # if not then don't run validation / processing
    process_functions.checkUrl(genie_config["oncotreeLink"])

    # HACK: Add genie annotation package to config
    if process == "mutation" and genie_annotation_pkg is None:
        raise ValueError("Must define genie annotation pkg if mutation processing")
    genie_config["genie_annotation_pkg"] = genie_annotation_pkg

    # HACK: This is essential, because Synapse has concurrency update issues
    center_mapping_ent = syn.get(genie_config["centerMapping"])
    if center_mapping_ent.get("isProcessing", ["True"])[0] == "True":
        raise Exception(
            "Processing/validation is currently happening.  Please change/add the "
            f"'isProcessing' annotation on {genie_config['centerMapping']} "
            "to False to enable processing"
        )
    else:
        center_mapping_ent.isProcessing = "True"
        center_mapping_ent = syn.store(center_mapping_ent)

    # HACK: Create new maf database, should only happen once if its specified
    # Will modify genie configuration
    if create_new_maf_database:
        today = date.today()
        table_name = f"Narrow MAF Database - {today}"
        # filetype = "vcf2maf"
        # syn7208886 is the GENIE staging project to archive maf table
        new_tables = process_functions.create_new_fileformat_table(
            syn, "vcf2maf", table_name, project_id, "syn7208886"
        )
        syn.setPermissions(new_tables["newdb_ent"].id, 3326313, [])
        genie_config["vcf2maf"] = new_tables["newdb_ent"].id

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
            delete_old=delete_old,
            format_registry=format_registry,
            genie_config=genie_config,
        )

    # HACK: To ensure that this is the new entity
    center_mapping_ent = syn.get(genie_config["centerMapping"])
    center_mapping_ent.isProcessing = "False"
    # No need to return ent variable because it is unused
    syn.store(center_mapping_ent)

    error_tracker_synid = genie_config["errorTracker"]
    # Only write out invalid reasons if the center
    # isnt specified and if only validate
    if center is None and only_validate:
        logger.info("WRITING INVALID REASONS TO CENTER STAGING DIRS")
        write_invalid_reasons.write(
            syn, genie_config["centerMapping"], error_tracker_synid
        )
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
