#! /usr/bin/env python3
"""Workflow  crawl Synapse folder for a center, validate, and update database tables.

Known issues:

* All data per center is downloaded every time (even if not changed)

    This increases the length of validation because not all data is changed.
    The complexity is around cross file validation.  There needs to be
    a way to determine relationships between files so that if a file depends on another
    only those files are downloaded.  This is harder for VCFs as there may be
    thousands of them

* Retraction from the sample/patient table is inefficient.  The sample and patient tables will
first contain the data, and it will shortly be removed from the synapse tables.
* Newer implementation of cross-file validation is available but clinical data is still using the
older implementation which unnecesssarily complicates the code.

``` mermaid
flowchart TD
    A["For each center"] --> B["Extract center input files"]
    B --> D{"Has the center uploaded any data"}
    D -- No --> E["Log: No files uploaded"]
    E --> F["End"]
    D -- Yes --> G["Validate files"]
    G --> H{"Are there valid files?"}
    H -- No --> I["Log: No valid files"]
    I --> F
    H -- Yes --> J{"only_validate flag set?"}
    J -- Yes --> K["Log: Validation only"]
    K --> F
    J -- No --> L["Process valid files"]
    %% L --> M["Update process tracker with start time"]
    L --> N["Process files based on fileType"]
    %% N --> O["Update process tracker with end time"]
    N --> P["Upload processed data per file type into internal Synapse Tables and Center staging folders"]
    P --> Q["Retract samples and patients based on retraction tables"]
    Q --> F

    subgraph "File Type Processing"
        N --> N1["clinical"]
        N1 --> N1a["Remap NAACCR and other values to cbioportal accepted values"]
        N1a --> N1b["Map SEQ_DATE to SEQ_YEAR"]
        N1b --> N1d["Redact PHI"]
        N1d --> N1e["Exclude all samples with invalid oncotree codes"]

        N --> N2["maf"]
        N2 --> N2a["Preprocessing: light edits to maf"]
        N2a --> N2b["Re-annotate using Genome Nexus"]

        N --> N3["vcf"]
        N3 --> N3a["Preprocessing: Convert VCF to MAF"]
        N3a --> N2b

        N --> N5["cna"]
        N5 --> N5a["Remap gene symbols using processed BED data"]
        N5a --> N5b["Duplicate genes are merged if possible after remapping"]

        N --> N9["BED"]
        N9 --> N9a["Remap genes to hg19 positions"]
        N9a --> N9b["Gene panel files are created and uploaded to staging folder"]

        N --> N6["assay information, mutation in cis, sample/patient retraction, seg, Structural Variants"]
        N6 --> N6a["Extract data and mild transforms"]

    end

    %% %% Subgraph for get_center_input_files function
    %% subgraph "get_center_input_files Function"
    %%     B --> B1["Iterate over all files uploaded by center"]
    %%     B1 --> B6{"Does file name end with '.vcf' and process != 'mutation'?"}
    %%     B6 -- Yes --> B7["Skip this file"]
    %%     B6 -- No --> B8["Download file"]
    %%     B8 --> B10{"Is it a clinical file?"}
    %%     B10 -- Yes --> B11["Append entity to clinicalpair_entities"]
    %%     B10 -- No --> B12["Append [entity] to prepared_center_file_list"]
    %%     B11 --> B13
    %%     B12 --> B13
    %%     B13{"Is clinicalpair_entities not empty after loop?"}
    %%     B13 -- Yes --> B14["Append clinicalpair_entities to prepared_center_file_list"]
    %%     B14 --> B15["Return prepared_center_file_list"]
    %%     B13 -- No --> B15
    %% end
```
"""
import argparse
import logging
from datetime import date

from genie import (
    config,
    extract,
    input_to_database,
    process_functions,
    write_invalid_reasons,
)

logger = logging.getLogger(__name__)


def main(
    process: str,
    project_id: str,
    center=None,
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

    syn = process_functions.synapse_login(debug=debug)

    # Get project GENIE configurations
    genie_config = extract.get_genie_config(syn=syn, project_id=project_id)

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
        # save maf table to testing or production project as the mode
        new_tables = process_functions.create_new_fileformat_table(
            syn, "vcf2maf", table_name, project_id, project_id
        )
        syn.setPermissions(new_tables["newdb_ent"].id, 3326313, [])
        genie_config["vcf2maf"] = new_tables["newdb_ent"].id

    # Get file format classes
    format_registry = config.collect_format_types(args.format_registry_packages)

    # Start GENIE processing
    for process_center in centers:
        # Check if the genie genome nexus is up, if not then don't run
        # processing
        process_functions.checkUrl("http://genie.genomenexus.org/")
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
        delete_old=args.deleteOld,
        only_validate=args.onlyValidate,
        oncotree_link=args.oncotree_link,
        create_new_maf_database=args.createNewMafDatabase,
        debug=args.debug,
        genie_annotation_pkg=args.genie_annotation_pkg,
        format_registry=args.format_registry_packages,
    )
