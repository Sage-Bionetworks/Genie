![genie banner](https://raw.githubusercontent.com/Sage-Bionetworks/Genie/master/genie_banner.png)

# AACR Project GENIE

[![PyPi](https://img.shields.io/pypi/v/aacrgenie.svg?style=for-the-badge&label=PyPi&logo=PyPi)](https://pypi.org/project/aacrgenie)
[![Docker Automated](https://img.shields.io/docker/automated/sagebionetworks/genie.svg?style=for-the-badge&logo=docker)](https://hub.docker.com/r/sagebionetworks/genie)
[![GitHub CI](https://img.shields.io/github/actions/workflow/status/Sage-Bionetworks/Genie/ci.yml?branch=develop&style=for-the-badge&logo=github)](https://github.com/Sage-Bionetworks/Genie)


## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by institutes participating in AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange). 

## Dependencies

This package contains both R, Python and cli tools.  These are tools or packages you will need, to be able to reproduce these results:
- Python >3.7 or <3.10
    - `pip install -r requirements.txt`
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- R 4.2.2
    - `renv::install()`
    - Follow instructions [here](https://r-docs.synapse.org/#note-for-windows-and-mac-users) to install synapser
- [Java > 8](https://www.java.com/en/download/)
    - For mac users, it seems to work better to run `brew install java`
- [wget](https://www.gnu.org/software/wget/)
    - For mac users, have to run `brew install wget`

## File Validator

One of the features of the `aacrgenie` package is that is provides a local validation tool that GENIE data contributors and install and use to validate their files locally prior to uploading to Synapse.

```
pip install aacrgenie
genie -v
```

This will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.  Please view the help to see how to run to validator.

```
genie validate -h
genie validate data_clinical_supp_SAGE.txt SAGE
```



## Contributing

Please view [contributing guide](CONTRIBUTING.md) to learn how to contribute to the GENIE package.


# Sage Bionetworks Only

## Developing locally

These are instructions on how you would develop and test the pipeline locally.

1. Make sure you have read through the [GENIE Onboarding Docs](https://sagebionetworks.jira.com/wiki/spaces/APGD/pages/2163344270/Onboarding) and have access to all of the required repositories, resources and synapse projects for Main GENIE.
1. Be sure you are invited to the Synapse GENIE Admin team.
1. Make sure you are a Synapse certified user: [Certified User - Synapse User Account Types](https://help.synapse.org/docs/Synapse-User-Account-Types.2007072795.html#SynapseUserAccountTypes-CertifiedUser)
1. Clone this repo and install the package locally.

    ```
    pip install -e .
    pip install -r requirements.txt
    pip install -r requirements-dev.txt
    ```

1. Run the different pipelines on the test project.  The `--project_id syn7208886` points to the test project.

    1. Validate all the files.

        ```
        python bin/input_to_database.py main --project_id syn7208886 --onlyValidate
        ```

    1. Process all the files aside from the mutation (maf, vcf) files.  The mutation processing was split because it takes at least 2 days to process all the production mutation data.  Ideally, there is a parameter to exclude or include file types to process/validate, but that is not implemented.

        ```
        python bin/input_to_database.py main --project_id syn7208886 --deleteOld
        ```

    1. Process the mutation data.  Be sure to clone this repo: https://github.com/Sage-Bionetworks/annotation-tools.  This repo houses the code that re-annotates the mutation data with genome nexus.  The `--createNewMafDatabase` will create a new mutation tables in the test project.  This flag is necessary for production data for two main reasons:
        * During processing of mutation data, the data is appended to the data, so without creating an empty table, there will be duplicated data uploaded.
        * By design, Synapse Tables were meant to be appended to.  When a Synapse Tables is updated, it takes time to index the table and return results. This can cause problems for the pipeline when trying to query the mutation table.  It is actually faster to create an entire new table than updating or deleting all rows and appending new rows when dealing with millions of rows.

        ```
        python bin/input_to_database.py mutation --project_id syn7208886 --deleteOld --genie_annotation_pkg ../annotation-tools --createNewMafDatabase
        ```

    1. Create a consortium release.  Be sure to add the `--test` parameter. Be sure to clone the cbioportal repo: https://github.com/cBioPortal/cbioportal

        ```
        python bin/database_to_staging.py Jan-2017 ../cbioportal TEST --test
        ```

    1. Create a public release.  Be sure to add the `--test` parameter.  Be sure to clone the cbioportal repo: https://github.com/cBioPortal/cbioportal

        ```
        python bin/consortium_to_public.py Jan-2017 ../cbioportal TEST --test
        ```

## Production

The production pipeline is run on Nextflow Tower and the Nextflow workflow is captured in [nf-genie](https://github.com/Sage-Bionetworks-Workflows/nf-genie).  It is wise to create an ec2 via the Sage Bionetworks service catalog to work with the production data, because there is limited PHI in GENIE.
