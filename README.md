![genie banner](https://raw.githubusercontent.com/Sage-Bionetworks/Genie/master/genie_banner.png)

# AACR Project GENIE

[![PyPi](https://img.shields.io/pypi/v/aacrgenie.svg?style=for-the-badge&label=PyPi&logo=PyPi)](https://pypi.org/project/aacrgenie)
[![Docker Automated](https://img.shields.io/docker/automated/sagebionetworks/genie.svg?style=for-the-badge&logo=docker)](https://hub.docker.com/r/sagebionetworks/genie)
[![GitHub CI](https://img.shields.io/github/actions/workflow/status/Sage-Bionetworks/Genie/ci.yml?branch=develop&style=for-the-badge&logo=github)](https://github.com/Sage-Bionetworks/Genie)


## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by institutes participating in AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange). 

## Dependencies

This package contains both R, Python and cli tools.  These are tools or packages you will need, to be able to reproduce these results:
- Python 3.7 or higher
    - `pip install -r requirements.txt`
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- R 4.0.2 or higher
    - `renv::install()`
    - If you have trouble installing synapser, follow instructions [here](https://r-docs.synapse.org/#note-for-windows-and-mac-users)

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

1. Be sure you are invited to the Synapse GENIE Admin team.
1. Install the package locally

    ```
    pip install -e .
    ```

1. Run the pipelines on the test project

    ```
    python bin/input_to_database.py main --project_id syn7208886 --deleteOld
    python bin/input_to_database.py mutation --project_id syn7208886 --deleteOld
    python bin/database_to_staging.py Jan-2017 ./cbioportal TEST --test
    python bin/consortium_to_public.py Jan-2017 ./cbioportal TEST --test
    ```
