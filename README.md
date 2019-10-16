![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/veoibd/genie)

# VEO-IBD Consortium Data Processing

## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by centers participating in the VEO-IBD consortium.

## Dependencies

These are tools or packages you will need, to be able to reproduce these results:

- Python 3.5 or higher

## Installation

```console
pip install git+https://github.com/veo-ibd/genie.git
```

This will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.  Please view the help to see how to run to validator.

```console
genie validate -h
```

# Development

## Versioning

1. Update the version in [genie/__version__.py](genie/__version__.py) based on semantic versioning. Use the suffix `-dev` for development branch versions.
2. When releasing, remove the `-dev` from the version.
3. Add a tag and release named the same as the version.

# SAGE BIONETWORKS USE ONLY

## Batch Processing instructions

1. Check docker hub builds to see if theres any failures
2. Log into AWS Batch
3. Run `genie-job-mainprocess`
4. Run `genie-job-mafprocess` (Make sure to add `--createdMafDatabase` flag)
5. Run `genie-job-vcfprocess`
6. Run `genie-job-release` (Make sure to update release version and number)

## Processing on EC2

1. Input to database: `input_to_database.py -h`
2. Create Files

## Example Releases

1. release 4.1-consortium and 4.0-public

    ```console
    python database_to_staging.py Jan-2018 ~/cbioportal/ 4.1-consortium --skipMutationsInCis
    python consortium_to_public.py Jul-2018 ~/cbioportal/ 4.0-public
    ```

1. release 5.1-consortium and 5.0-public

    ```console
    python database_to_staging.py Jul-2018 ~/cbioportal/ 5.1-consortium
    python consortium_to_public.py Jan-2019 ~/cbioportal/ 5.0-public
    ```

## Instructions to setup batch

1. Build an AMI that can run batch jobs! Start from [this page](https://console.aws.amazon.com/batch/home?region=us-east-1#/first-run) and follow instructions and specify your docker image.  It is important at this stage that you time the building of your AMI, or your AMI will not be able to start batch jobs.  After doing so, you will have to start an instance with the AMI and run these 2 commands:

    ```console
    sudo stop ecs
    sudo rm -rf /var/lib/ecs/data/ecs_agent_data.json
    ```

2. Rebuild the AMI above, specify the size of the image and put whatever you want in the instance that you would want to bind

## Adding VEO-IBD sites

1. Invite users to VEO-IBD participant Team
1. Creates CENTER (input/staging) folder (Set up ACLs)
1. Update Center Mapping table https://www.synapse.org/#!Synapse:syn10061452/tables/
1. Add center to distribution tables: https://www.synapse.org/#!Synapse:syn10627220/tables/, https://www.synapse.org/#!Synapse:syn7268822/tables/
1. Add users to their VEO-IBD folder
