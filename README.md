![genie banner](https://raw.githubusercontent.com/Sage-Bionetworks/Genie/master/genie_banner.png)

# AACR Project GENIE

[![Docker Automated](https://img.shields.io/docker/automated/sagebionetworks/genie.svg)](https://hub.docker.com/r/sagebionetworks/genie) ![Docker Build](https://img.shields.io/docker/build/sagebionetworks/genie.svg)


## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by institutes participating in AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange). 

## Dependencies

These are tools or packages you will need, to be able to reproduce these results:
- Python 3.6 or higher
- Synapse [command-line client](http://python-docs.synapse.org/CommandLineClient.html) (`pip install synapseclient`)
- Python [pandas](http://pandas.pydata.org/) (`pip install pandas`)
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

## File Validator
```
pip install aacrgenie
genie -v
```

This will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.  Please view the help to see how to run to validator.  
```
genie validate -h
genie validate data_clinical_supp_SAGE.txt SAGE
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
5. Run `genie-job-release` (Make sure to update release version and number)

## Processing on EC2

1. Input to database: `python input_to_database.py -h`
2. Create GENIE Files
**Example Releases**
a. release 4.1-consortium and 4.0-public
```
python database_to_staging.py Jan-2018 ~/cbioportal/ 4.1-consortium --skipMutationsInCis
python consortium_to_public.py Jul-2018 ~/cbioportal/ 4.0-public
```
b. release 5.1-consortium and 5.0-public
```
python database_to_staging.py Jul-2018 ~/cbioportal/ 5.1-consortium
python consortium_to_public.py Jan-2019 ~/cbioportal/ 5.0-public
```


## Instructions to setup batch
1. Build an AMI that can run batch jobs! Start from [this page](https://console.aws.amazon.com/batch/home?region=us-east-1#/first-run) and follow instructions and specify your docker image. It is important at this stage that you time the building of your AMI, or your AMI will not be able to start batch jobs.  After doing so, you will have to start an instance with the AMI and run these 2 commands:

```
sudo stop ecs
sudo rm -rf /var/lib/ecs/data/ecs_agent_data.json
```

2. Rebuild the AMI above, specify the size of the image and put whatever you want in the instance that you would want to bind 

## Adding GENIE sites

1. Invite users to GENIE participant Team 
2. Creates CENTER (input/staging) folder (Set up ACLs) 
3. Update Center Mapping table https://www.synapse.org/#!Synapse:syn10061452/tables/
4. Add center to distribution tables: https://www.synapse.org/#!Synapse:syn10627220/tables/, https://www.synapse.org/#!Synapse:syn7268822/tables/
3. Add users to their GENIE folder
