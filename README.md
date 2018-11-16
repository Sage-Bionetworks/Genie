# AACR Project GENIE

## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by institutes participating in AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange). 

## Dependencies

These are tools or packages you will need, to be able to reproduce these results:
- A Linux machine or a cluster of compute nodes with a job scheduler like LSF (`bsub`) or SGE (`qsub`)
- Python 2.7.10 or higher
- Sage Synapse's [command-line client](http://python-docs.synapse.org/CommandLineClient.html) (`pip install synapseclient`)
- Python [pandas](http://pandas.pydata.org/) (`pip install pandas`)
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

## File Validator
```
pip install git+https://github.com/Sage-Bionetworks/Genie.git
```

This will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.  Please view the help to see how to run to validator.  
```
genie validate -h
genie validate clinical data_clincal_supp_SAGE.txt SAGE
```

# SAGE USE ONLY
## Processing instructions
1. Go onto GENIE ec2

2. Run shell script - `genie/processGENIE.sh` 

**Releases**

a. release 1.1.0 - 2.0.0
```
python database_to_staging.py Jan-2017 ~/cbioportal/ 1.1.0 --skipMutationsInCis
python consortium_to_public.py Jul-2017 ~/cbioportal/ 2.0.0
```
b. release 2.1.0 - 3.0.0
```
python database_to_staging.py Jul-2017 ~/cbioportal/ 2.1.0
python consortium_to_public.py Jan-2018 ~/cbioportal/ 3.0.0
```
c. release 4.1-consortium and 4.0-public
```
python database_to_staging.py Jan-2018 ~/cbioportal/ 4.1-consortium
python consortium_to_public.py Jul-2018 ~/cbioportal/ 4.0-public
```
d. release 5.1-consortium and 5.0-public
```
python database_to_staging.py Jul-2018 ~/cbioportal/ 5.1-consortium
python consortium_to_public.py Jan-2019 ~/cbioportal/ 5.0-public
```

3. Run dashboard scripts and check dashboard page to make sure numbers are correct

4. If numbers don't match up, please check the log files to make sure no errors came up during processing.

## Instructions for batch
1. Build an AMI that can run batch jobs! Start from [this page](https://console.aws.amazon.com/batch/home?region=us-east-1#/first-run) and follow instructions and specify your docker image.  It is important at this stage that you time the building of your AMI, or your AMI will not be able to start batch jobs.  After doing so, you will have to start an instance with the AMI and run these 2 commands:

```
sudo stop ecs
sudo rm -rf /var/lib/ecs/data/ecs_agent_data.json
```

2. Rebuild the AMI above, specify the size of the image and put whatever you want in the instance that you would want to bind 

## Adding GENIE sites

1. Invite users to GENIE participant Team 
2. Run python utility_sage/addCenter.py https://github.com/Sage-Bionetworks/Genie/blob/master/utility_sage/addCenter.py 
* Creates CENTER (input/staging) folder (Set up ACLs) 
* Update Center Mapping table https://www.synapse.org/#!Synapse:syn10061452/tables/
* Add center to distribution tables: https://www.synapse.org/#!Synapse:syn10627220/tables/, https://www.synapse.org/#!Synapse:syn7268822/tables/
3. Add users to their GENIE folder

