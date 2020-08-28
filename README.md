# Synapse Genie

[![Docker Automated](https://img.shields.io/docker/automated/sagebionetworks/genie.svg)](https://hub.docker.com/r/sagebionetworks/genie) ![Docker Build](https://img.shields.io/docker/build/sagebionetworks/genie.svg)


## Introduction

This package can deploy a AACR GENIE like project on Synapse and perform validation and processing of flat files.

## Dependencies

These are tools or packages you will need, to be able to reproduce these results:
- Python 3.6 or higher
- Synapse [command-line client](http://python-docs.synapse.org/CommandLineClient.html) (`pip install synapseclient`)
- Python [pandas](http://pandas.pydata.org/) (`pip install pandas`)


## File Validator
```
pip install synapsegenie
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
