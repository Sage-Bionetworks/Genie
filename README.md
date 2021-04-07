![genie banner](https://raw.githubusercontent.com/Sage-Bionetworks/Genie/master/genie_banner.png)

# AACR Project GENIE

[![PyPi](https://img.shields.io/pypi/v/aacrgenie.svg?style=for-the-badge&label=PyPi&logo=PyPi)](https://pypi.org/project/aacrgenie)
[![Docker Automated](https://img.shields.io/docker/automated/sagebionetworks/genie.svg?style=for-the-badge&logo=docker)](https://hub.docker.com/r/sagebionetworks/genie)
[![GitHub CI](https://img.shields.io/github/workflow/status/Sage-Bionetworks/Genie/build.svg?&style=for-the-badge&logo=github)](https://github.com/nlpsandbox/nlpsandbox-client)


## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by institutes participating in AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange). 

## Dependencies

These are tools or packages you will need, to be able to reproduce these results:
- Python 3.6 or higher
- `pip install -r requirements.txt`
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

The GENIE project follows the standard [git flow](https://guides.github.com/introduction/flow/) development strategy.

## Versioning

1. Update the version in [genie/__version__.py](genie/__version__.py) based on semantic versioning. Use the suffix `-dev` for development branch versions.
1. When releasing, remove the `-dev` from the version.
1. Add a tag and release named the same as the version.
