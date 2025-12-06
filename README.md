![genie banner](https://raw.githubusercontent.com/Sage-Bionetworks/Genie/master/genie_banner.png)

# AACR Project GENIE

[![PyPi](https://img.shields.io/pypi/v/aacrgenie.svg?style=for-the-badge&label=PyPi&logo=PyPi)](https://pypi.org/project/aacrgenie)
[![GHCR Docker Package](https://img.shields.io/badge/ghcr.io-sage--bionetworks%2Fgenie-blue?style=for-the-badge&logo=github)](https://github.com/orgs/sage-bionetworks/packages/container/package/genie)
[![GitHub CI](https://img.shields.io/github/actions/workflow/status/Sage-Bionetworks/Genie/ci.yml?branch=develop&style=for-the-badge&logo=github)](https://github.com/Sage-Bionetworks/Genie)


## Table of Contents

- [Introduction](#introduction)
- [Documentation](#documentation)
- [Dependencies](#dependencies)
- [File Validator](#file-validator)
  - [Setting up your environment](#setting-up-your-environment)
  - [Running the validator](#running-the-validator)
  - [Example commands](#example-commands)
- [Contributing](#contributing)
- [Sage Bionetworks Only](#sage-bionetworks-only)
  - [Running locally](#running-locally)
    - [Using conda](#using-conda)
    - [Using pipenv](#using-pipenv)
    - [Using docker (**HIGHLY** Recommended)](#using-docker-highly-recommended)
    - [Setting up](#setting-up)
  - [Developing](#developing)
    - [Developing with Docker](#developing-with-docker)
    - [Modifying Docker](#modifying-docker)
- [Testing](#testing)
  - [Running unit tests](#running-unit-tests)
  - [Running integration tests](#running-integration-tests)
- [Production](#production)
- [Github Workflows](#github-workflows)


## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by institutes participating in AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange).

## Documentation

For more information about the AACR genie repository, [visit the GitHub Pages site.](https://sage-bionetworks.github.io/Genie/)

## Dependencies

This package contains both R, Python and cli tools.  These are tools or packages you will need, to be able to reproduce these results:
- Python >=3.10 or <3.12
    - `pip install -r requirements.txt`
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- R 4.3.3
    - `renv::install()`
    - Follow instructions [here](https://r-docs.synapse.org/#note-for-windows-and-mac-users) to install synapser
- [Java = 21](https://www.java.com/en/download/)
    - For mac users, it seems to work better to run `brew install java`
- [wget](https://www.gnu.org/software/wget/)
    - For mac users, have to run `brew install wget`

## File Validator

One of the features of the `aacrgenie` package is that is provides a local validation tool that GENIE data contributors and install and use to validate their files locally prior to uploading to Synapse.


### Setting up your environment

These instructions will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.

1. Create a virtual environment using package manager of your choice (e.g: `conda`, `pipenv`, `pip`)

Example of creating a simple python environment

```
python3 -m venv <env_name>
source <env_name>/bin/activate
```

2. Install the genie package

```
pip install aacrgenie
```

3. Verify the installation

```
genie -v
```

4. Set up authentication with Synapse through the [local .synapseConfig](https://python-docs.synapse.org/tutorials/authentication/#use-synapseconfig) or using an [environment variable](https://python-docs.synapse.org/tutorials/authentication/#use-environment-variable)

### Running the validator

Get help of all available commands

```
genie validate -h
```

### Example commands

Running validator on clinical file

```
genie validate data_clinical_supp_SAGE.txt SAGE
```

Running validator on cna file. **Note** that the flag `--skip-database-checks` is **REQUIRED** when running the validator for cna and assay information files because you would need access to internal bed and clinical database tables respectively without it. For DEVELOPERS this is not required.

```
genie validate data_cna_SAGE.txt SAGE --skip-database-checks
```

Running validator on assay_information file.

```
genie validate assay_information.yaml SAGE --skip-database-checks
```

## Contributing

Please view [contributing guide](CONTRIBUTING.md) to learn how to contribute to the GENIE package.


# Sage Bionetworks Only

## Running locally

These are instructions on how you would setup your environment and run the pipeline locally.

1. Make sure you have read through the [GENIE Onboarding Docs](https://sagebionetworks.jira.com/wiki/spaces/APGD/pages/2163344270/Onboarding) and have access to all of the required repositories, resources and synapse projects for Main GENIE.
1. Be sure you are invited to the Synapse GENIE Admin team.
1. Make sure you are a Synapse certified user: [Certified User - Synapse User Account Types](https://help.synapse.org/docs/Synapse-User-Account-Types.2007072795.html#SynapseUserAccountTypes-CertifiedUser)
1. Be sure to clone the cbioportal repo: https://github.com/cBioPortal/cbioportal and `git checkout` the version of the repo pinned to the [Dockerfile](https://github.com/Sage-Bionetworks/Genie/blob/main/Dockerfile)
1. Be sure to clone the annotation-tools repo: https://github.com/Sage-Bionetworks/annotation-tools and `git checkout` the version of the repo pinned to the [Dockerfile](https://github.com/Sage-Bionetworks/Genie/blob/main/Dockerfile).

### Using `conda`

Follow instructions to install conda on your computer:

Install `conda-forge` and [`mamba`](https://github.com/mamba-org/mamba)
```
conda install -n base -c conda-forge mamba
```

Install Python and R versions via `mamba`
```
mamba create -n genie_dev -c conda-forge python=3.10 r-base=4.3  
```

### Using `pipenv`

Installing via [pipenv](https://pipenv.pypa.io/en/latest/installation.html)

1. Specify a python version that is supported by this repo:

    ```
    pipenv --python <python_version>
    ```

1. [pipenv install from requirements file](https://docs.pipenv.org/en/latest/advanced.html#importing-from-requirements-txt)

1. Activate your `pipenv`:

    ```
    pipenv shell
    ```

### Using `docker` (**HIGHLY** Recommended)

This is the most reproducible method even though it will be the most tedious to develop with. See [CONTRIBUTING docs for how to locally develop with docker.](/CONTRIBUTING.md). This will setup the docker image in your environment.

1. Pull pre-existing docker image or build from Dockerfile:
    Pull pre-existing docker image. You can find the list of images [from here.](https://github.com/Sage-Bionetworks/Genie/pkgs/container/genie)
    ```
    docker pull <some_docker_image_name>
    ```

    Build from Dockerfile
    ```
    docker build -f Dockerfile -t <some_docker_image_name> .
    ```

1. Run docker image:
    ```
    docker run --rm -it -e SYNAPSE_AUTH_TOKEN=$YOUR_SYNAPSE_TOKEN <some_docker_image_name>
    ```

### Setting up

1. Clone this repo and install the package locally.

    Install Python packages. This is the more traditional way of installing dependencies. Follow instructions [here](https://pip.pypa.io/en/stable/installation/) to learn how to install pip.

    ```
    pip install -e .
    pip install -r requirements.txt
    pip install -r requirements-dev.txt
    ```

    Install R packages. Note that the R package setup of this is the most unpredictable so it's likely you have to manually install specific packages first before the rest of it will install.
    ```
    Rscript R/install_packages.R
    ```

1. Configure the Synapse client to authenticate to Synapse.
    1. Create a Synapse [Personal Access token (PAT)](https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens).
    1. Add a `~/.synapseConfig` file
        ```
        [authentication]
        authtoken = <PAT here>
        ```
    1. OR set an environmental variable
        ```
        export SYNAPSE_AUTH_TOKEN=<PAT here>
        ```
    1. Confirm you can log in your terminal.
        ```shell
        synapse login
        ```

1. Run the different steps of the pipeline on the test project.  The `--project_id syn7208886` points to the test project. You should always be using the test project when developing, testing and running locally.

    1. Validate all the files **excluding vcf files**:

        ```
        python3 bin/input_to_database.py main --project_id syn7208886 --onlyValidate
        ```

    1. Validate **all** the files:

        ```
        python3 bin/input_to_database.py mutation --project_id syn7208886 --onlyValidate --genie_annotation_pkg ../annotation-tools
        ```

    1. Process all the files aside from the mutation (maf, vcf) files.  The mutation processing was split because it takes at least 2 days to process all the production mutation data.  Ideally, there is a parameter to exclude or include file types to process/validate, but that is not implemented.

        ```
        python3 bin/input_to_database.py main --project_id syn7208886 --deleteOld
        ```

    1. Process the mutation data. This command uses the `annotation-tools` repo that you cloned previously which houses the code that standardizes/merges the mutation (both maf and vcf) files and re-annotates the mutation data with genome nexus.  The `--createNewMafDatabase` will create a new mutation tables in the test project.  This flag is necessary for production data for two main reasons:
        * During processing of mutation data, the data is appended to the data, so without creating an empty table, there will be duplicated data uploaded.
        * By design, Synapse Tables were meant to be appended to.  When a Synapse Tables is updated, it takes time to index the table and return results. This can cause problems for the pipeline when trying to query the mutation table.  It is actually faster to create an entire new table than updating or deleting all rows and appending new rows when dealing with millions of rows.
        * If you run this more than once on the same day, you'll run into an issue with overwriting the narrow maf table as it already exists. Be sure to rename the current narrow maf database under `Tables` in the test synapse project and try again.

        ```
        python3 bin/input_to_database.py mutation --project_id syn7208886 --deleteOld --genie_annotation_pkg ../annotation-tools --createNewMafDatabase
        ```

    1. Create a consortium release.  Be sure to add the `--test` parameter. For consistency, the `processingDate` specified here should match the one used in the `consortium_map` for the `TEST` key [nf-genie.](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/main.nf)

        ```
        python3 bin/database_to_staging.py <processingDate> ../cbioportal TEST --test
        ```

    1. Create a public release.  Be sure to add the `--test` parameter. For consistency, the `processingDate` specified here should match the one used in the `public_map` for the `TEST` key [nf-genie.](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/main.nf)

        ```
        python3 bin/consortium_to_public.py <processingDate> ../cbioportal TEST --test
        ```

## Developing

1. Navigate to your cloned repository on your computer/server.
1. Make sure your `develop` branch is up to date with the `Sage-Bionetworks/Genie` `develop` branch.

    ```
    cd Genie
    git checkout develop
    git pull
    ```

1. Create a feature branch which off the `develop` branch. If there is a GitHub/JIRA issue that you are addressing, name the branch after the issue with some more detail (like `{GH|GEN}-123-add-some-new-feature`).

    ```
    git checkout -b GEN-123-new-feature
    ```

1. At this point, you have only created the branch locally, you need to push this remotely to Github.

    ```
    git push -u origin GEN-123-new-feature
    ```

1. Add your code changes and push them via useful commit message
    ```
    git add
    git commit changed_file.txt -m "Remove X parameter because it was unused"
    git push
    ```

1. Once you have completed all the steps above, in Github, create a pull request (PR) from your feature branch to the `develop` branch of Sage-Bionetworks/Genie.


### Developing with Docker

See [using `docker`](#using-docker-highly-recommended) for setting up the initial docker environment.

A docker build will be created for your feature branch every time you have an open PR on github and add the label `run_integration_tests` to it.

It is recommended to develop with docker. You can either write the code changes locally, push it to your remote and wait for docker to rebuild OR do the following:

1. Make any code changes. These cannot be dependency changes - those would require a docker rebuild.
1. Create a running docker container with the image that you pulled down or created earlier

    ```
    docker run -d <docker_image_name> /bin/bash -c "while true; do sleep 1; done"
    ```

1. Copy your code changes to the docker image:

    ```
    docker cp <folder or name of file> <docker_image_name>:/root/Genie/<folder or name of files>
    ```

1. Run your image in interactive mode:

    ```
    docker exec -it -e SYNAPSE_AUTH_TOKEN=$YOUR_SYNAPSE_TOKEN <docker_image_name> /bin/bash
    ```

1. Do any commands or tests you need to do

### Modifying Docker

Follow this section when modifying the [Dockerfile](https://github.com/Sage-Bionetworks/Genie/blob/main/Dockerfile):

1. Have your synapse authentication token handy
1. ```docker build -f Dockerfile -t <some_docker_image_name> .```
1. ```docker run --rm -it -e SYNAPSE_AUTH_TOKEN=$YOUR_SYNAPSE_TOKEN <some_docker_image_name>```
1. Run [test code](README.md#developing-locally) relevant to the dockerfile changes to make sure changes are present and working
1. Once changes are tested, follow [genie contributing guidelines](#developing) for adding it to the repo
1. Once deployed to main, make sure the CI/CD build successfully completed (our docker image gets automatically deployed via Github Actions CI/CD) [here](https://github.com/Sage-Bionetworks/Genie/actions/workflows/ci.yml)
1. Check that your docker image got successfully deployed [here](https://github.com/Sage-Bionetworks/Genie/pkgs/container/genie)


## Testing

Currently our Github Actions will run unit tests from our test suite `/tests` and run integration tests - each of the [pipeline steps here](README.md#developing-locally) on the test pipeline.

These are all triggered by adding the Github label `run_integration_tests` on your open PR.

To trigger `run_integration_tests`:

- Add `run_integration_tests` for the first time when you just open your PR
- Remove `run_integration_tests` label and re-add it
- Make any commit and pushes when the PR is still open

If you are developing with docker, docker images for your feature branch also gets build via the `run_integration_tests` trigger so check that your docker image got successfully deployed[here](https://github.com/Sage-Bionetworks/Genie/pkgs/container/genie).

### Running unit tests

Unit tests in Python are also run automatically by Github Actions on any PR and are required to pass before merging.

Otherwise, if you want to add tests and run tests outside of the CI/CD, see [how to run tests and general test development](./CONTRIBUTING.md#testing)

### Running integration tests

See [running pipeline steps here](README.md#developing-locally) if you want to run the integration tests locally.

You can also run them in nextflow via [nf-genie](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/README.md)


## Production

The production pipeline is run on Nextflow Tower and the Nextflow workflow is captured in [nf-genie](https://github.com/Sage-Bionetworks-Workflows/nf-genie).  It is wise to create an ec2 via the Sage Bionetworks service catalog to work with the production data, because there is limited PHI in GENIE.

## Github Workflows

For technical details about our CI/CD, please see [the github workflows README](.github/workflows/README.md)