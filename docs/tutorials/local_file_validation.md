# Local File Validation

One of the features of the `aacrgenie` package is that is provides a local validation tool that GENIE data contributors and install and use to validate their files locally prior to uploading to Synapse.

## Setting up your environment

These instructions will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.

1. Create a virtual environment using package manager of your choice (e.g: `conda`, `pipenv`, `pip`)

Example of creating a simple python environment

```bash
python3 -m venv <env_name>
source <env_name>/bin/activate
```

2. Install the genie package

```bash
pip install aacrgenie
```

3. Verify the installation

```bash
genie -v
```

4. Set up authentication with Synapse through the [local .synapseConfig](https://python-docs.synapse.org/tutorials/authentication/#use-synapseconfig) or using an [environment variable](https://python-docs.synapse.org/tutorials/authentication/#use-environment-variable)


This will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.

## Running the validator

Please view the help to see how to run the validator.

```
genie validate -h
```

Validate a file

```
genie validate data_clinical_supp_SAGE.txt SAGE
```

### Special Consideration

The flag `--skip-database-checks` is **REQUIRED** when running the validator for cna and assay information files because you would need access to internal bed and clinical database tables respectively without it. Without the flag, you will hit an Synapse `READ access` error.


#### Examples

Running validator on cna file.

```
genie validate data_cna_SAGE.txt SAGE --skip-database-checks
```

Running validator on assay_information file.

```
genie validate assay_information.yaml SAGE --skip-database-checks
```
