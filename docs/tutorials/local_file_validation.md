# Local File Validation

One of the features of the `aacrgenie` package is that is provides a local validation tool that GENIE data contributors and install and use to validate their files locally prior to uploading to Synapse.

```
pip install aacrgenie
genie -v
```

This will install all the necessary components for you to run the validator locally on all of your files, including the Synapse client.  Please view the help to see how to run to validator.

```
genie validate -h
```

Validate a file

```
genie validate data_clinical_supp_SAGE.txt SAGE
```
