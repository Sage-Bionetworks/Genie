# Getting Started

## Dependencies

This package contains both R, Python and cli tools.  These are tools or packages you will need, to be able to reproduce these results:

- Python >=3.8 or <3.10
    - `pip install -r requirements.txt`
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- R 4.2.2
    - `renv::install()`
    - Follow instructions [here](https://r-docs.synapse.org/#note-for-windows-and-mac-users) to install synapser
- [Java > 8](https://www.java.com/en/download/)
    - For mac users, it seems to work better to run `brew install java`
- [wget](https://www.gnu.org/software/wget/)
    - For mac users, have to run `brew install wget`

## Installation

### PyPi

The [aacrgenie](https://pypi.org/project/aacrgenie/) package is available from PyPI. It can be installed or upgraded with pip.

```
pip install aacrgenie
```

### Local

Source code and development versions are available on Github. Installing from source:

```
git clone https://github.com/Sage-Bionetworks/Genie.git
cd Genie
```

Install the packages locally.

```
pip install -e .
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

You can stay on the main branch to get the latest stable release or check out the develop branch or a tagged revision:

```
git checkout <branch or tag>
```


## Configuration

Configure the Synapse client to authenticate to Synapse.

1. Create a Synapse [Personal Access token (PAT)](https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens).
1. Add a `~/.synapseConfig` file with the following information:

    ```
    [authentication]
    authtoken = <PAT here>
    ```

1. OR set an environmental variable

    ```
    export SYNAPSE_AUTH_TOKEN=<PAT here>
    ```

1. Confirm you can log in to synapse in your terminal.

    ```shell
    synapse login
    ```