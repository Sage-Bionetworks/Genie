
## Contributing

We welcome all contributions!  Please head to [issues](https://github.com/Sage-Bionetworks/Genie/issues) to either file any bugs/feature requests or find a task you want to assist with.  Make sure to assign yourself the task if you decide to work on it.

## Coding Style

This package uses `flake8` - it's settings are described in [setup.cfg](setup.cfg).  The code in this package is also automatically formatted by `black` for consistency.

## The Development Life Cycle

### Fork and clone this repository

1. See the [Github docs](https://help.github.com/articles/fork-a-repo/) for how to make a copy (a fork) of a repository to your own Github account.
1. Then, [clone the repository](https://help.github.com/articles/cloning-a-repository/) to your local machine so you can begin making changes.
1. Add this repository as an [upstream remote](https://help.github.com/en/articles/configuring-a-remote-for-a-fork) on your local git repository so that you are able to fetch the latest commits.
1. On your local machine make sure you have the latest version of the `develop` branch:

    ```
    git checkout develop
    git pull upstream develop
    ```

### Install development dependencies

This will install all the dependencies of the package including the active branch of `Genie`.  We highly recommend that you leverage some form of python version management like [pyenv](https://github.com/pyenv/pyenv) or [anaconda](https://www.anaconda.com/products/individual). There are two ways you can install the dependencies for this package.

#### pip
This is the more traditional way of installing dependencies. Follow instructions [here](https://pip.pypa.io/en/stable/installation/) to learn how to install pip.

```
pip install -r requirements-dev.txt
pip install -r requirements.txt
```

#### pipenv
`pipenv` is a Python package manager.  Learn more about [pipenv](https://pipenv.pypa.io/en/latest/) and how to install it.

```
# Coming soon
```

### Developing

The GENIE project follows the standard [git flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) development strategy.
> To ensure the most fluid development, try not to push to your `develop` or `main` branch.

1. (Assuming you have followed all 4 steps above in the "fork and clone this repository" section). Navigate to your cloned repository on your computer/server.
1. Make sure your `develop` branch is up to date with the `Sage-Bionetworks/Genie` `develop` branch.

    ```
    cd {your-github-username}/Genie
    git checkout develop
    git pull upstream develop
    ```

1. Create a feature branch which off the `develop` branch. If there is a GitHub/JIRA issue that you are addressing, name the branch after the issue with some more detail (like `{GH|JIRA}-123-add-some-new-feature`).

    ```
    git checkout develop
    git checkout -b JIRA-123-new-feature
    ```

1. At this point, you have only created the branch locally, you need to push this to your fork on GitHub.

    ```
    git push --set-upstream origin JIRA-123-new-feature
    ```

    You should now be able to see the branch on GitHub. Make commits as you deem necessary. It helps to provide useful commit messages - a commit message saying 'Update' is a lot less helpful than saying 'Remove X parameter because it was unused'.

    ```
    git commit changed_file.txt -m "Remove X parameter because it was unused"
    git push
    ```

1. (Make sure you have follow instructions in "Install development dependencies") Once you have made your additions or changes, make sure you write tests and run the test suite.  More information on testing below.

    ```
    pytest ./test
    ```

1. (Make sure you have follow instructions in "Install development dependencies") Make sure to run the auto python code formatter, black.

    ```
    black ./
    ```

1. Once you have completed all the steps above, in Github, create a pull request from the feature branch of your fork to the `develop` branch of Sage-Bionetworks/Genie.

> *A code maintainer must review and accept your pull request.* A code review ideally happens with both the contributor and the reviewer present, but is not strictly required for contributing. This can be performed remotely (e.g., Zoom, Hangout, or other video or phone conference).

This package uses [semantic versioning](https://semver.org/) for releasing new versions. The version should be updated on the `develop` branch as changes are reviewed and merged in by a code maintainer. The version for the package is maintained in the [genie/__init__.py](genie/__init__.py) file.  A github release should also occur every time `develop` is pushed into `main` and it should match the version for the package.

### Testing

Testing in GENIE is split between unit tests usually written in Python and R and integration tests.

#### Unit tests

##### Tests in Python

This package uses [`pytest`](https://pytest.org/en/latest/) to run tests. The test code is located in the [tests](./tests) subdirectory.

Here's how to run the test suite:

```shell
pytest -vs tests/
```

Tests in Python are also run automatically by Github Actions on any pull request and are required to pass before merging.

##### Tests in R

This package uses [`testthat`](https://testthat.r-lib.org/) to run tests in R. The test code is located in the [testthat](./R/tests/testthat) subdirectory.

Here's how to run the test suite:

```shell
Rscript -e "testthat::test_dir('R/tests/testthat/')"
```

##### Test Development

Please add tests for new code. These might include unit tests (to test specific functionality of code that was added to support fixing the bug or feature), integration tests (to test that the feature is usable - e.g., it should have complete the expected behavior as reported in the feature request or bug report), or both.

###### Mock Testing

It is recommended to use the following style (see example below) for mock testing across this package:

```python
from unittest.mock import patch
...
patch.object(MODULE_NAME, "FUNCTION_TO_MOCK_NAME".return_value=SOME_RETURN_VALUE)
```

#### Integration tests

Integration tests in Genie involve running all parts of the **test** pipeline under the following conditions:

- locally. Run each of the [pipeline steps here](README.md#developing-locally) on the test pipeline and verify that your pipeline runs as expected.
- in an ec2 instance using nextflow. [See here for how to run nextflow locally](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/README.md#process-and-developing-locally)
- in Nextflow Tower. [See here for how to run on Nextflow Tower.](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/README.md#processing-on-nextflow-tower)

These are __not__ automatically run by Github Actions and have to be manually run.

#### Modifying Docker

Follow this section when modifying the [Dockerfile](https://github.com/Sage-Bionetworks/Genie/blob/main/Dockerfile):

There are usually three common conditions when docker needs to be updated:

- The tag version of the `cbioportal` or `annotation-tools` repos needs to be updated (most common)
- The python version needs to be updated
- The java version needs to be updated

##### Developing and testing docker updates

1. Have your synapse authentication token handy
1. Make relevant changes to `Dockerfile`
1. ```docker build -f Dockerfile -t <some_docker_image_name> .```
1. ```docker run --rm -it -e SYNAPSE_AUTH_TOKEN=$YOUR_SYNAPSE_TOKEN <some_docker_image_name>```
1. Run [test code](README.md#developing-locally) relevant to the dockerfile changes to make sure changes are present and working
1. Push image to dockerhub under your branch of this repo
1. Follow instructions [here](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/README.md#process-and-developing-locally) for setting up your environment for running the nextflow workflow on [`nf-genie` repo](https://github.com/Sage-Bionetworks-Workflows/nf-genie).
1. Create new branch off the `main` branch in `nf-genie` repo
1. Update relevant module docker paths in `nf-genie` to specify your docker image (use `docker images` to see list)
1. Run `nf-genie` on ec2 under that docker image. [See here for how to run nextflow locally.](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/README.md#process-and-developing-locally)
1. Push your updates to your branch in `nf-genie`
1. Run the test pipeline on Nextflow Tower, be sure to specify your branch in the `nf-genie` pipeline as the `revision` when running. [See here for how to run on Nextflow Tower.](https://github.com/Sage-Bionetworks-Workflows/nf-genie/blob/main/README.md#processing-on-nextflow-tower)
1. Be sure to delete your branch on `nf-genie`
1. Once changes are tested, follow [genie contributing guidelines](#developing) for adding the docker updates to this repo
1. Once deployed to main, make sure docker image was successfully deployed remotely (our docker image gets automatically deployed) [here](https://hub.docker.com/repository/docker/sagebionetworks/genie/builds)


### Release Procedure (For Package Maintainers)

Follow gitflow best practices as linked above.

1. Always merge all new features into `develop` branch first (unless it is a documentation, readme, or github action patch into `main`)
1. After initial features are ready in the `develop` branch, create a `release-X.X` branch (do not need to push this branch to remote) to prepare for the release.
    1. update the `__version__` parameter in `genie/__init__.py`
1. Merge `release-X.X` branch into `main` - Not by pull request!
1. Create release tag (`v...`) and a brief message
1. Push tag and change(s) from `main`
1. Create a new release on the repo. Include release notes.  Also include any known bugs for each release here. Wait for the CI/CD to finish.
1. Merge `main` back into `develop`
1. Push `develop`

### Docker

This repository does not use github actions to push docker images.  By adding the `sagebiodockerhub` github user as an Admin to this GitHub repository, we can configure an automated build in DockerHub.  You can view the builds [here](https://hub.docker.com/repository/docker/sagebionetworks/genie/builds).  To get admin access to the DockerHub repository, ask Sage IT to be added to the `genieadmin` DockerHub team.
