
## Contributing

We welcome all contributions!  Please head to [issues](https://github.com/Sage-Bionetworks/Genie/issues) to either file any bugs/feature requests or find a task you want to assist with.  Make sure to assign yourself the task if you decide to work on it.

## Table of Contents

- [Coding Style](#coding-style)
- [The Development Life Cycle](#the-development-life-cycle)
  - [Fork and clone this repository](#fork-and-clone-this-repository)
  - [Install development dependencies](#install-development-dependencies)
  - [Developing](#developing)
  - [Testing](#testing)
    - [Running tests](#running-tests)
      - [Tests in Python](#tests-in-python)
      - [Tests in R](#tests-in-r)
    - [Test Development](#test-development)
      - [Mock Testing](#mock-testing)
  - [Release Procedure (For Package Maintainers)](#release-procedure-for-package-maintainers)
  - [Contributing to the docs](#contributing-to-the-docs)

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

This will install all the dependencies of the package including the active branch of `Genie`.  We highly recommend that you leverage some form of python version management like [pyenv](https://github.com/pyenv/pyenv) or [anaconda](https://www.anaconda.com/products/individual). Follow [dependencies installation instruction here](./README.md#running-locally)

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

1. Create a feature branch which off the `develop` branch. If there is a GitHub/JIRA issue that you are addressing, name the branch after the issue with some more detail (like `{GH|GEN}-123-add-some-new-feature`).

    ```
    git checkout develop
    git checkout -b GEN-123-new-feature
    ```

1. At this point, you have only created the branch locally, you need to push this remotely to Github.

    ```
    git push
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

#### Running tests

##### Tests in Python

This package uses [`pytest`](https://pytest.org/en/latest/) to run tests. The test code is located in the [tests](./tests) subdirectory.

Here's how to run the test suite:

```shell
pytest -vs tests/
```

##### Tests in R

This package uses [`testthat`](https://testthat.r-lib.org/) to run tests in R. The test code is located in the [testthat](./R/tests/testthat) subdirectory.

Here's how to run the test suite:

```shell
Rscript -e "testthat::test_dir('R/tests/testthat/')"
```

#### Test Development

Please add tests for new code. These might include unit tests (to test specific functionality of code that was added to support fixing the bug or feature), integration tests (to test that the feature is usable - e.g., it should have complete the expected behavior as reported in the feature request or bug report), or both.

##### Mock Testing

It is recommended to use the following style (see example below) for mock testing across this package:

```python
from unittest.mock import patch
...
patch.object(MODULE_NAME, "FUNCTION_TO_MOCK_NAME".return_value=SOME_RETURN_VALUE)
```

### Release Procedure (For Package Maintainers)

Follow gitflow best practices as linked above.

1. Always merge all new features into `develop` branch first (unless it is a documentation, readme, or github action patch into `main`)
2. After initial features are ready in the `develop` branch, create a `release-X.X.X` branch (**do not** push this branch to remote) to prepare for the release.
3. Update the `__version__` parameter in `genie/__init__.py`. Create a commit with a message like `"Release <X.X.X>"`
4. Merge `release-X.X.X` branch into `main` - Not by pull request! Resolve any merge conflicts.
5. Prior to pushing changes into `main`, reset any [annotations on the test project](https://www.synapse.org/Synapse:syn11601248/tables/)/[rename maf database table](https://www.synapse.org/Synapse:syn7208886/tables/) as integration tests will run.
6. Push change(s) from `main`.
7. You can do the following either through git or through the Github UI
   1. Create release tag (`v...`) and a brief message
  
    ```bash
    git tag -a v<X.X.X> -m "Release <X.X.X>"
    ```
   2. Push tag from `main` to remote.

   OR

   1. Create your tag via the release notes UI on github by creating the tag name `vX.X.X`. **Be sure to set `target branch` as `main`!**
8. Create a new release on the repo (if you created your tag through git prior, select your tag here). Generate Github release notes. See [writing release notes for project](https://sagebionetworks.jira.com/wiki/spaces/DPE/pages/3892543489/Writing+release+notes+for+projects) for creating informative technical release notes using the Github auto-generated release notes. **RECOMMENDED:** Use the [genie specific release notes prompt here.](https://sagebionetworks.jira.com/wiki/spaces/DPE/pages/3892543489/Writing+release+notes+for+projects#Genie-project-specific-Release-Notes-Prompt)
    - See [example release](https://github.com/Sage-Bionetworks/Genie/releases/tag/v16.5.0) for typical sections to include in a Genie release.
9. Wait for the CI/CD to finish. You should see a new pypi release here: https://pypi.org/project/aacrgenie/.
10. Reset any [annotations on the test project](https://www.synapse.org/Synapse:syn11601248/tables/)/[rename maf database table](https://www.synapse.org/Synapse:syn7208886/tables/) as integration tests will run.
11. Merge `main` changes back into `develop` locally
12. Push changes to `develop`
13. Wait for the CI/CD to finish.

### Contributing to the docs

This [documentation](https://sagebionetworks.jira.com/wiki/spaces/APGD/pages/3369631808/Contributing+to+Main+GENIE+repository+docs) is internal to Sage employees.
