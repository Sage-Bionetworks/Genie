
## Contributing

We welcome all contributions!  Please head to [issues](https://github.com/Sage-Bionetworks/challengeutils/issues) to either file any bugs/feature requests or find a task you want to assist with.  Make sure to assign yourself the task if you decide to work on it.


### Fork and clone this repository

See the [Github docs](https://help.github.com/articles/fork-a-repo/) for how to make a copy (a fork) of a repository to your own Github account.

Then, [clone the repository](https://help.github.com/articles/cloning-a-repository/) to your local machine so you can begin making changes.

Add this repository as an [upstream remote](https://help.github.com/en/articles/configuring-a-remote-for-a-fork) on your local git repository so that you are able to fetch the latest commits.

On your local machine make sure you have the latest version of the `develop` branch:

```
git checkout develop
git pull upstream develop
```

### Install development dependencies
This will install all the dependencies of the package including the active branch of `aacrgenie`.

```
pip install -r requirements-dev.txt
```


### The development life cycle

The GENIE project follows the standard [git flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) development strategy.

1. Pull the latest content from the `develop` branch of this central repository (not your fork).
1. Create a feature branch which off the `develop` branch. If there is a GitHub issue that you are addressing, name the branch after the issue with some more detail (like `issue-123-add-some-new-feature`).
1. After completing work and testing locally (see below), push to your fork.
1. In Github, create a pull request from the feature branch of your fork to the `develop` branch of the central repository.

> *A code maintainer must review and accept your pull request.* A code review (which happens with both the contributor and the reviewer present) is required for contributing. This can be performed remotely (e.g., Skype, Hangout, or other video or phone conference).

This package uses [semantic versioning](https://semver.org/) for releasing new versions. The version should be updated on the `develop` branch as changes are reviewed and merged in by a code maintainer. The version for the package is maintained in the [genie/__version__.py](genie/__version__.py) file.  A github release should also occur every time `develop` is pushed into `master` and it should match the version for the package.

### Testing

Please add tests for new code. These might include unit tests (to test specific functionality of code that was added to support fixing the bug or feature), integration tests (to test that the feature is usable - e.g., it should have complete the expected behavior as reported in the feature request or bug report), or both.

This package uses [`pytest`](https://pytest.org/en/latest/) to run tests. The test code is located in the [test](./test) subdirectory.

Here's how to run the test suite:

```
pytest -vs tests/
```

Tests are also run automatically by Github Actions on any pull request and are required to pass before merging.


### Release Procedure (For Package Maintainers)

Follow gitflow best practices as linked above.

* Always merge all new features into `develop` branch first (unless it is a documentation, readme, or github action patch into `main`)
* After initial features are ready in the `develop` branch, create a `release-X.X` branch to prepare for the release.
    * update `genie/__version__.py`
* Merge `release-X.X` branch into `main`
* Create release tag (`v...`) and include release notes.  Also include any known bugs for each release here.
* Merge `master` back into `develop`
