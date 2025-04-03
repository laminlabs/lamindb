# Contributing

Contributions are generally welcome. Please make an issue to discuss proposals.

## Installation

### PyPI

For installation from PyPI, see [docs.lamin.ai/setup](https://docs.lamin.ai/setup).

### Github

For installation from GitHub, call:

```bash
git clone --recursive https://github.com/laminlabs/lamindb
pip install laminci
python -m venv .venv
source .venv/bin/activate
nox -s install
```

This will install a few dependencies from the git submodules linked [here](https://github.com/laminlabs/lamindb/tree/main/sub), as well as packages
like `pytest` and `pre-commit` that you'll need when developing.

lamindb depends on several other packages that may require modifications for pull requests to successfully pass the continuous integration build.
We suggest the following workflow if commits to any of the submodules are essential for the current modifications in lamindb:

1. Change directory into the submodule that you want to modify: `cd sub/SUBMODULE`.
2. Switch to a new feature branch: `git switch -c feature/NEWFEATURE`.
3. Make a pull request with your changes to the `SUBMODULE` and ensure that the CI passes.
4. In the repository root of lamindb, create a new commit and push:

```bash
cd ..
git add -u
git commit -m "Upgraded SUBMODULE"
git push
```

Any pull request of yours should now also have the changes of the submodule included allowing you to test that changes in the submodule and lamindb are compatible.

## Running and writing tests

This package uses the [pytest][] for automated testing.
Please add a test for every function added to the package.

Running tests requires the [Docker daemon][] up, then run at the root of the repository:

```bash
pytest --ignore=tests/storage --ignore=tests/permission
```

in the root of the repository.
We exclude specific directories in local `pytest` runs because they directly access external resources such as AWS, which require specific access keys.
Continuous integration will automatically run **all** tests on pull requests.

## Code-style

This project uses [pre-commit][] to enforce consistent code-styles. On every commit, pre-commit checks will either
automatically fix issues with the code, or raise an error message.

To enable pre-commit locally, simply run

```bash
pre-commit install
```

in the root of the repository. Pre-commit will automatically download all dependencies when it is run for the first time.

We further use [gitmoji][] to add emoticons to commits.
These allow us to more easily categorize them allowing for faster visual filtering.

It can be installed by running:

```bash
npm i -g gitmoji-cli
```

and enabled for the repository via:

```bash
gitmoji -i
```

If you don't have `sudo` in your working environment, follow [these instructions](https://github.com/sindresorhus/guides/blob/main/npm-global-without-sudo.md).

## Documentation

We build our documentation with an internal tool called `lndocs`.
We have not made it public yet and therefore external contributors need to rely on the Github Actions `docs` job to build the documentation.
If the `docs` job succeeds, a preview URL will be posted automatically as a comment to your pull request.

## Releases

Currently only lamin employees have release rights.

[Docker daemon]: https://docs.docker.com/engine/install/
[gitmoji]: https://gitmoji.dev/
[pre-commit]: https://pre-commit.com/
[pytest]: https://docs.pytest.org/
