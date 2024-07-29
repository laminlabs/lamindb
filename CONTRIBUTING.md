# Contributing

Contributions are generally welcome. Please make an issue to discuss proposals.

For installation from PyPI, see [docs.lamin.ai/setup](https://docs.lamin.ai/setup).

For installation from GitHub, call:

```
git clone --recursive https://github.com/laminlabs/lamindb
pip install laminci
nox -s install
```

This will install a few dependencies from the git submodules linked [here](https://github.com/laminlabs/lamindb/tree/main/sub).

Running tests requires the Docker daemon up, then run at the root of the repository:

```
pytest
```

We use `pre-commit` and `gitmoji`. At the root of your repo, please call

```
gitmoji -i
pre-commit install
```
