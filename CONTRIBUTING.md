# Contributing

Contributions are generally welcome. Please make an issue to discuss proposals.

## Code-style

Enable `pre-commit` and `gitmoji` in the root of this repo:

```bash
pre-commit install
gitmoji -i
```

## Coverage

All test suites only run on the `release` and `coverage` branches. Only on these branches you'll recapitulate the full code coverage of the combined suites.

## Releases

Only maintainers have release rights.

Release publishing is managed via `laminci release --pypi`. For `lamindb`, the
release flow publishes two distributions:

- `lamindb-core`: installs the `lamindb/` namespace package
- `lamindb`: installs additional recommended dependencies on top of `lamindb-core`

Start the release flow by pushing to the `release` branch, which will trigger integration tests. Once you're happy, make the release.
