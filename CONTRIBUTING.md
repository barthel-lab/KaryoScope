# Contributing to KaryoScope

Thank you for your interest in contributing to KaryoScope! Contributions of all kinds are welcome — bug reports, feature requests, documentation improvements, and code.

This document covers the practical details. For community norms, please also read our [Code of Conduct](CODE_OF_CONDUCT.md).

## Ways to contribute

- **Report bugs** by opening an [issue](https://github.com/barthel-lab/KaryoScope/issues/new?template=bug_report.md). Include a minimal reproduction if possible.
- **Suggest features** via the [feature request template](https://github.com/barthel-lab/KaryoScope/issues/new?template=feature_request.md).
- **Ask questions** via the [question template](https://github.com/barthel-lab/KaryoScope/issues/new?template=question.md) or [GitHub Discussions](https://github.com/barthel-lab/KaryoScope/discussions) once enabled.
- **Improve documentation** — typo fixes and clarifications are always welcome, no issue needed first.
- **Contribute code** — see below.

## Getting set up for development

KaryoScope requires Python ≥3.10 and several external command-line tools (`KMC`, `bgzip`, `tabix`, `seqtk`). We recommend using a dedicated [conda](https://docs.conda.io/) (or [mamba](https://mamba.readthedocs.io/)) environment for everything, which lets you manage Python and the bioinformatics tools together.

```bash
git clone git@github.com:barthel-lab/KaryoScope.git
cd KaryoScope

# Create a dedicated environment with Python and the required external tools
conda create -n karyoscope-dev -c conda-forge -c bioconda \
    python=3.12 pip kmc htslib seqtk
conda activate karyoscope-dev

# Install KaryoScope in editable mode with dev dependencies
pip install -e ".[dev]"

# Verify the install
karyoscope --version
```

If you prefer to manage Python separately (e.g., a system `venv` with conda only providing the external tools), that works too, but you'll need to activate both environments. The single-conda-env approach above is simpler.

> **Note:** Do not install KaryoScope or its bioinformatics dependencies into your `base` conda environment. Always use a dedicated env. This avoids conflicts with other tools and keeps your setup reproducible.

## Running tests and linters

```bash
# Run the full test suite
pytest

# Skip slow / integration tests during local development
pytest -m "not slow and not integration"

# Lint and format
ruff check .
ruff format --check .

# Auto-fix what can be auto-fixed
ruff check --fix .
ruff format .
```

CI runs all of the above on every push and pull request.

## Contributing code

### Before you start

For non-trivial changes (anything more than a small bug fix or doc tweak), please open an issue first to discuss the approach. This avoids wasted work if the change isn't a fit, and lets us point you at relevant context.

### Branching and pull requests

1. Fork the repository and create a branch off `main`.
2. Use a descriptive branch name: `fix/issue-42`, `feature/scaffold-improvements`, etc.
3. Make your changes. Keep commits focused and write clear messages.
4. Run linters and tests locally before pushing.
5. Open a pull request against `main`. Reference any related issues.
6. CI will run automatically. Address any failures.
7. A maintainer will review. We aim to respond within a week.

### Style

- Code is formatted by `ruff format` (Black-compatible style, 100-char lines). CI enforces this.
- Lint with `ruff check`. Configured rules are in `pyproject.toml`.
- Type annotations are encouraged on public APIs but not currently required throughout.
- Docstrings: short and direct. Document the *why* in comments; let function signatures and names communicate the *what*.

### Tests

- New features should include tests. Tests live in `tests/`.
- Unit tests should be fast (<1s each). Mark slow tests with `@pytest.mark.slow`.
- Tests that shell out to external tools should be marked `@pytest.mark.integration`.

### Commit messages

We don't enforce a strict format, but we appreciate:

- A concise subject line (≤72 chars), imperative mood ("Fix bug" not "Fixed bug").
- A blank line, then a body explaining the *why* if non-obvious.
- References to issues (`Fixes #42`, `Refs #17`) when relevant.

## License of contributions

KaryoScope is currently licensed under [GPL-3.0-or-later](LICENSE). By submitting code, you agree that your contributions will be licensed under the same terms.

> **Important note on a future license change:** KaryoScope is GPL-3.0 because it bundles the GPL-3.0 KMC API for k-mer indexing. We plan to relicense KaryoScope under the MIT license once the KMC dependency is replaced with the [HKS](https://github.com/jnalanko/HKS) library (which has a more permissive license). By contributing, you also consent to having your contribution relicensed under MIT at that time. If you are unable or unwilling to consent to this future relicensing, please mention it in your pull request, and we will discuss how to incorporate your work compatibly.

## Reporting security issues

Please **do not** report security vulnerabilities through public GitHub issues. Instead, email [rranallo-benavidez@tgen.org](mailto:rranallo-benavidez@tgen.org) and we will respond within 7 days.

## Maintainers

- Timothy Rhyker Ranallo-Benavidez ([@tbenavi1](https://github.com/tbenavi1))
- Floris P. Barthel ([@fpbarthel](https://github.com/fpbarthel))

## Thank you!

KaryoScope is a research tool built and maintained alongside other lab work. Contributions — large or small — directly help the genomics community. We genuinely appreciate your time.
