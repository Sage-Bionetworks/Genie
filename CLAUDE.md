<!-- Last reviewed: 2026-03 -->

## Project

AACR Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange) ETL pipeline.
Validates, processes, and releases genomic cancer data submitted by participating centers.
Data flows: center files → validation → Synapse tables → consortium release → public release.
Published as `aacrgenie` on PyPI. Jira project: GEN.

## Stack

- Python >=3.10, <3.12 (primary)
- R 4.3.3 (analytics, renv-managed)
- Synapse (data repository and table backend via `synapseclient`)
- Docker (ubuntu:jammy base with Python, R, Java 21, cBioPortal, bedtools)
- pytest (testing), black (formatting), ruff + flake8 (linting)
- GitHub Actions CI/CD, published to ghcr.io and PyPI

## Commands

```bash
# Run tests with coverage
pytest tests/ --cov=genie --cov=genie_registry --cov-report=html

# Format code
black ./

# Build package
python -m build

# Install for development
pip install -e .

# Local file validation
genie validate <file> <center> [--filetype TYPE] [--skip-database-checks]
```

## Data Models

The pipeline uses a plugin registry pattern for file formats:

- **`FileTypeFormat`** (`genie/example_filetype_format.py`) — abstract base class (ABCMeta). All file formats inherit from it. Defines the contract: `_validateFilename()`, `_get_dataframe()`, `_validate()`, `_cross_validate()`, `process_steps()`.
- **`ValidationResults`** (dataclass) — holds `errors: str`, `warnings: str`, `detailed: Optional[str]`. Empty `errors` string = valid.
- **Format classes** in `genie_registry/` — each sets `_fileType`, `_validation_kwargs`, `_process_kwargs`. Discovered at runtime via `genie/config.py` using `__subclasses__()` reflection.
- **Config registry** (`genie/config.py`) — `collect_format_types(package_names)` imports packages and finds all `FileTypeFormat` subclasses. Extensible via `--format_registry_packages` CLI arg.

Data flows through pandas DataFrames. All column names are uppercased on read. Files are TSV with `comment="#"` support. Synapse tables store processed data via `synapseclient.tableQuery()`.

## Conventions

- **Validation error message format**: `<filename>: <description>. <valid values if applicable>. You have N row(s) where <rule fails>. The row(s) this occurs in are: <indices>. Please correct.` — per Confluence SOP.
- **Cross-file validation always goes in `clinical.py`** `_cross_validate()` — because clinical is the source of truth for cross-file checks.
- **Clinical files come as pairs** — two files merged together. See the HACK comment in `genie/extract.py` around line 64.
- **Validation functions need detailed docstrings** with examples of passing and failing data — this is a Confluence-documented requirement.
- **Branch naming**: `gen-XXXX-short-description` (lowercase, Jira ticket prefix).
- **Commit messages**: `[GEN-XXXX] Brief description (#PR_number)`.
- **PR template**: Problem / Solution / Testing sections. Add label `run_integration_tests` to trigger integration test CI on feature branches.
- **Test data in Synapse**: GOLD center files must PASS validation. TEST/SAGE center files should FAIL for new rules. Update Synapse files with version comments.
- **`check_col_and_values()` returns `(warning, error)`** — reversed from all other validation functions that return `(errors, warnings)`. Do not "fix" this without updating all call sites.

## Architecture

```
bin/ (CLI entry scripts)
  → genie/ (core ETL: validate, extract, transform, load, process_mutation)
    → genie_registry/ (pluggable file format classes)
      → Synapse (data storage via synapseclient)

Pipeline stages:
1. input_to_database — validate center files, process valid ones, load to Synapse tables
2. database_to_staging — consortium release: redact PHI, filter germline, format for cBioPortal
3. consortium_to_public — public release: filter to approved samples
```

Plugin discovery: `genie/config.py` uses `importlib.import_module()` + `__subclasses__()` to find all `FileTypeFormat` subclasses at runtime. Add new formats by creating a class in a new package — no core changes needed.

## Constraints

- **Never expose PHI.** Ages >89 are masked, pediatric <18 redacted, germline variants filtered (gnomAD AF >0.01). This is a compliance requirement — because GENIE handles real patient genomic data.
- **New validation rules require 6-month advance notice to sites before release.** Exceptions: urgent, <3 sites affected, or low-complexity fix. This is documented in Confluence SOP.
- **Synapse auth via `SYNAPSE_AUTH_TOKEN` env var.** Never commit tokens. The `.synapseConfig` file is an alternative for local dev.
- **Version lives in `genie/__init__.py`** — single source of truth. `setup.cfg` reads it via `attr:`. Update only there.
- **Gitflow workflow.** All features merge to `develop` first. Release branches merge to `main`. Never push directly to `main` except doc/CI patches.

## Testing

- pytest with session-scoped fixtures for `syn` (mocked Synapse client) and `genie_config` (config dict).
- Mock pattern: `from unittest.mock import patch` then `patch.object(MODULE, "FUNC", return_value=...)`.
- Test parametrization via `@pytest.mark.parametrize` — include ALL pass/fail cases for validation rules.
- Integration tests triggered by `run_integration_tests` label on PRs. Run against test Synapse project `syn7208886`.
- R tests via `testthat`: `Rscript -e "testthat::test_dir('R/tests/testthat/')"`.

## Related Systems

- **nf-genie** (Sage-Bionetworks-Workflows/nf-genie) — Nextflow pipeline that orchestrates GENIE runs on Seqera Platform (Nextflow Tower).
- **annotation-tools** (Sage-Bionetworks/annotation-tools) — Genome Nexus annotation wrapper. Built as `annotator.jar` (Java 21). Used for MAF/VCF re-annotation.
- **Synapse** (`synapseclient`) — all data stored in Synapse projects. Production: `syn3380222`. Test: `syn7208886`.
- **cBioPortal** (v5.3.19) — output format target for releases. Case lists and gene matrices generated for it.
- **Confluence space APGD** — "AACR Project GENIE Documentation". Contains SOPs, validation rule docs, troubleshooting.
- **Jira project GEN** — all tickets prefixed `GEN-`. Board at sagebionetworks.jira.com.
- **R scripts** (`R/`) use the same three hardcoded dbMapping Synapse IDs as `bin/` scripts (test: `syn11600968`, staging: `syn12094210`, prod: `syn10967259`). When updating IDs, grep across both Python and R.

## Anti-Patterns — Do NOT

- **Do NOT simplify validation error messages.** Error messages must include the specific invalid values found (e.g., "Found: 2, 3"). A revert (7135f18) happened when error detail was removed — because sites use these messages to debug their data submissions.
- **Do NOT deprecate a feature without full cleanup in the same PR.** When deprecating validation/processing for a field, remove ALL related test cases, docstrings, and processing functions together (PR #638/#639 lesson).
- **Do NOT refactor `DataFrame.append()` to `pd.concat()` without extensive testing.** A previous attempt was reverted (da18b5f) because it broke subtle dataframe ordering behavior in dashboard_table_updater.py.
- **Do NOT use bare `set()` for DataFrame column/index construction.** Use `sorted(set(...))` — because set ordering is non-deterministic and broke tests (PR #621, GEN-2377).
- **Do NOT hardcode Synapse IDs inline.** They are already scattered across Python and R files. When adding new ID references, add them near existing ones and document the duplication.
- **Do NOT declare Python version support without CI coverage.** Python 3.12 is in setup.cfg range but fails CI — this was flagged in PR #559 review.
- **Do NOT add dependency version upper bounds without documenting why.** Reviewer asked "Is there a reason <4.0.0?" in PR #559 — unclear bounds create unnecessary conflicts.
