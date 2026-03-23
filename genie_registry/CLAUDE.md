<!-- Last reviewed: 2026-03 -->

## Project

Pluggable file format registry for the GENIE pipeline. Each file defines a format class that inherits from `FileTypeFormat` and implements validation and processing logic for one genomic data format.

## Conventions

- Every class must set `_fileType` (string identifier), `_validation_kwargs`, and `_process_kwargs` class attributes.
- Override `_validateFilename()` using `assert` — the pipeline catches `AssertionError` to try the next format.
- Override `_validate()` to return `(error_string, warning_string)`. Empty error string = valid.
- Override `_cross_validate()` only if this format needs cross-file checks. Cross-validation for clinical data stays in `clinical.py`.
- Override `process_steps()` for processing logic. Call `preprocess()` for setup.
- Default file reading is TSV via `pd.read_csv(path, sep="\t", comment="#")`. Override `_get_dataframe()` only if format differs (e.g., YAML for assay).
- Class naming is inconsistent (lowercase `bed`, `maf`, `cna` vs PascalCase `Clinical`, `StructuralVariant`). Match the existing convention for the class you're modifying.

## Constraints

- Validation error messages must follow the format documented in root CLAUDE.md — because sites receive these messages directly.

## Do NOT

- **Do NOT replace `assert` statements in `_validateFilename()` with `raise ValueError`.** The pipeline catches `AssertionError` to try the next file format in the plugin registry. Changing to ValueError breaks format discovery.
- **Do NOT put cross-file validation in format classes other than `clinical.py`.** All cross-file checks go in `clinical.py._cross_validate()` — because clinical is the source of truth.
- **Do NOT rename `_fileType` values.** They are registry keys referenced throughout the pipeline and in Synapse table mappings.
- **Do NOT "fix" the `patientRetraction` → `sampleRetraction` inheritance chain.** `patientRetraction` intentionally inherits from `sampleRetraction`, not directly from `FileTypeFormat`.
- **Do NOT write custom column existence checks.** Use `process_functions.checkColExist()` — reviewer explicitly flagged this in PR #622.
- **Do NOT omit specific invalid values from error messages.** Messages must include "Found: <values>" — a revert happened when this was removed (sites need it for debugging).
- **Do NOT skip NA/empty string edge cases in validation tests.** Reviewer required NA test cases in PR #622 and PR #634. Always test: NA values, empty strings, and empty columns.
