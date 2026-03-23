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

- Do not rename `_fileType` values — they are used as registry keys and referenced throughout the pipeline and in Synapse table mappings.
- `patientRetraction` inherits from `sampleRetraction`, not directly from `FileTypeFormat`. Preserve this chain.
- Validation error messages must follow the format documented in root CLAUDE.md — because sites receive these messages directly.
