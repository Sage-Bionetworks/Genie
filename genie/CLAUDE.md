<!-- Last reviewed: 2026-03 -->

## Project

Core ETL pipeline modules for GENIE. Handles validation, extraction, transformation, loading, and mutation processing. All file format logic lives in `genie_registry/` — this package provides the framework and orchestration.

## Template Method Pattern

- `validate()` in `example_filetype_format.py` calls `_validate()` first. Cross-validation via `_cross_validate()` runs ONLY if `_validate()` returned no errors AND `self.ancillary_files` is a non-empty list. Do not assume `_cross_validate()` always runs.
- `process()` reads the file differently by `_fileType`: clinical receives the raw `filePath` list, vcf/maf/mafSP/md receive the filepath string directly (no read), all others receive `read_file([filePath])` with path wrapped in a list. Match the convention for your file type.
- `_validation_kwargs` and `_process_kwargs` are enforced via `assert required_parameter in kwargs.keys()`. Missing kwargs crash with AssertionError — not a graceful error. Always register required kwargs in the class attribute.

## Return Conventions

- `_validate()` and `_cross_validate()` return `(errors: str, warnings: str)`. Empty string = valid. NOT booleans.
- `check_col_and_values()` in `process_functions.py` returns `(warning, error)` — **reversed order** from all other validation functions. Always destructure carefully.
- `checkColExist()` returns `bool`, not a tuple.
- `ValidationResults.is_valid()` checks `errors == ""` — use this method, not truthiness checks.

## Reusable Utilities — Do Not Reinvent

Use these existing functions instead of writing new ones:

- `process_functions.checkColExist(df, col)` — column existence check.
- `process_functions.check_col_and_values(df, col, possible_values, filename, ...)` — column existence + allowed values. Returns `(warning, error)`.
- `process_functions.validate_genie_identifier(identifiers, center, filename, col)` — validate GENIE-{CENTER}-... ID format.
- `process_functions.get_clinical_dataframe(filePathList)` — read and merge clinical file pair on PATIENT_ID.
- `process_functions.get_row_indices_for_invalid_column_values(df, col, possible_values)` — returns pd.Index of invalid rows.
- `process_functions.get_message_for_invalid_column_value(df, col, possible_values, filename)` — generates Confluence-format error message with row indices.
- `extract.get_syntabledf(syn, query_string)` — query Synapse table, return DataFrame.
- `validate._validate_chromosome(df, col, fileformat)` — chromosome validation with automatic "chr" prefix removal.
- `validate.check_values_between_two_df(df1, df2, col, ...)` — cross-file ID existence check with auto-uppercasing.
- `validate.get_invalid_allele_rows(df, col, allowed_comb_alleles, ...)` — allele validation with regex.

## DataFrame Conventions

- Use `fillna("")` before string comparisons on columns that may contain NAs — the primary key generation in `load.py` depends on this.
- Use `pd.Int64Dtype()` for nullable integer columns — pandas casts int columns with NAs to float otherwise.
- The `.0` stripping hack in `load.store_database()` does `.replace(".0,", ",").replace(".0\n", "\n")` on CSV output before Synapse upload. Do not introduce column values that legitimately end in `.0`.
- Column order matters: `load._reorder_new_dataset()` forces new data to match existing Synapse table column order. Never assume column order is arbitrary.
- All column names are uppercased on read: `df.columns = [col.upper() for col in df.columns]`.
- `pd.options.mode.chained_assignment = None` is set globally in `process_functions.py`.

## Synapse Interaction Patterns

- `load.store_file()` auto-adds GitHub tag provenance via `executed=f"https://github.com/Sage-Bionetworks/Genie/tree/v{__version__}"`. Do not add separate provenance.
- `load.store_table()` catches `SynapseTimeoutError` silently — intentional because timeout occurs during indexing, not actual failure.
- `load.update_table()` generates a `UNIQUE_KEY` column by space-joining primary key columns. This is the primary key convention for all table comparisons.
- Set `syn.table_query_timeout = 50000` before large table queries (done in `process_mutation.py` but not universally).

## Constraints

- `config.py` discovers subclasses via `importlib.import_module()` + `__subclasses__()`. Subclasses are NOT found unless their module is imported first. If adding a new format package, it must be passed via `--format_registry_packages`.
- The concurrency lock in `bin/input_to_database.py` uses `isProcessing` annotation on the centerMapping entity. If the script crashes, the lock stays True and must be manually reset in Synapse.
