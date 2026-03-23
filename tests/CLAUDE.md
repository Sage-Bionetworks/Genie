<!-- Last reviewed: 2026-03 -->

## Project

pytest test suite for `genie` and `genie_registry` packages. Tests mirror source modules: `test_clinical.py` tests `genie_registry/clinical.py`, etc.

## Fixtures (conftest.py)

- `syn` — `mock.create_autospec(synapseclient.Synapse)`. Session-scoped. Shared across all tests. Do not modify its default state in tests — add `side_effect` or `return_value` per test instead.
- `genie_config` — dict mapping table type names to Synapse IDs, with `center_config` containing SAGE/TEST/GOLD centers. Session-scoped.
- Format class fixtures: instantiate as `FormatClass(syn, "SAGE")` for validation-only, or `FormatClass(syn, "SAGE", genie_config)` when config is needed. Some accept `ancillary_files` as fourth arg for cross-validation.

## Mock Table Pattern

`createMockTable(dataframe)` is copy-pasted across test files (test_clinical.py, test_cna.py, test_bed.py, etc.) — it is NOT in conftest.py. Pattern:
```python
def createMockTable(dataframe):
    table = mock.create_autospec(synapseclient.table.CsvFileTable)
    table.asDataFrame.return_value = dataframe
    return table
```

Use with query dispatch:
```python
table_query_results_map = {
    ("select * from syn123",): createMockTable(some_df),
    ("select * from syn456",): createMockTable(other_df),
}
syn.tableQuery.side_effect = lambda *args: table_query_results_map[args]
```

## Test Data

- Create inline using `pd.DataFrame(dict(COL1=[...], COL2=[...]))` with **UPPERCASE** column names — because the pipeline uppercases on read.
- No external test data fixture files. All test data lives in the test module.
- Module-level DataFrames are fine for read-only test data shared across functions.
- Synapse entities: `synapseclient.File(name="...", path="...", parentId="syn123")` with attributes set after creation.

## Naming Conventions

- Format class validation: `test_<adjective>_<method>` (e.g., `test_perfect__validate`, `test_missingcols__validate`). Double underscore before private method names.
- Behavior-driven: `test_that_<behavior>` (e.g., `test_that_to_unix_epoch_time_utc_gives_expected_time`).
- Always use `ids=` parameter with `@pytest.mark.parametrize` for clear test output.
- Include ALL pass and fail cases in the same parametrize block.

## Mock Pattern

- Prefer `with patch.object(module, "function_name") as mock_fn:` — place assertions INSIDE the `with` block.
- For multiple mocks: chain in single `with` statement: `with patch.object(...) as m1, patch.object(...) as m2:`.
- For Synapse queries, use `syn.tableQuery.side_effect` pointing to a dispatch function, not `return_value`.
- Verify calls with `assert_called_once_with()` or `assert_has_calls()`.
