<!-- Last reviewed: 2026-03 -->

## Project

Three CLI scripts mapping to the three GENIE pipeline stages: `input_to_database.py` (validation + ingestion), `database_to_staging.py` (consortium release), `consortium_to_public.py` (public release). Each wraps core `genie/` functions with argparse.

## Argument Naming

CLI arguments use camelCase (`--onlyValidate`, `--deleteOld`, `--createNewMafDatabase`). argparse stores them as-is on the namespace, but `main()` functions receive snake_case parameter names. Match existing naming conventions when adding new arguments.

## Concurrency Lock

`input_to_database.py` uses an `isProcessing` annotation on the `centerMapping` Synapse entity as a mutex lock. Flow: check if True → raise if locked → set True → process → set False. If the script crashes mid-run, the lock stays True and must be manually reset in Synapse. This prevents parallel pipeline runs from corrupting data.

## Hardcoded Database Mapping IDs

Three dbMapping Synapse IDs are duplicated across `database_to_staging.py`, `consortium_to_public.py`, all R scripts, and `dashboard_table_updater.py`:
- Test: `syn11600968`
- Staging: `syn12094210`
- Production: `syn10967259`

When updating these IDs, grep the entire repo (Python AND R files). There is no single constant.

## Environment Routing

`database_to_staging.py` and `consortium_to_public.py` use mutually exclusive `--test` and `--staging` flags to select the dbMapping ID. Production is the default (neither flag). `input_to_database.py` uses `--project_id` instead.

## cBioPortal Validator Hack

Both `database_to_staging.py` and `consortium_to_public.py` append `"; exit 0"` to the cBioPortal validator shell command. This forces a success exit code so `subprocess.check_output()` captures output even when the validator reports errors. Do not remove this without changing the subprocess handling.

## Process Tracking

`load.update_process_trackingdf()` is called only in production — never in test or staging mode. Condition: `if not args.test and not args.staging`. Do not add process tracking calls without this guard.

## Constraints

- R scripts are invoked via `subprocess` calls to `Rscript` (e.g., dashboard generation). Templates in `/templates/` are parameterized via Python string replacement BEFORE being passed to R.
- The data guide template (`data_guide_template.Rnw`) uses LaTeX/Sweave. Underscore characters in template values must be escaped as `\_` for LaTeX.
