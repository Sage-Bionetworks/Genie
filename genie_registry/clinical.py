"""Clinical file format validation and processing"""
# from __future__ import annotations
import datetime
from io import StringIO
import logging
import os

import pandas as pd
import synapseclient

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions
from genie.database_to_staging import redact_phi

logger = logging.getLogger(__name__)


def _check_year(
    clinicaldf: pd.DataFrame,
    year_col: int,
    filename: str,
    allowed_string_values: list = None,
) -> str:
    """Check year columns

    Args:
        clinicaldf: Clinical dataframe
        year_col: YEAR column
        filename: Name of file
        allowed_string_values: list of other allowed string values

    Returns:
        Error message
    """
    error = ""
    if allowed_string_values is None:
        allowed_string_values = []
    if process_functions.checkColExist(clinicaldf, year_col):
        # Deal with pre-redacted values and other allowed strings
        # first because can't int(text) because there are
        # instances that have <YYYY
        year_series = clinicaldf[year_col][
            ~clinicaldf[year_col].isin(allowed_string_values)
        ]
        year_now = datetime.datetime.utcnow().year
        try:
            years = year_series.apply(
                lambda x: datetime.datetime.strptime(str(int(x)), "%Y").year > year_now
            )
            # Make sure that none of the years are greater than the current
            # year.  It can be the same, but can't future years.
            assert not years.any()
        except Exception:
            error = (
                f"{filename}: Please double check your {year_col} "
                "column, it must be an integer in YYYY format "
                f"<= {year_now}"
            )
            # Tack on allowed string values
            if allowed_string_values:
                error += " or '{}'.\n".format("', '".join(allowed_string_values))
            else:
                error += ".\n"
    else:
        error = f"{filename}: Must have {year_col} column.\n"

    return error


def _check_int_dead_consistency(clinicaldf: pd.DataFrame) -> str:
    """Check if vital status interval and dead column are consistent

    Args:
        clinicaldf: Clinical Data Frame

    Returns:
        Error message if values and inconsistent or blank string
    """
    cols = ["INT_DOD", "DEAD"]
    for col in cols:
        # Return empty string is columns don't exist because this error
        # is already handled.
        if not process_functions.checkColExist(clinicaldf, col):
            return ""
    is_dead = clinicaldf["DEAD"].astype(str) == "True"
    is_alive = clinicaldf["DEAD"].astype(str) == "False"
    allowed_str = [
        "Unknown",
        "Not Collected",
        "Not Applicable",
        "Not Released",
    ]
    is_str = clinicaldf["DEAD"].isin(allowed_str)
    # Check that all string values are equal each other
    is_equal = all(clinicaldf.loc[is_str, "DEAD"] == clinicaldf.loc[is_str, "INT_DOD"])
    # If dead, int column can't be Not Applicable
    # If alive, int column must be Not Applicable
    if (
        any(clinicaldf.loc[is_dead, "INT_DOD"] == "Not Applicable")
        or not all(clinicaldf.loc[is_alive, "INT_DOD"] == "Not Applicable")
        or not is_equal
    ):
        return (
            "Patient Clinical File: DEAD value is inconsistent with INT_DOD "
            "for at least one patient.\n"
        )
    return ""


def _check_int_year_consistency(
    clinicaldf: pd.DataFrame, cols: list, string_vals: list
) -> str:
    """
    Check if vital status interval and year columns are consistent in
    their values

    Args:
        clinicaldf: Clinical Data Frame
        cols: Columns in the clinical data frame
        string_vals: String values that aren't integers

    Returns:
        Error message if values and inconsistent or blank string
    """
    interval_col = ""
    year_col = ""
    for col in cols:
        # This is assuming that interval and year columns start with
        # INT/YEAR
        interval_col = col if col.startswith("INT") else interval_col
        year_col = col if col.startswith("YEAR") else year_col
        # Return empty string is columns don't exist because this error
        # is already handled.
        if not process_functions.checkColExist(clinicaldf, col):
            return ""

    is_text_inconsistent = False
    # Get index of all rows that have 'missing' values
    for str_val in string_vals:
        # n string values per row
        n_str = (clinicaldf[cols] == str_val).sum(axis=1)
        # year can be known with unknown interval value
        # otherwise must be all numeric or the same text value
        if str_val == "Unknown":
            if ((n_str == 1) & (clinicaldf[interval_col] != "Unknown")).any():
                is_text_inconsistent = True
        else:
            if n_str.between(0, len(cols), inclusive="neither").any():
                is_text_inconsistent = True

    is_redaction_inconsistent = False
    # Check that the redacted values are consistent
    is_redacted_int_89 = clinicaldf[interval_col] == ">32485"
    is_redacted_year_89 = clinicaldf[year_col] == ">89"
    is_redacted_int = clinicaldf[interval_col] == "<6570"
    is_redacted_year = clinicaldf[year_col] == "<18"
    if any(is_redacted_int != is_redacted_year) or any(
        is_redacted_int_89 != is_redacted_year_89
    ):
        is_redaction_inconsistent = True

    col_strs = ", ".join(cols)
    if is_text_inconsistent and is_redaction_inconsistent:
        return (
            "Patient: you have inconsistent redaction and text "
            f"values in {col_strs}.\n"
        )
    if is_redaction_inconsistent:
        return f"Patient: you have inconsistent redaction values in {col_strs}.\n"
    if is_text_inconsistent:
        return f"Patient: you have inconsistent text values in {col_strs}.\n"

    return ""


# PROCESSING
def remap_clinical_values(
    clinicaldf: pd.DataFrame,
    sex_mapping: pd.DataFrame,
    race_mapping: pd.DataFrame,
    ethnicity_mapping: pd.DataFrame,
    sampletype_mapping: pd.DataFrame,
) -> pd.DataFrame:
    """Remap clinical attributes from integer to string values

    Args:
        clinicaldf: Clinical data
        sex_mapping: Sex mapping data
        race_mapping: Race mapping data
        ethnicity_mapping: Ethnicity mapping data
        sample_type: Sample type mapping data

    Returns:
        Mapped clinical dataframe
    """

    race_mapping.index = race_mapping["CODE"]
    race_dict = race_mapping.to_dict()

    ethnicity_mapping.index = ethnicity_mapping["CODE"]
    ethnicity_dict = ethnicity_mapping.to_dict()

    sex_mapping.index = sex_mapping["CODE"]
    sex_dict = sex_mapping.to_dict()

    sampletype_mapping.index = sampletype_mapping["CODE"]
    sampletype_dict = sampletype_mapping.to_dict()

    if clinicaldf.get("SAMPLE_TYPE") is not None:
        clinicaldf["SAMPLE_TYPE_DETAILED"] = clinicaldf["SAMPLE_TYPE"]

    # Use pandas mapping feature
    clinicaldf = clinicaldf.replace(
        {
            "PRIMARY_RACE": race_dict["CBIO_LABEL"],
            "SECONDARY_RACE": race_dict["CBIO_LABEL"],
            "TERTIARY_RACE": race_dict["CBIO_LABEL"],
            "SAMPLE_TYPE": sampletype_dict["CBIO_LABEL"],
            "SAMPLE_TYPE_DETAILED": sampletype_dict["DESCRIPTION"],
            "SEX": sex_dict["CBIO_LABEL"],
            "ETHNICITY": ethnicity_dict["CBIO_LABEL"],
        }
    )

    return clinicaldf


class Clinical(FileTypeFormat):

    _fileType = "clinical"

    # _process_kwargs = [
    #     "newPath", "patientSynId", "sampleSynId",
    #     "parentId", "retractedSampleSynId", "retractedPatientSynId"]
    _process_kwargs = [
        "newPath",
        "parentId",
        "clinicalTemplate",
        "sample",
        "patient",
        "patientCols",
        "sampleCols",
    ]

    # VALIDATE FILE NAME
    def _validateFilename(self, filePath):
        if len(filePath) == 1:
            assert os.path.basename(filePath[0]) == "data_clinical_supp_{}.txt".format(
                self.center
            )
        else:
            required = pd.Series(
                [
                    "data_clinical_supp_sample_{}.txt".format(self.center),
                    "data_clinical_supp_patient_{}.txt".format(self.center),
                ]
            )
            assert all(required.isin([os.path.basename(i) for i in filePath]))

    # Update clinical file with the correct mappings
    def update_clinical(self, row):
        """Transform the values of each row of the clinical file"""
        # Must create copy or else it will overwrite the original row
        x = row.copy()

        # BIRTH YEAR
        if x.get("BIRTH_YEAR") is not None:
            # BIRTH YEAR (Check if integer)
            if process_functions.checkInt(x["BIRTH_YEAR"]):
                x["BIRTH_YEAR"] = int(x["BIRTH_YEAR"])

        # AGE AT SEQ REPORT
        if x.get("AGE_AT_SEQ_REPORT") is not None:
            if process_functions.checkInt(x["AGE_AT_SEQ_REPORT"]):
                x["AGE_AT_SEQ_REPORT"] = int(x["AGE_AT_SEQ_REPORT"])

        # SEQ ASSAY ID
        if x.get("SEQ_ASSAY_ID") is not None:
            x["SEQ_ASSAY_ID"] = x["SEQ_ASSAY_ID"].replace("_", "-")
            # standardize all SEQ_ASSAY_ID with uppercase
            x["SEQ_ASSAY_ID"] = x["SEQ_ASSAY_ID"].upper()

        if x.get("SEQ_DATE") is not None:
            x["SEQ_DATE"] = x["SEQ_DATE"].title()
            x["SEQ_YEAR"] = (
                int(str(x["SEQ_DATE"]).split("-")[1])
                if str(x["SEQ_DATE"]) != "Release"
                else float("nan")
            )

        if x.get("YEAR_CONTACT") is not None:
            if process_functions.checkInt(x["YEAR_CONTACT"]):
                x["YEAR_CONTACT"] = int(x["YEAR_CONTACT"])

        if x.get("YEAR_DEATH") is not None:
            if process_functions.checkInt(x["YEAR_DEATH"]):
                x["YEAR_DEATH"] = int(x["YEAR_DEATH"])

        # TRIM EVERY COLUMN MAKE ALL DASHES
        for i in x.keys():
            if isinstance(x[i], str):
                x[i] = x[i].strip(" ")
        return x

    def uploadMissingData(
        self, df: pd.DataFrame, col: str, dbSynId: str, stagingSynId: str
    ):
        """Uploads missing clinical samples / patients

        Args:
            df (pd.DataFrame): dataframe with clinical data
            col (str): column in dataframe. Usually SAMPLE_ID or PATIENT_ID.
            dbSynId (str): Synapse table Synapse id
            stagingSynId (str): Center Synapse staging Id
        """
        path = os.path.join(
            process_functions.SCRIPT_DIR, f"{self._fileType}_missing_{col}.csv"
        )
        # PLFM-7428 - there are limits on a "not in" function on Synapse tables
        center_samples = self.syn.tableQuery(
            f"select {col} from {dbSynId} where " f"CENTER='{self.center}'"
        )
        center_samples_df = center_samples.asDataFrame()
        # Get all the samples that are in the database but missing from
        # the input file
        missing_df = center_samples_df[col][~center_samples_df[col].isin(df[col])]
        missing_df.to_csv(path, index=False)
        self.syn.store(synapseclient.File(path, parent=stagingSynId))
        os.remove(path)

    def _process(self, clinical, clinicalTemplate):
        # Capitalize all clinical dataframe columns
        clinical.columns = [col.upper() for col in clinical.columns]
        clinical = clinical.fillna("")
        # clinicalMerged = clinical.merge(clinicalTemplate,how='outer')
        # Remove unwanted clinical columns prior to update
        # clinicalMerged = clinicalMerged.drop(clinicalMerged.columns[
        #    ~clinicalMerged.columns.isin(clinicalTemplate.columns)],1)
        ethnicity_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['ethnicity_mapping']}"
        )
        race_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['race_mapping']}"
        )
        sex_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['sex_mapping']}"
        )
        sampletype_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['sampletype_mapping']}"
        )
        # Attach MSK to centers
        # clinicalMerged = clinicalMerged.fillna("")
        clinical = remap_clinical_values(
            clinical,
            sex_mapping,
            race_mapping,
            ethnicity_mapping,
            sampletype_mapping,
        )
        remapped_clindf = clinical.apply(self.update_clinical, axis=1)
        # Some columns may have been added during update,
        # remove unwanted columns again
        keep_cols_idx = remapped_clindf.columns.isin(clinicalTemplate.columns)
        remapped_clindf = remapped_clindf.drop(
            columns=remapped_clindf.columns[~keep_cols_idx]
        )
        remapped_clindf["CENTER"] = self.center
        return remapped_clindf

    def preprocess(self, newpath):
        """
        Gather preprocess parameters

        Args:
            filePath: Path to file

        Returns:
            dict with keys - 'clinicalTemplate', 'sample', 'patient',
                             'patientCols', 'sampleCols'
        """
        # These synapse ids for the clinical tier release scope is
        # hardcoded because it never changes
        # TODO: Add clinical tier release scope to GENIE config
        patient_cols_table = self.syn.tableQuery(
            "select fieldName from syn8545211 where "
            "patient is True and inClinicalDb is True"
        )
        patient_cols = patient_cols_table.asDataFrame()["fieldName"].tolist()
        sample_cols_table = self.syn.tableQuery(
            "select fieldName from syn8545211 where "
            "sample is True and inClinicalDb is True"
        )
        sample_cols = sample_cols_table.asDataFrame()["fieldName"].tolist()
        clinicalTemplate = pd.DataFrame(columns=set(patient_cols + sample_cols))
        sample = True
        patient = True

        return {
            "clinicalTemplate": clinicalTemplate,
            "sample": sample,
            "patient": patient,
            "patientCols": patient_cols,
            "sampleCols": sample_cols,
        }

    def process_steps(
        self,
        clinicalDf,
        newPath,
        parentId,
        clinicalTemplate,
        sample,
        patient,
        patientCols,
        sampleCols,
    ):
        """Process clincial file, redact PHI values, upload to clinical
        database
        """
        patient_synid = self.genie_config["patient"]
        sample_synid = self.genie_config["sample"]

        newClinicalDf = self._process(clinicalDf, clinicalTemplate)
        newClinicalDf = redact_phi(newClinicalDf)

        if patient:
            cols = newClinicalDf.columns[newClinicalDf.columns.isin(patientCols)]
            patientClinical = newClinicalDf[cols].drop_duplicates("PATIENT_ID")
            self.uploadMissingData(
                patientClinical, "PATIENT_ID", patient_synid, parentId
            )

            process_functions.updateData(
                self.syn,
                patient_synid,
                patientClinical,
                self.center,
                col=cols.tolist(),
                toDelete=True,
            )
        if sample:
            cols = newClinicalDf.columns[newClinicalDf.columns.isin(sampleCols)]
            if sum(newClinicalDf["SAMPLE_ID"].duplicated()) > 0:
                logger.error(
                    "There are duplicated samples, " "and the duplicates are removed"
                )
            sampleClinical = newClinicalDf[cols].drop_duplicates("SAMPLE_ID")
            # Exclude all clinical samples with wrong oncotree codes
            oncotree_mapping = pd.DataFrame()
            oncotree_mapping_dict = process_functions.get_oncotree_code_mappings(
                self.genie_config["oncotreeLink"]
            )
            # Add in unknown key for oncotree code
            oncotree_mapping_dict["UNKNOWN"] = {}
            oncotree_mapping["ONCOTREE_CODE"] = list(oncotree_mapping_dict.keys())
            # Make oncotree codes uppercase (SpCC/SPCC)
            sampleClinical["ONCOTREE_CODE"] = (
                sampleClinical["ONCOTREE_CODE"].astype(str).str.upper()
            )
            sampleClinical = sampleClinical[
                sampleClinical["ONCOTREE_CODE"].isin(oncotree_mapping["ONCOTREE_CODE"])
            ]
            self.uploadMissingData(sampleClinical, "SAMPLE_ID", sample_synid, parentId)
            process_functions.updateData(
                self.syn,
                sample_synid,
                sampleClinical,
                self.center,
                col=cols.tolist(),
                toDelete=True,
            )

        newClinicalDf.to_csv(newPath, sep="\t", index=False)
        return newPath

    # VALIDATION
    def _validate(self, clinicaldf):
        """
        This function validates the clinical file to make sure it adhere
        to the clinical SOP.

        Args:
            clinicalDF: Merged clinical file with patient and sample
                        information

        Returns:
            Error message
        """
        total_error = StringIO()
        warning = StringIO()

        clinicaldf.columns = [col.upper() for col in clinicaldf.columns]
        # CHECK: for empty rows
        empty_rows = clinicaldf.isnull().values.all(axis=1)
        if empty_rows.any():
            total_error.write("Clinical file(s): No empty rows allowed.\n")
            # Remove completely empty rows to speed up processing
            clinicaldf = clinicaldf[~empty_rows]

        clinicaldf = clinicaldf.fillna("")

        oncotree_mapping_dict = process_functions.get_oncotree_code_mappings(
            self.genie_config["oncotreeLink"]
        )
        oncotree_mapping = pd.DataFrame(
            {"ONCOTREE_CODE": list(oncotree_mapping_dict.keys())}
        )

        ethnicity_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['ethnicity_mapping']}"
        )
        race_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['race_mapping']}"
        )
        sex_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['sex_mapping']}"
        )
        sampletype_mapping = process_functions.get_syntabledf(
            self.syn, f"select * from {self.genie_config['sampletype_mapping']}"
        )
        # CHECK: SAMPLE_ID
        sample_id = "SAMPLE_ID"
        haveSampleColumn = process_functions.checkColExist(clinicaldf, sample_id)

        if not haveSampleColumn:
            total_error.write("Sample Clinical File: Must have SAMPLE_ID column.\n")
        else:
            if sum(clinicaldf[sample_id].duplicated()) > 0:
                total_error.write(
                    "Sample Clinical File: No duplicated SAMPLE_ID "
                    "allowed.\nIf there are no duplicated "
                    "SAMPLE_IDs, and both sample and patient files are "
                    "uploaded, then please check to make sure no duplicated "
                    "PATIENT_IDs exist in the patient clinical file.\n"
                )
            error = process_functions.validate_genie_identifier(
                identifiers=clinicaldf[sample_id],
                center=self.center,
                filename="Sample Clinical File",
                col="SAMPLE_ID",
            )
            total_error.write(error)

        # CHECK: PATIENT_ID
        patientId = "PATIENT_ID"
        # #CHECK: PATIENT_ID IN SAMPLE FILE
        havePatientColumn = process_functions.checkColExist(clinicaldf, patientId)

        if not havePatientColumn:
            total_error.write("Patient Clinical File: Must have PATIENT_ID column.\n")
        else:
            if not all(clinicaldf[patientId].str.startswith(f"GENIE-{self.center}")):
                total_error.write(
                    "Patient Clinical File: "
                    f"PATIENT_ID must start with GENIE-{self.center}\n"
                )
            if any(clinicaldf[patientId].str.len() >= 50):
                total_error.write(
                    "Patient Clinical File: PATIENT_ID must have less than "
                    "50 characters.\n"
                )
        # CHECK: within the sample file that the sample ids match
        # the patient ids
        if haveSampleColumn and havePatientColumn:
            # Make sure sample and patient ids are string cols
            clinicaldf[sample_id] = clinicaldf[sample_id].astype(str)
            clinicaldf[patientId] = clinicaldf[patientId].astype(str)
            if not all(
                [
                    patient in sample
                    for sample, patient in zip(
                        clinicaldf[sample_id], clinicaldf[patientId]
                    )
                ]
            ):

                total_error.write(
                    "Sample Clinical File: PATIENT_ID's much be contained in "
                    "the SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n"
                )
            # #CHECK: All samples must have associated patient data
            # (GENIE requires patient data)
            if not all(clinicaldf[patientId] != ""):
                total_error.write(
                    "Patient Clinical File: All samples must have associated "
                    "patient information and no null patient ids allowed. "
                    "These samples are missing patient data: {}\n".format(
                        ", ".join(
                            clinicaldf[sample_id][clinicaldf[patientId] == ""].unique()
                        )
                    )
                )

            # CHECK: All patients should have associated sample data
            if not all(clinicaldf[sample_id] != ""):
                # ## MAKE WARNING FOR NOW###
                warning.write(
                    "Sample Clinical File: All patients must have associated "
                    "sample information. These patients are missing sample "
                    "data: {}\n".format(
                        ", ".join(
                            clinicaldf[patientId][clinicaldf[sample_id] == ""].unique()
                        )
                    )
                )

        # CHECK: AGE_AT_SEQ_REPORT
        age = "AGE_AT_SEQ_REPORT"
        haveColumn = process_functions.checkColExist(clinicaldf, age)
        if haveColumn:
            # Deal with HIPAA converted rows from DFCI
            # First for loop can't int(text) because there
            # are instances that have <3435
            age_seq_report_df = clinicaldf[
                ~clinicaldf[age].isin(["Unknown", ">32485", "<6570"])
            ]

            # age_seq_report_df[age] = \
            #     remove_greaterthan_lessthan_str(age_seq_report_df[age])

            if not all([process_functions.checkInt(i) for i in age_seq_report_df[age]]):
                total_error.write(
                    "Sample Clinical File: Please double check your "
                    "AGE_AT_SEQ_REPORT. It must be an integer, 'Unknown', "
                    "'>32485', '<6570'.\n"
                )
            else:
                age_seq_report_df[age] = age_seq_report_df[age].astype(int)
                median_age = age_seq_report_df[age].median()
                if median_age < 100:
                    total_error.write(
                        "Sample Clinical File: Please double check your "
                        "AGE_AT_SEQ_REPORT. You may be reporting this value "
                        "in YEARS, please report in DAYS.\n"
                    )
        else:
            total_error.write(
                "Sample Clinical File: Must have AGE_AT_SEQ_REPORT column.\n"
            )

        # CHECK: ONCOTREE_CODE
        haveColumn = process_functions.checkColExist(clinicaldf, "ONCOTREE_CODE")
        maleOncoCodes = ["TESTIS", "PROSTATE", "PENIS"]
        womenOncoCodes = ["CERVIX", "VULVA", "UTERUS", "OVARY"]
        if haveColumn:
            # Make oncotree codes uppercase (SpCC/SPCC)
            clinicaldf["ONCOTREE_CODE"] = (
                clinicaldf["ONCOTREE_CODE"].astype(str).str.upper()
            )

            oncotree_codes = clinicaldf["ONCOTREE_CODE"][
                clinicaldf["ONCOTREE_CODE"] != "UNKNOWN"
            ]

            if not all(oncotree_codes.isin(oncotree_mapping["ONCOTREE_CODE"])):
                unmapped_oncotrees = oncotree_codes[
                    ~oncotree_codes.isin(oncotree_mapping["ONCOTREE_CODE"])
                ]
                total_error.write(
                    "Sample Clinical File: Please double check that all your "
                    "ONCOTREE CODES exist in the mapping. You have {} samples "
                    "that don't map. These are the codes that "
                    "don't map: {}\n".format(
                        len(unmapped_oncotrees),
                        ",".join(set(unmapped_oncotrees)),
                    )
                )
            # Should add the SEX mismatch into the dashboard file
            if (
                process_functions.checkColExist(clinicaldf, "SEX")
                and "oncotree_mapping_dict" in locals()
                and havePatientColumn
                and haveSampleColumn
            ):

                wrongCodeSamples = []
                # This is to check if oncotree codes match the sex,
                # returns list of samples that have conflicting codes and sex
                for code, patient, sample in zip(
                    clinicaldf["ONCOTREE_CODE"],
                    clinicaldf["PATIENT_ID"],
                    clinicaldf["SAMPLE_ID"],
                ):

                    if (
                        oncotree_mapping_dict.get(code) is not None
                        and sum(clinicaldf["PATIENT_ID"] == patient) > 0
                    ):

                        primaryCode = oncotree_mapping_dict[code][
                            "ONCOTREE_PRIMARY_NODE"
                        ]

                        sex = clinicaldf["SEX"][
                            clinicaldf["PATIENT_ID"] == patient
                        ].values[0]
                        sex = float("nan") if sex == "" else float(sex)
                        if (
                            oncotree_mapping_dict[code]["ONCOTREE_PRIMARY_NODE"]
                            in maleOncoCodes
                            and sex != 1.0
                        ):

                            wrongCodeSamples.append(sample)
                        if (
                            oncotree_mapping_dict[code]["ONCOTREE_PRIMARY_NODE"]
                            in womenOncoCodes
                            and sex != 2.0
                        ):

                            wrongCodeSamples.append(sample)
                if len(wrongCodeSamples) > 0:
                    warning.write(
                        "Sample Clinical File: Some SAMPLE_IDs have "
                        "conflicting SEX and ONCOTREE_CODES: {}\n".format(
                            ",".join(wrongCodeSamples)
                        )
                    )
        else:
            total_error.write("Sample Clinical File: Must have ONCOTREE_CODE column.\n")

        warn, error = process_functions.check_col_and_values(
            clinicaldf,
            "SAMPLE_TYPE",
            sampletype_mapping["CODE"].tolist(),
            "Sample Clinical File",
            required=True,
        )
        total_error.write(error)

        # CHECK: SEQ_ASSAY_ID
        haveColumn = process_functions.checkColExist(clinicaldf, "SEQ_ASSAY_ID")
        if haveColumn:
            if not all([i != "" for i in clinicaldf["SEQ_ASSAY_ID"]]):
                total_error.write(
                    "Sample Clinical File: Please double check your "
                    "SEQ_ASSAY_ID columns, there are empty rows.\n"
                )
            # must remove empty seq assay ids first
            # Checking if seq assay ids start with the center name
            empty_seq_idx = clinicaldf.SEQ_ASSAY_ID != ""
            seqassay_ids = clinicaldf.SEQ_ASSAY_ID[empty_seq_idx]
            uniq_seqassay_ids = seqassay_ids.unique()
            invalid_seqassay = []
            for seqassay in uniq_seqassay_ids:
                # SEQ Ids are all capitalized now, so no need to check
                # for differences in case
                if not seqassay.upper().startswith(self.center):
                    invalid_seqassay.append(seqassay)
            if invalid_seqassay:
                total_error.write(
                    "Sample Clinical File: Please make sure your "
                    "SEQ_ASSAY_IDs start with your center "
                    "abbreviation: {}.\n".format(", ".join(invalid_seqassay))
                )
        else:
            total_error.write("Sample Clinical File: Must have SEQ_ASSAY_ID column.\n")

        haveColumn = process_functions.checkColExist(clinicaldf, "SEQ_DATE")
        seq_date_error = (
            "Sample Clinical File: SEQ_DATE must be one of five values- "
            "For Jan-March: use Jan-YEAR. "
            "For Apr-June: use Apr-YEAR. "
            "For July-Sep: use Jul-YEAR. "
            "For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) "
            "For values that don't have SEQ_DATES that "
            "you want released use 'release'.\n"
        )

        if haveColumn:
            clinicaldf["SEQ_DATE"] = [
                i.title() for i in clinicaldf["SEQ_DATE"].astype(str)
            ]

            seqdate = clinicaldf["SEQ_DATE"][clinicaldf["SEQ_DATE"] != "Release"]
            if sum(clinicaldf["SEQ_DATE"] == "") > 0:
                total_error.write(
                    "Sample Clinical File: Samples without SEQ_DATEs will "
                    "NOT be released.\n"
                )
            try:
                if not seqdate.empty:
                    seqdate.apply(
                        lambda date: datetime.datetime.strptime(date, "%b-%Y")
                    )
                    if not seqdate.str.startswith(("Jan", "Apr", "Jul", "Oct")).all():
                        total_error.write(seq_date_error)
            except ValueError:
                total_error.write(seq_date_error)
        else:
            total_error.write("Sample Clinical File: Must have SEQ_DATE column.\n")

        # CHECK: BIRTH_YEAR
        error = _check_year(
            clinicaldf=clinicaldf,
            year_col="BIRTH_YEAR",
            filename="Patient Clinical File",
            allowed_string_values=["Unknown", ">89", "<18"],
        )
        total_error.write(error)

        # CHECK: YEAR DEATH
        error = _check_year(
            clinicaldf=clinicaldf,
            year_col="YEAR_DEATH",
            filename="Patient Clinical File",
            allowed_string_values=[
                "Unknown",
                "Not Collected",
                "Not Applicable",
                "Not Released",
                ">89",
                "<18",
            ],
        )
        total_error.write(error)

        # CHECK: YEAR CONTACT
        error = _check_year(
            clinicaldf=clinicaldf,
            year_col="YEAR_CONTACT",
            filename="Patient Clinical File",
            allowed_string_values=[
                "Unknown",
                "Not Collected",
                "Not Released",
                ">89",
                "<18",
            ],
        )
        total_error.write(error)

        # CHECK: INT CONTACT
        haveColumn = process_functions.checkColExist(clinicaldf, "INT_CONTACT")
        if haveColumn:
            if not all(
                [
                    process_functions.checkInt(i)
                    for i in clinicaldf.INT_CONTACT
                    if i
                    not in [
                        ">32485",
                        "<6570",
                        "Unknown",
                        "Not Collected",
                        "Not Released",
                    ]
                ]
            ):

                total_error.write(
                    "Patient Clinical File: Please double check your "
                    "INT_CONTACT column, it must be an integer, '>32485', "
                    "'<6570', 'Unknown', 'Not Released' or 'Not Collected'.\n"
                )
        else:
            total_error.write("Patient Clinical File: Must have INT_CONTACT column.\n")

        # INT DOD
        haveColumn = process_functions.checkColExist(clinicaldf, "INT_DOD")
        if haveColumn:
            if not all(
                [
                    process_functions.checkInt(i)
                    for i in clinicaldf.INT_DOD
                    if i
                    not in [
                        ">32485",
                        "<6570",
                        "Unknown",
                        "Not Collected",
                        "Not Applicable",
                        "Not Released",
                    ]
                ]
            ):

                total_error.write(
                    "Patient Clinical File: Please double check your INT_DOD "
                    "column, it must be an integer, '>32485', '<6570', "
                    "'Unknown', 'Not Collected', 'Not Released' or "
                    "'Not Applicable'.\n"
                )
        else:
            total_error.write("Patient Clinical File: Must have INT_DOD column.\n")

        haveColumn = process_functions.checkColExist(clinicaldf, "DEAD")
        if haveColumn:
            # Need to have check_bool function
            if not all(
                [
                    str(i).upper() in ["TRUE", "FALSE"]
                    for i in clinicaldf.DEAD
                    if i not in ["Unknown", "Not Collected", "Not Released"]
                ]
            ):
                total_error.write(
                    "Patient Clinical File: Please double check your "
                    "DEAD column, it must be True, False, 'Unknown', "
                    "'Not Released' or 'Not Collected'.\n"
                )
        else:
            total_error.write("Patient Clinical File: Must have DEAD column.\n")
        # CHECK: contact vital status value consistency
        contact_error = _check_int_year_consistency(
            clinicaldf=clinicaldf,
            cols=["YEAR_CONTACT", "INT_CONTACT"],
            string_vals=["Not Collected", "Unknown", "Not Released"],
        )
        total_error.write(contact_error)

        # CHECK: death vital status value consistency
        death_error = _check_int_year_consistency(
            clinicaldf=clinicaldf,
            cols=["YEAR_DEATH", "INT_DOD"],
            string_vals=[
                "Not Collected",
                "Unknown",
                "Not Applicable",
                "Not Released",
            ],
        )
        total_error.write(death_error)
        death_error = _check_int_dead_consistency(clinicaldf=clinicaldf)
        total_error.write(death_error)

        # CHECK: SAMPLE_CLASS is optional attribute
        have_column = process_functions.checkColExist(clinicaldf, "SAMPLE_CLASS")
        if have_column:
            sample_class_vals = pd.Series(clinicaldf["SAMPLE_CLASS"].unique().tolist())
            if not sample_class_vals.isin(["Tumor", "cfDNA"]).all():
                total_error.write(
                    "Sample Clinical File: SAMPLE_CLASS column must "
                    "be 'Tumor', or 'cfDNA'\n"
                )

        # CHECK: PRIMARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicaldf,
            "PRIMARY_RACE",
            race_mapping["CODE"].tolist(),
            "Patient Clinical File",
        )
        warning.write(warn)
        total_error.write(error)

        # CHECK: SECONDARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicaldf,
            "SECONDARY_RACE",
            race_mapping["CODE"].tolist(),
            "Patient Clinical File",
        )
        warning.write(warn)
        total_error.write(error)

        # CHECK: TERTIARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicaldf,
            "TERTIARY_RACE",
            race_mapping["CODE"].tolist(),
            "Patient Clinical File",
        )
        warning.write(warn)
        total_error.write(error)

        # CHECK: SEX
        warn, error = process_functions.check_col_and_values(
            clinicaldf,
            "SEX",
            sex_mapping["CODE"].tolist(),
            "Patient Clinical File",
            required=True,
        )
        warning.write(warn)
        total_error.write(error)

        # CHECK: ETHNICITY
        warn, error = process_functions.check_col_and_values(
            clinicaldf,
            "ETHNICITY",
            ethnicity_mapping["CODE"].tolist(),
            "Patient Clinical File",
        )
        warning.write(warn)
        total_error.write(error)

        return total_error.getvalue(), warning.getvalue()

    def _get_dataframe(self, filePathList):
        clinicaldf = pd.read_csv(filePathList[0], sep="\t", comment="#")
        clinicaldf.columns = [col.upper() for col in clinicaldf.columns]

        if len(filePathList) > 1:
            other_clinicaldf = pd.read_csv(filePathList[1], sep="\t", comment="#")
            other_clinicaldf.columns = [col.upper() for col in other_clinicaldf.columns]

            try:
                clinicaldf = clinicaldf.merge(other_clinicaldf, on="PATIENT_ID")
            except Exception:
                raise ValueError(
                    (
                        "If submitting separate patient and sample files, "
                        "they both must have the PATIENT_ID column"
                    )
                )
            # Must figure out which is sample and which is patient
            if "sample" in filePathList[0]:
                sample = clinicaldf
                patient = other_clinicaldf
            else:
                sample = other_clinicaldf
                patient = clinicaldf

            if not all(sample["PATIENT_ID"].isin(patient["PATIENT_ID"])):
                raise ValueError(
                    (
                        "Patient Clinical File: All samples must have associated "
                        "patient information"
                    )
                )

        return clinicaldf
