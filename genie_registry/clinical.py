"""Clinical file format validation and processing"""
import datetime
import os
import logging
import subprocess
import yaml

import pandas as pd
import synapseclient

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions
from genie.database_to_staging import redact_phi

logger = logging.getLogger(__name__)


# def remove_greaterthan_lessthan_str(col):
#     '''
#     In clinical file, there are redacted value such as >89 and <17.
#     These < and > signs must be removed
#     '''
#     try:
#         col = [
#             text.replace(">", "")
#             if isinstance(text, str) else text for text in col]
#         col = [
#             int(text.replace("<", ""))
#             if isinstance(text, str) and text != "" else text
#             for text in col]
#     except ValueError:
#         pass
#     return(col)
def _check_year(clinicaldf: pd.DataFrame, year_col: int, filename: str,
                allowed_string_values: list = []) -> str:
    """Check year columns

    Args:
        clinicaldf: Clinical dataframe
        year_col: YEAR column
        filename: Name of file
        allowed_string_values: list of other allowed string values

    Returns:
        Error message
    """
    error = ''
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
                lambda x: datetime.datetime.strptime(
                    str(int(x)), '%Y').year > year_now
            )
            # Make sure that none of the years are greater than the current
            # year.  It can be the same, but can't future years.
            assert not years.any()
        except Exception:
            error = (f"{filename}: Please double check your {year_col} "
                     "column, it must be an integer in YYYY format "
                     f"<= {year_now}")
            # Tack on allowed string values
            if allowed_string_values:
                error += " or '{}'.\n".format(
                    "', '".join(allowed_string_values)
                )
            else:
                error += ".\n"
    else:
        error = f"{filename}: Must have {year_col} column.\n"

    return error


# PROCESSING
def remap_clinical_values(clinicaldf: pd.DataFrame, sex_mapping: pd.DataFrame,
                          race_mapping: pd.DataFrame,
                          ethnicity_mapping: pd.DataFrame,
                          sampletype_mapping: pd.DataFrame) -> pd.DataFrame:
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

    race_mapping.index = race_mapping['CODE']
    race_dict = race_mapping.to_dict()

    ethnicity_mapping.index = ethnicity_mapping['CODE']
    ethnicity_dict = ethnicity_mapping.to_dict()

    sex_mapping.index = sex_mapping['CODE']
    sex_dict = sex_mapping.to_dict()

    sampletype_mapping.index = sampletype_mapping['CODE']
    sampletype_dict = sampletype_mapping.to_dict()

    if clinicaldf.get("SAMPLE_TYPE") is not None:
        clinicaldf['SAMPLE_TYPE_DETAILED'] = clinicaldf['SAMPLE_TYPE']

    # Use pandas mapping feature
    clinicaldf = clinicaldf.replace({
        "PRIMARY_RACE": race_dict['CBIO_LABEL'],
        "SECONDARY_RACE": race_dict['CBIO_LABEL'],
        "TERTIARY_RACE": race_dict['CBIO_LABEL'],
        "SAMPLE_TYPE": sampletype_dict['CBIO_LABEL'],
        "SAMPLE_TYPE_DETAILED": sampletype_dict['DESCRIPTION'],
        "SEX": sex_dict['CBIO_LABEL'],
        'ETHNICITY': ethnicity_dict['CBIO_LABEL']
    })

    return clinicaldf


class clinical(FileTypeFormat):

    _fileType = "clinical"

    # _process_kwargs = [
    #     "newPath", "patientSynId", "sampleSynId",
    #     "parentId", "retractedSampleSynId", "retractedPatientSynId"]
    _process_kwargs = [
        "newPath", "parentId", "databaseToSynIdMappingDf", "oncotree_link",
        'clinicalTemplate', 'sample', 'patient', 'patientCols', 'sampleCols']

    _validation_kwargs = ["oncotree_link"]

    # VALIDATE FILE NAME
    def _validateFilename(self, filePath):
        if len(filePath) == 1:
            assert os.path.basename(filePath[0]) == \
                "data_clinical_supp_{}.txt".format(self.center)
        else:
            required = pd.Series([
                "data_clinical_supp_sample_{}.txt".format(self.center),
                "data_clinical_supp_patient_{}.txt".format(self.center)])
            assert all(required.isin([os.path.basename(i) for i in filePath]))

    # Update clinical file with the correct mappings
    def update_clinical(self, row):
        """Transform the values of each row of the clinical file"""
        # Must create copy or else it will overwrite the original row
        x = row.copy()
        # # PATIENT ID
        if x.get("PATIENT_ID") is not None:
            x['PATIENT_ID'] = process_functions.checkGenieId(
                x['PATIENT_ID'], self.center)
        # # RACE
        if x.get('PRIMARY_RACE') is None:
            x['PRIMARY_RACE'] = "Not Collected"

        if x.get('SECONDARY_RACE') is None:
            x['SECONDARY_RACE'] = "Not Collected"

        if x.get('TERTIARY_RACE') is None:
            x['TERTIARY_RACE'] = "Not Collected"
        # ETHNICITY
        if x.get('ETHNICITY') is None:
            x['ETHNICITY'] = "Not Collected"
        # BIRTH YEAR
        if x.get("BIRTH_YEAR") is not None:
            # BIRTH YEAR (Check if integer)
            if process_functions.checkInt(x['BIRTH_YEAR']):
                x['BIRTH_YEAR'] = int(x['BIRTH_YEAR'])

        # SAMPLE ID
        if x.get('SAMPLE_ID') is not None:
            x['SAMPLE_ID'] = process_functions.checkGenieId(
                x['SAMPLE_ID'], self.center)

        # AGE AT SEQ REPORT
        if x.get('AGE_AT_SEQ_REPORT') is not None:
            if process_functions.checkInt(x['AGE_AT_SEQ_REPORT']):
                x['AGE_AT_SEQ_REPORT'] = int(x['AGE_AT_SEQ_REPORT'])

        # SEQ ASSAY ID
        if x.get('SEQ_ASSAY_ID') is not None:
            x['SEQ_ASSAY_ID'] = x['SEQ_ASSAY_ID'].replace('_', '-')
            # standardize all SEQ_ASSAY_ID with uppercase
            x['SEQ_ASSAY_ID'] = x['SEQ_ASSAY_ID'].upper()

        if x.get('SEQ_DATE') is not None:
            x['SEQ_DATE'] = x['SEQ_DATE'].title()
            x['SEQ_YEAR'] = \
                int(str(x['SEQ_DATE']).split("-")[1]) \
                if str(x['SEQ_DATE']) != "Release" else float('nan')

        if x.get('YEAR_CONTACT') is None:
            x['YEAR_CONTACT'] = 'Not Collected'
        else:
            if process_functions.checkInt(x['YEAR_CONTACT']):
                x['YEAR_CONTACT'] = int(x['YEAR_CONTACT'])

        if x.get('YEAR_DEATH') is None:
            x['YEAR_DEATH'] = 'Not Collected'
        else:
            if process_functions.checkInt(x['YEAR_DEATH']):
                x['YEAR_DEATH'] = int(x['YEAR_DEATH'])

        if x.get('INT_CONTACT') is None:
            x['INT_CONTACT'] = 'Not Collected'

        if x.get('INT_DOD') is None:
            x['INT_DOD'] = 'Not Collected'

        if x.get('DEAD') is None:
            x['DEAD'] = 'Not Collected'

        # TRIM EVERY COLUMN MAKE ALL DASHES
        for i in x.keys():
            if isinstance(x[i], str):
                x[i] = x[i].strip(" ")
        return x

    def uploadMissingData(self, df, col, dbSynId, stagingSynId,
                          retractionSynId=None):
        """Uploads missing clinical samples / patients"""
        samples = "','".join(df[col])
        path = os.path.join(
            process_functions.SCRIPT_DIR,
            f"{self._fileType}_missing_{col}.csv"
        )
        missing = self.syn.tableQuery(
            f"select {col} from {dbSynId} where "
            f"CENTER='{self.center}' and {col} not in ('{samples}')"
        )
        missing.asDataFrame().to_csv(path, index=False)
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
        ethnicity_mapping = process_functions.getGenieMapping(
            self.syn, "syn7434242")
        race_mapping = process_functions.getGenieMapping(
            self.syn, "syn7434236")
        sex_mapping = process_functions.getGenieMapping(self.syn, "syn7434222")
        sampletype_mapping = process_functions.getGenieMapping(
            self.syn, "syn7434273")
        # Attach MSK to centers
        # clinicalMerged = clinicalMerged.fillna("")
        clinical = remap_clinical_values(
            clinical, sex_mapping, race_mapping, ethnicity_mapping,
            sampletype_mapping
        )
        remapped_clindf = clinical.apply(self.update_clinical, axis=1)
        # Some columns may have been added during update,
        # remove unwanted columns again
        keep_cols_idx = remapped_clindf.columns.isin(clinicalTemplate.columns)
        remapped_clindf = remapped_clindf.drop(
            remapped_clindf.columns[~keep_cols_idx], 1
        )
        remapped_clindf['CENTER'] = self.center
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
        patient_cols_table = self.syn.tableQuery(
            'select fieldName from syn8545211 where '
            'patient is True and inClinicalDb is True'
        )
        patient_cols = patient_cols_table.asDataFrame()['fieldName'].tolist()
        sample_cols_table = self.syn.tableQuery(
            'select fieldName from syn8545211 where '
            'sample is True and inClinicalDb is True'
        )
        sample_cols = sample_cols_table.asDataFrame()['fieldName'].tolist()
        clinicalTemplate = pd.DataFrame(
            columns=set(patient_cols + sample_cols)
        )
        sample = True
        patient = True

        return({'clinicalTemplate': clinicalTemplate,
                'sample': sample,
                'patient': patient,
                'patientCols': patient_cols,
                'sampleCols': sample_cols})

    def process_steps(self, clinicalDf,
                      databaseToSynIdMappingDf, newPath,
                      parentId, oncotree_link, clinicalTemplate,
                      sample, patient, patientCols, sampleCols):
        """Process clincial file, redact PHI values, upload to clinical
        database
        """
        patientdb_idx = databaseToSynIdMappingDf['Database'] == "patient"
        patient_synid = databaseToSynIdMappingDf.Id[patientdb_idx][0]
        sampledb_idx = databaseToSynIdMappingDf['Database'] == "sample"
        sample_synid = databaseToSynIdMappingDf.Id[sampledb_idx][0]

        newClinicalDf = self._process(clinicalDf, clinicalTemplate)
        newClinicalDf = redact_phi(newClinicalDf)

        if patient:
            patientClinical = newClinicalDf[
                patientCols].drop_duplicates("PATIENT_ID")
            self.uploadMissingData(patientClinical, "PATIENT_ID",
                                   patient_synid, parentId)

            process_functions.updateData(self.syn, patient_synid,
                                         patientClinical, self.center,
                                         col=patientCols, toDelete=True)
        if sample:
            if sum(newClinicalDf["SAMPLE_ID"].duplicated()) > 0:
                logger.error("There are duplicated samples, "
                             "and the duplicates are removed")
            sampleClinical = newClinicalDf[sampleCols].drop_duplicates(
                "SAMPLE_ID")
            # Exclude all clinical samples with wrong oncotree codes
            oncotree_mapping = pd.DataFrame()
            oncotree_mapping_dict = \
                process_functions.get_oncotree_code_mappings(oncotree_link)
            # Add in unknown key for oncotree code
            oncotree_mapping_dict['UNKNOWN'] = {}
            oncotree_mapping['ONCOTREE_CODE'] = list(oncotree_mapping_dict.keys())
            # Make oncotree codes uppercase (SpCC/SPCC)
            sampleClinical['ONCOTREE_CODE'] = sampleClinical[
                'ONCOTREE_CODE'].astype(str).str.upper()
            sampleClinical = sampleClinical[sampleClinical[
                'ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])]
            self.uploadMissingData(
                sampleClinical, "SAMPLE_ID", sample_synid, parentId)
            # ,retractedSampleSynId)
            process_functions.updateData(
                self.syn, sample_synid, sampleClinical,
                self.center, col=sampleCols, toDelete=True)

        newClinicalDf.to_csv(newPath, sep="\t", index=False)
        return(newPath)

    # VALIDATION
    def _validate(self, clinicaldf, oncotree_link):
        """
        This function validates the clinical file to make sure it adhere
        to the clinical SOP.

        Args:
            clinicalDF: Merged clinical file with patient and sample
                        information
            oncotree_link: Link to oncotree

        Returns:
            Error message
        """
        total_error = ""
        warning = ""

        clinicaldf.columns = [col.upper() for col in clinicaldf.columns]
        clinicaldf = clinicaldf.fillna("")

        # oncotree_mapping = process_functions.get_oncotree_codes(oncotree_link)
        # if oncotree_mapping.empty:
        oncotree_mapping_dict = \
            process_functions.get_oncotree_code_mappings(oncotree_link)
        oncotree_mapping = pd.DataFrame(
            {"ONCOTREE_CODE": list(oncotree_mapping_dict.keys())}
        )

        sampletype_mapping = \
            process_functions.getGenieMapping(self.syn, "syn7434273")

        ethnicity_mapping = \
            process_functions.getGenieMapping(self.syn, "syn7434242")

        race_mapping = \
            process_functions.getGenieMapping(self.syn, "syn7434236")

        sex_mapping = \
            process_functions.getGenieMapping(self.syn, "syn7434222")

        # CHECK: SAMPLE_ID
        sampleId = 'SAMPLE_ID'
        haveSampleColumn = \
            process_functions.checkColExist(clinicaldf, sampleId)

        if not haveSampleColumn:
            total_error += \
                "Sample Clinical File: Must have SAMPLE_ID column.\n"
        else:
            if sum(clinicaldf[sampleId].duplicated()) > 0:
                total_error += (
                    "Sample Clinical File: No duplicated SAMPLE_ID "
                    "allowed.\nIf there are no duplicated "
                    "SAMPLE_IDs, and both sample and patient files are "
                    "uploaded, then please check to make sure no duplicated "
                    "PATIENT_IDs exist in the patient clinical file.\n")
        # CHECK: PATIENT_ID
        patientId = "PATIENT_ID"
        # #CHECK: PATIENT_ID IN SAMPLE FILE
        havePatientColumn = \
            process_functions.checkColExist(clinicaldf, patientId)

        if not havePatientColumn:
            total_error += \
                "Patient Clinical File: Must have PATIENT_ID column.\n"

        # CHECK: within the sample file that the sample ids match
        # the patient ids
        if haveSampleColumn and havePatientColumn:
            # Make sure sample and patient ids are string cols
            clinicaldf[sampleId] = clinicaldf[sampleId].astype(str)
            clinicaldf[patientId] = clinicaldf[patientId].astype(str)
            if not all([patient in sample
                        for sample, patient in
                        zip(clinicaldf[sampleId], clinicaldf[patientId])]):

                total_error += (
                    "Sample Clinical File: PATIENT_ID's much be contained in "
                    "the SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n")
            # #CHECK: All samples must have associated patient data
            # (GENIE requires patient data)
            if not all(clinicaldf[patientId] != ""):
                total_error += (
                    "Patient Clinical File: All samples must have associated "
                    "patient information and no null patient ids allowed. "
                    "These samples are missing patient data: {}\n".format(
                        ", ".join(clinicaldf[sampleId][
                                  clinicaldf[patientId] == ""])))
            # CHECK: All patients should have associated sample data
            if not all(clinicaldf[sampleId] != ""):
                # ## MAKE WARNING FOR NOW###
                warning += (
                    "Sample Clinical File: All patients must have associated "
                    "sample information. These patients are missing sample "
                    "data: {}\n".format(
                        ", ".join(clinicaldf[patientId][
                                  clinicaldf[sampleId] == ""])))

        # CHECK: AGE_AT_SEQ_REPORT
        age = "AGE_AT_SEQ_REPORT"
        haveColumn = process_functions.checkColExist(clinicaldf, age)
        if haveColumn:
            # Deal with HIPAA converted rows from DFCI
            # First for loop can't int(text) because there
            # are instances that have <3435
            age_seq_report_df = clinicaldf[
                ~clinicaldf[age].isin(['Unknown', '>32485', '<6570'])
            ]

            # age_seq_report_df[age] = \
            #     remove_greaterthan_lessthan_str(age_seq_report_df[age])

            if not all([process_functions.checkInt(i)
                        for i in age_seq_report_df[age]]):
                total_error += (
                    "Sample Clinical File: Please double check your "
                    "AGE_AT_SEQ_REPORT. It must be an integer, 'Unknown', "
                    "'>32485', '<6570'.\n")
            else:
                age_seq_report_df[age] = age_seq_report_df[age].astype(int)
                median_age = age_seq_report_df[age].median()
                if median_age < 100:
                    total_error += (
                        "Sample Clinical File: Please double check your "
                        "AGE_AT_SEQ_REPORT. You may be reporting this value "
                        "in YEARS, please report in DAYS.\n")
        else:
            total_error += \
                "Sample Clinical File: Must have AGE_AT_SEQ_REPORT column.\n"

        # CHECK: ONCOTREE_CODE
        haveColumn = \
            process_functions.checkColExist(clinicaldf, "ONCOTREE_CODE")
        maleOncoCodes = ["TESTIS", "PROSTATE", "PENIS"]
        womenOncoCodes = ["CERVIX", "VULVA", "UTERUS", "OVARY"]
        if haveColumn:
            # Make oncotree codes uppercase (SpCC/SPCC)
            clinicaldf['ONCOTREE_CODE'] = \
                clinicaldf['ONCOTREE_CODE'].astype(str).str.upper()

            oncotree_codes = clinicaldf['ONCOTREE_CODE'][
                clinicaldf['ONCOTREE_CODE'] != "UNKNOWN"]

            if not all(oncotree_codes.isin(oncotree_mapping['ONCOTREE_CODE'])):
                unmapped_oncotrees = oncotree_codes[
                    ~oncotree_codes.isin(oncotree_mapping['ONCOTREE_CODE'])]
                total_error += (
                    "Sample Clinical File: Please double check that all your "
                    "ONCOTREE CODES exist in the mapping. You have {} samples "
                    "that don't map. These are the codes that "
                    "don't map: {}\n".format(
                        len(unmapped_oncotrees),
                        ",".join(set(unmapped_oncotrees))))
            # Should add the SEX mismatch into the dashboard file
            if process_functions.checkColExist(clinicaldf, "SEX") and \
               'oncotree_mapping_dict' in locals() and \
               havePatientColumn and \
               haveSampleColumn:

                wrongCodeSamples = []
                # This is to check if oncotree codes match the sex,
                # returns list of samples that have conflicting codes and sex
                for code, patient, sample in zip(
                        clinicaldf['ONCOTREE_CODE'],
                        clinicaldf['PATIENT_ID'],
                        clinicaldf['SAMPLE_ID']):

                    if oncotree_mapping_dict.get(code) is not None and \
                       sum(clinicaldf['PATIENT_ID'] == patient) > 0:

                        primaryCode = oncotree_mapping_dict[code][
                            'ONCOTREE_PRIMARY_NODE']

                        sex = clinicaldf['SEX'][
                            clinicaldf['PATIENT_ID'] == patient].values[0]
                        sex = float('nan') if sex == '' else float(sex)
                        if oncotree_mapping_dict[code][
                                'ONCOTREE_PRIMARY_NODE'] in maleOncoCodes and \
                           sex != 1.0:

                            wrongCodeSamples.append(sample)
                        if oncotree_mapping_dict[code][
                                'ONCOTREE_PRIMARY_NODE'] in womenOncoCodes and\
                           sex != 2.0:

                            wrongCodeSamples.append(sample)
                if len(wrongCodeSamples) > 0:
                    warning += (
                        "Sample Clinical File: Some SAMPLE_IDs have "
                        "conflicting SEX and ONCOTREE_CODES: {}\n".format(
                            ",".join(wrongCodeSamples)))
        else:
            total_error += \
                "Sample Clinical File: Must have ONCOTREE_CODE column.\n"

        warn, error = process_functions.check_col_and_values(
            clinicaldf, "SAMPLE_TYPE", sampletype_mapping['CODE'].tolist(),
            "Sample Clinical File", required=True)
        total_error += error

        # CHECK: SEQ_ASSAY_ID
        haveColumn = process_functions.checkColExist(clinicaldf,
                                                     "SEQ_ASSAY_ID")
        if haveColumn:
            if not all([i != "" for i in clinicaldf['SEQ_ASSAY_ID']]):
                total_error += (
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
                total_error += (
                    "Sample Clinical File: Please make sure your "
                    "SEQ_ASSAY_IDs start with your center "
                    "abbreviation: {}.\n".format(
                        ", ".join(invalid_seqassay)))
        else:
            total_error += \
                "Sample Clinical File: Must have SEQ_ASSAY_ID column.\n"

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
            clinicaldf['SEQ_DATE'] = [
                i.title() for i in clinicaldf['SEQ_DATE'].astype(str)
            ]

            seqdate = clinicaldf['SEQ_DATE'][
                clinicaldf['SEQ_DATE'] != 'Release']
            if sum(clinicaldf['SEQ_DATE'] == '') > 0:
                total_error += (
                    "Sample Clinical File: Samples without SEQ_DATEs will "
                    "NOT be released.\n")
            try:
                if not seqdate.empty:
                    seqdate.apply(
                        lambda date: datetime.datetime.strptime(date, '%b-%Y'))
                    if not seqdate.str.startswith(
                            ("Jan", "Apr", "Jul", "Oct")).all():
                        total_error += seq_date_error
            except ValueError:
                total_error += seq_date_error
        else:
            total_error += "Sample Clinical File: Must have SEQ_DATE column.\n"

        # CHECK: BIRTH_YEAR
        error = _check_year(clinicaldf=clinicaldf,
                            year_col="BIRTH_YEAR",
                            filename="Patient Clinical File",
                            allowed_string_values=['Unknown', '>89', '<18'])
        total_error += error

        # CHECK: YEAR DEATH
        error = _check_year(clinicaldf=clinicaldf,
                            year_col="YEAR_DEATH",
                            filename="Patient Clinical File",
                            allowed_string_values=['Unknown', 'Not Collected',
                                                   'Not Applicable',
                                                   '>89', '<18'])
        total_error += error

        # CHECK: YEAR CONTACT
        error = _check_year(clinicaldf=clinicaldf,
                            year_col="YEAR_CONTACT",
                            filename="Patient Clinical File",
                            allowed_string_values=['Unknown', 'Not Collected',
                                                   '>89', '<18'])
        total_error += error

        # CHECK: INT CONTACT
        haveColumn = process_functions.checkColExist(clinicaldf, "INT_CONTACT")
        if haveColumn:
            if not all([process_functions.checkInt(i)
                        for i in clinicaldf.INT_CONTACT if i not in
                        ['>32485', '<6570', 'Unknown', 'Not Collected']]):

                total_error += (
                    "Patient Clinical File: Please double check your "
                    "INT_CONTACT column, it must be an integer, '>32485', "
                    "'<6570', 'Unknown' or 'Not Collected'.\n")
        else:
            total_error += \
                "Patient Clinical File: Must have INT_CONTACT column.\n"

        # INT DOD
        haveColumn = process_functions.checkColExist(clinicaldf, "INT_DOD")
        if haveColumn:
            if not all([process_functions.checkInt(i)
                        for i in clinicaldf.INT_DOD if i not in
                        ['>32485', '<6570', 'Unknown',
                         'Not Collected', 'Not Applicable']]):

                total_error += (
                    "Patient Clinical File: Please double check your INT_DOD "
                    "column, it must be an integer, '>32485', '<6570', "
                    "'Unknown', 'Not Collected' or 'Not Applicable'.\n")
        else:
            total_error += \
                "Patient Clinical File: Must have INT_DOD column.\n"

        haveColumn = process_functions.checkColExist(clinicaldf, "DEAD")
        if haveColumn:
            # Need to have check_bool function
            if not all([
                    str(i).upper() in ['TRUE', 'FALSE']
                    for i in clinicaldf.DEAD if i not in
                    ['Unknown', 'Not Collected']]):
                total_error += (
                    "Patient Clinical File: Please double check your "
                    "DEAD column, it must be True, False, 'Unknown' or "
                    "'Not Collected'.\n")
        else:
            total_error += \
                "Patient Clinical File: Must have DEAD column.\n"

        # CHECK: PRIMARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicaldf, "PRIMARY_RACE", race_mapping['CODE'].tolist(),
            "Patient Clinical File"
        )
        warning += warn
        total_error += error

        # CHECK: SECONDARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicaldf, "SECONDARY_RACE", race_mapping['CODE'].tolist(),
            "Patient Clinical File"
        )
        warning += warn
        total_error += error

        # CHECK: TERTIARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicaldf, "TERTIARY_RACE", race_mapping['CODE'].tolist(),
            "Patient Clinical File"
        )
        warning += warn
        total_error += error

        # CHECK: SEX
        warn, error = process_functions.check_col_and_values(
            clinicaldf, "SEX", sex_mapping['CODE'].tolist(),
            "Patient Clinical File", required=True
        )
        warning += warn
        total_error += error

        # CHECK: ETHNICITY
        warn, error = process_functions.check_col_and_values(
            clinicaldf, "ETHNICITY", ethnicity_mapping['CODE'].tolist(),
            "Patient Clinical File"
        )
        warning += warn
        total_error += error

        return total_error, warning

    def _get_dataframe(self, filePathList):
        clinicaldf = pd.read_csv(filePathList[0], sep="\t", comment="#")
        if len(filePathList) > 1:
            other_clinicaldf = pd.read_csv(filePathList[1], sep="\t",
                                           comment="#")
            try:
                clinicaldf = clinicaldf.merge(other_clinicaldf, on="PATIENT_ID")
            except Exception:
                raise ValueError((
                    "If submitting separate patient and sample files, "
                    "they both must have the PATIENT_ID column"))
            # Must figure out which is sample and which is patient
            if "sample" in filePathList[0]:
                sample = clinicaldf
                patient = other_clinicaldf
            else:
                sample = other_clinicaldf
                patient = clinicaldf

            if not all(sample['PATIENT_ID'].isin(patient['PATIENT_ID'])):
                raise ValueError((
                    "Patient Clinical File: All samples must have associated "
                    "patient information"))

        return clinicaldf
