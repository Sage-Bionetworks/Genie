import datetime
import os
import logging
import subprocess
import yaml

import pandas as pd
import synapseclient

from .example_filetype_format import FileTypeFormat
from . import process_functions
from .database_to_staging import redact_phi

logger = logging.getLogger(__name__)


def remove_greaterthan_lessthan_str(col):
    '''
    In clinical file, there are redacted value such as >89 and <17.
    These < and > signs must be removed
    '''
    try:
        col = [
            text.replace(">", "")
            if isinstance(text, str) else text for text in col]
        col = [
            int(text.replace("<", ""))
            if isinstance(text, str) and text != "" else text
            for text in col]
    except ValueError:
        pass
    return(col)


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

    # PROCESSING
    # Update clinical file with the correct mappings
    def update_clinical(self, x, sex_mapping,
                        race_mapping, ethnicity_mapping, sample_type):
        # PATIENT ID
        if x.get("PATIENT_ID") is not None:
            x['PATIENT_ID'] = process_functions.checkGenieId(
                x['PATIENT_ID'], self.center)
        # RACE
        if x.get('PRIMARY_RACE') is not None:
            x['PRIMARY_RACE'] = process_functions.getCODE(
                race_mapping, x['PRIMARY_RACE'])
        else:
            x['PRIMARY_RACE'] = "Not Collected"

        if x.get('SECONDARY_RACE') is not None:
            x['SECONDARY_RACE'] = process_functions.getCODE(
                race_mapping, x['SECONDARY_RACE'])
        else:
            x['SECONDARY_RACE'] = "Not Collected"

        if x.get('TERTIARY_RACE') is not None:
            x['TERTIARY_RACE'] = process_functions.getCODE(
                race_mapping, x['TERTIARY_RACE'])
        else:
            x['TERTIARY_RACE'] = "Not Collected"
        # ETHNICITY
        if x.get('ETHNICITY') is not None:
            x['ETHNICITY'] = process_functions.getCODE(
                ethnicity_mapping, x['ETHNICITY'])
        else:
            x['ETHNICITY'] = "Not Collected"
        # BIRTH YEAR
        if x.get("BIRTH_YEAR") is not None:
            # BIRTH YEAR (Check if integer)
            if process_functions.checkInt(x['BIRTH_YEAR']):
                x['BIRTH_YEAR'] = int(x['BIRTH_YEAR'])
        # SEX
        if x.get("SEX") is not None:
            x['SEX'] = process_functions.getCODE(sex_mapping, x['SEX'])
        # TRIM EVERY COLUMN MAKE ALL DASHES
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

        # SAMPLE_TYPE
        if x.get('SAMPLE_TYPE') is not None:
            sampleType = x['SAMPLE_TYPE']
            x['SAMPLE_TYPE'] = process_functions.getCODE(
                sample_type, sampleType)
            # Trim spaces
            x['SAMPLE_TYPE_DETAILED'] = process_functions.getCODE(
                sample_type, sampleType, useDescription=True)

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
        return(x)

    def uploadMissingData(self, df, col, dbSynId, stagingSynId,
                          retractionSynId=None):
        samples = "','".join(df[col])
        path = os.path.join(
            process_functions.SCRIPT_DIR,
            "{}_missing_{}.csv".format(self._fileType, col))
        missing = self.syn.tableQuery(
            "select {} from {} where CENTER='{}' and {} not in ('{}')".format(
                col, dbSynId, self.center, col, samples))
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
        sampleType_mapping = process_functions.getGenieMapping(
            self.syn, "syn7434273")

        # Attach MSK to centers
        # clinicalMerged = clinicalMerged.fillna("")
        clinicalRemapped = clinical.apply(lambda x: self.update_clinical(
            x, sex_mapping, race_mapping,
            ethnicity_mapping, sampleType_mapping), 1)
        # Some columns may have been added during update,
        # remove unwanted columns again
        clinicalRemapped = clinicalRemapped.drop(clinicalRemapped.columns[
            ~clinicalRemapped.columns.isin(clinicalTemplate.columns)], 1)

        clinicalRemapped['CENTER'] = self.center

        return(clinicalRemapped)

    def preprocess(self, newpath):
        '''
        Gather preprocess parameters

        Args:
            filePath: Path to file

        Returns:
            dict with keys - 'clinicalTemplate', 'sample', 'patient',
                             'patientCols', 'sampleCols'
        '''
        entity_name = os.path.basename(newpath)
        # These synapse ids for the clinical tier release scope is
        # hardcoded because it never changes
        patientColsTable = self.syn.tableQuery(
            'select fieldName from syn8545211 where patient is '
            'True and inClinicalDb is True')
        patientCols = patientColsTable.asDataFrame()['fieldName'].tolist()
        sampleColsTable = self.syn.tableQuery(
            'select fieldName from syn8545211 where sample is True '
            'and inClinicalDb is True')
        sampleCols = sampleColsTable.asDataFrame()['fieldName'].tolist()

        # if "patient" in filepath.lower():
        #     clinicalTemplate = pd.DataFrame(columns=patientCols)
        #     sample = False
        #     patient = True
        # elif "sample" in filepath.lower():
        #     clinicalTemplate = pd.DataFrame(columns=sampleCols)
        #     sample = True
        #     patient = False
        # else:
        clinicalTemplate = pd.DataFrame(columns=set(patientCols + sampleCols))
        sample = True
        patient = True

        return({'clinicalTemplate': clinicalTemplate,
                'sample': sample,
                'patient': patient,
                'patientCols': patientCols,
                'sampleCols': sampleCols})

    def process_steps(self, clinicalDf,
                      databaseToSynIdMappingDf, newPath,
                      parentId, oncotree_link, clinicalTemplate,
                      sample, patient, patientCols, sampleCols):
        patientSynId = databaseToSynIdMappingDf.Id[
            databaseToSynIdMappingDf['Database'] == "patient"][0]
        sampleSynId = databaseToSynIdMappingDf.Id[
            databaseToSynIdMappingDf['Database'] == "sample"][0]

        newClinicalDf = self._process(clinicalDf, clinicalTemplate)
        newClinicalDf = redact_phi(newClinicalDf)
        if patient:
            patientClinical = newClinicalDf[
                patientCols].drop_duplicates("PATIENT_ID")
            self.uploadMissingData(
                patientClinical, "PATIENT_ID",
                patientSynId, parentId)
            # retractedPatientSynId)
            process_functions.updateData(
                self.syn, patientSynId, patientClinical,
                self.center, col=patientCols, toDelete=True)
        if sample:
            if sum(newClinicalDf["SAMPLE_ID"].duplicated()) > 0:
                logger.error(
                    "There are duplicated samples, "
                    "and the duplicates are removed")
            sampleClinical = newClinicalDf[sampleCols].drop_duplicates(
                "SAMPLE_ID")
            # Exclude all clinical samples with wrong oncotree codes
            oncotree_mapping = pd.DataFrame()
            oncotree_mapping_dict = \
                process_functions.get_oncotree_code_mappings(oncotree_link)
            # Add in unknown key for oncotree code
            oncotree_mapping_dict['UNKNOWN'] = {}
            oncotree_mapping['ONCOTREE_CODE'] = oncotree_mapping_dict.keys()
            # Make oncotree codes uppercase (SpCC/SPCC)
            sampleClinical['ONCOTREE_CODE'] = sampleClinical[
                'ONCOTREE_CODE'].astype(str).str.upper()
            sampleClinical = sampleClinical[sampleClinical[
                'ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])]
            self.uploadMissingData(
                sampleClinical, "SAMPLE_ID", sampleSynId, parentId)
            # ,retractedSampleSynId)
            process_functions.updateData(
                self.syn, sampleSynId, sampleClinical,
                self.center, col=sampleCols, toDelete=True)

        newClinicalDf.to_csv(newPath, sep="\t", index=False)
        return(newPath)

    # VALIDATION
    def _validate(self, clinicalDF, oncotree_link):
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

        clinicalDF.columns = [col.upper() for col in clinicalDF.columns]
        clinicalDF = clinicalDF.fillna("")

        # oncotree_mapping = process_functions.get_oncotree_codes(oncotree_link)
        # if oncotree_mapping.empty:
        oncotree_mapping = pd.DataFrame()
        oncotree_mapping_dict = \
            process_functions.get_oncotree_code_mappings(oncotree_link)
        oncotree_mapping['ONCOTREE_CODE'] = oncotree_mapping_dict.keys()

        sampleType_mapping = \
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
            process_functions.checkColExist(clinicalDF, sampleId)

        if not haveSampleColumn:
            total_error += \
                "Sample Clinical File: Must have SAMPLE_ID column.\n"
        else:
            if sum(clinicalDF[sampleId].duplicated()) > 0:
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
            process_functions.checkColExist(clinicalDF, patientId)

        if not havePatientColumn:
            total_error += \
                "Patient Clinical File: Must have PATIENT_ID column.\n"

        # CHECK: within the sample file that the sample ids match
        # the patient ids
        if haveSampleColumn and havePatientColumn:
            # Make sure sample and patient ids are string cols
            clinicalDF[sampleId] = clinicalDF[sampleId].astype(str)
            clinicalDF[patientId] = clinicalDF[patientId].astype(str)
            if not all([patient in sample
                        for sample, patient in
                        zip(clinicalDF[sampleId], clinicalDF[patientId])]):

                total_error += (
                    "Sample Clinical File: PATIENT_ID's much be contained in "
                    "the SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n")
            # #CHECK: All samples must have associated patient data
            # (GENIE requires patient data)
            if not all(clinicalDF[patientId] != ""):
                total_error += (
                    "Patient Clinical File: All samples must have associated "
                    "patient information and no null patient ids allowed. "
                    "These samples are missing patient data: {}\n".format(
                        ", ".join(clinicalDF[sampleId][
                                  clinicalDF[patientId] == ""])))
            # CHECK: All patients should have associated sample data
            if not all(clinicalDF[sampleId] != ""):
                # ## MAKE WARNING FOR NOW###
                warning += (
                    "Sample Clinical File: All patients must have associated "
                    "sample information. These patients are missing sample "
                    "data: {}\n".format(
                        ", ".join(clinicalDF[patientId][
                                  clinicalDF[sampleId] == ""])))

        # CHECK: AGE_AT_SEQ_REPORT
        age = "AGE_AT_SEQ_REPORT"
        haveColumn = process_functions.checkColExist(clinicalDF, age)
        if haveColumn:
            # Deal with HIPAA converted rows from DFCI
            # First for loop can't int(text) because there
            # are instances that have <3435
            age_seq_report_df = \
                clinicalDF[~clinicalDF[age].isin(['Unknown'])]

            age_seq_report_df[age] = \
                remove_greaterthan_lessthan_str(age_seq_report_df[age])

            if not all([process_functions.checkInt(i)
                        for i in age_seq_report_df[age]]):
                total_error += (
                    "Sample Clinical File: Please double check your "
                    "AGE_AT_SEQ_REPORT. It must be an integer or 'Unknown'.\n")
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
            process_functions.checkColExist(clinicalDF, "ONCOTREE_CODE")
        maleOncoCodes = ["TESTIS", "PROSTATE", "PENIS"]
        womenOncoCodes = ["CERVIX", "VULVA", "UTERUS", "OVARY"]
        if haveColumn:
            # Make oncotree codes uppercase (SpCC/SPCC)
            clinicalDF['ONCOTREE_CODE'] = \
                clinicalDF['ONCOTREE_CODE'].astype(str).str.upper()

            oncotree_codes = clinicalDF['ONCOTREE_CODE'][
                clinicalDF['ONCOTREE_CODE'] != "UNKNOWN"]

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

            if process_functions.checkColExist(clinicalDF, "SEX") and \
               'oncotree_mapping_dict' in locals() and \
               havePatientColumn and \
               haveSampleColumn:

                wrongCodeSamples = []
                # This is to check if oncotree codes match the sex,
                # returns list of samples that have conflicting codes and sex
                for code, patient, sample in zip(
                        clinicalDF['ONCOTREE_CODE'],
                        clinicalDF['PATIENT_ID'],
                        clinicalDF['SAMPLE_ID']):

                    if oncotree_mapping_dict.get(code) is not None and \
                       sum(clinicalDF['PATIENT_ID'] == patient) > 0:

                        primaryCode = oncotree_mapping_dict[code][
                            'ONCOTREE_PRIMARY_NODE']

                        sex = clinicalDF['SEX'][
                            clinicalDF['PATIENT_ID'] == patient].values[0]
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
            clinicalDF, "SAMPLE_TYPE", sampleType_mapping['CODE'].tolist(),
            "Sample Clinical File", required=True)
        total_error += error

        # CHECK: SEQ_ASSAY_ID
        haveColumn = \
            process_functions.checkColExist(clinicalDF, "SEQ_ASSAY_ID")
        if haveColumn:
            if not all([i != "" for i in clinicalDF['SEQ_ASSAY_ID']]):
                total_error += (
                    "Sample Clinical File: Please double check your "
                    "SEQ_ASSAY_ID columns, there are empty rows.\n")
            # must remove empty seq assay ids first
            # Checking if seq assay ids start with the center name
            seqAssayIds = \
                clinicalDF.SEQ_ASSAY_ID[clinicalDF.SEQ_ASSAY_ID != ""]
            allSeqAssays = seqAssayIds.unique()
            notNormalized = []
            not_caps = []
            for seqassay in allSeqAssays:
                # SEQ Ids are all capitalized now, so no need to check
                # for differences in case
                if not seqassay.upper().startswith(self.center):
                    not_caps.append(seqassay)
            if len(not_caps) > 0:
                total_error += (
                    "Sample Clinical File: Please make sure your "
                    "SEQ_ASSAY_IDs start with your center "
                    "abbreviation: {}.\n".format(
                        ", ".join(not_caps)))
        else:
            total_error += \
                "Sample Clinical File: Must have SEQ_ASSAY_ID column.\n"

        haveColumn = process_functions.checkColExist(clinicalDF, "SEQ_DATE")
        seq_date_error = (
            "Sample Clinical File: SEQ_DATE must be one of five values- "
            "For Jan-March: use Jan-YEAR. "
            "For Apr-June: use Apr-YEAR. "
            "For July-Sep: use Jul-YEAR. "
            "For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) "
            "For values that don't have SEQ_DATES that "
            "you want released use 'release'.\n")

        if haveColumn:
            clinicalDF['SEQ_DATE'] = [
                i.title() for i in clinicalDF['SEQ_DATE'].astype(str)]

            seqDate = clinicalDF['SEQ_DATE'][
                clinicalDF['SEQ_DATE'] != 'Release']
            if sum(clinicalDF['SEQ_DATE'] == '') > 0:
                total_error += (
                    "Sample Clinical File: Samples without SEQ_DATEs will "
                    "NOT be released.\n")
            try:
                if not seqDate.empty:
                    dates = seqDate.apply(
                        lambda date: datetime.datetime.strptime(date, '%b-%Y'))
                    # REMOVE JUN LATER
                    if not all([i.startswith(("Jul", "Jan", "Oct", "Apr"))
                                for i in seqDate]):
                        total_error += seq_date_error
            except ValueError:
                total_error += seq_date_error
        else:
            total_error += "Sample Clinical File: Must have SEQ_DATE column.\n"

        # CHECK: BIRTH_YEAR
        birth_year = "BIRTH_YEAR"
        haveColumn = process_functions.checkColExist(clinicalDF, birth_year)
        if haveColumn:
            birth_year_df = \
                clinicalDF[~clinicalDF[birth_year].isin(['Unknown'])]
            # Deal with HIPAA converted rows from DFCI
            # First for loop can't int(text) because there are
            # instances that have <YYYY
            birth_year_df[birth_year] = \
                remove_greaterthan_lessthan_str(birth_year_df[birth_year])

            try:
                years = birth_year_df[birth_year].apply(
                    lambda x: datetime.datetime.strptime(
                        str(int(x)), '%Y').year >
                    datetime.datetime.utcnow().year)

                assert not years.any()
            except Exception:
                total_error += (
                    "Patient Clinical File: Please double check your "
                    "BIRTH_YEAR column, it must be an integer in YYYY format "
                    "> {year} or 'Unknown'.\n".format(
                        year=datetime.datetime.utcnow().year))
        else:
            total_error += \
                "Patient Clinical File: Must have BIRTH_YEAR column.\n"

        # CHECK: VITAL_STATUS
        # YEAR DEATH
        haveColumn = process_functions.checkColExist(clinicalDF, "YEAR_DEATH")
        if haveColumn:
            notNullYears = clinicalDF.YEAR_DEATH[~clinicalDF.YEAR_DEATH.isin(
                ['Unknown', 'Not Collected', 'Not Applicable'])]
            try:
                notNullYears.apply(
                    lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
            except Exception:
                total_error += (
                    "Patient Clinical File: Please double check your "
                    "YEAR_DEATH column, it must be an integer in YYYY format, "
                    "'Unknown', 'Not Applicable' or 'Not Collected'.\n")
        else:
            total_error += \
                "Patient Clinical File: Must have YEAR_DEATH column.\n"

        # YEAR CONTACT
        haveColumn = process_functions.checkColExist(
            clinicalDF, "YEAR_CONTACT")
        if haveColumn:
            notNullYears = clinicalDF.YEAR_CONTACT[
                ~clinicalDF.YEAR_CONTACT.isin(['Unknown', 'Not Collected'])]
            try:
                notNullYears.apply(
                    lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
            except Exception:
                total_error += (
                    "Patient Clinical File: Please double check your "
                    "YEAR_CONTACT column, it must be an integer in YYYY "
                    "format, 'Unknown' or 'Not Collected'.\n")
        else:
            total_error += \
                "Patient Clinical File: Must have YEAR_CONTACT column.\n"

        # INT CONTACT
        haveColumn = process_functions.checkColExist(clinicalDF, "INT_CONTACT")
        if haveColumn:
            if not all([
                    process_functions.checkInt(i)
                    for i in clinicalDF.INT_CONTACT if i not in
                    ['>32485', '<6570', 'Unknown', 'Not Collected']]):

                total_error += (
                    "Patient Clinical File: Please double check your "
                    "INT_CONTACT column, it must be an integer, '>32485', "
                    "'<6570', 'Unknown' or 'Not Collected'.\n")
        else:
            total_error += \
                "Patient Clinical File: Must have INT_CONTACT column.\n"

        # INT DOD
        haveColumn = process_functions.checkColExist(clinicalDF, "INT_DOD")
        if haveColumn:
            if not all([
                    process_functions.checkInt(i)
                    for i in clinicalDF.INT_DOD if i not in
                    ['>32485', '<6570', 'Unknown',
                     'Not Collected', 'Not Applicable']]):

                total_error += (
                    "Patient Clinical File: Please double check your INT_DOD "
                    "column, it must be an integer, '>32485', '<6570', "
                    "'Unknown', 'Not Collected' or 'Not Applicable'.\n")
        else:
            total_error += \
                "Patient Clinical File: Must have INT_DOD column.\n"

        haveColumn = process_functions.checkColExist(clinicalDF, "DEAD")
        if haveColumn:
            # Need to have check_bool function
            if not all([
                    str(i).upper() in ['TRUE', 'FALSE']
                    for i in clinicalDF.DEAD if i not in
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
            clinicalDF, "PRIMARY_RACE", race_mapping['CODE'].tolist(),
            "Patient Clinical File")
        warning += warn
        total_error += error

        # CHECK: SECONDARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicalDF, "SECONDARY_RACE", race_mapping['CODE'].tolist(),
            "Patient Clinical File")
        warning += warn
        total_error += error

        # CHECK: TERTIARY_RACE
        warn, error = process_functions.check_col_and_values(
            clinicalDF, "TERTIARY_RACE", race_mapping['CODE'].tolist(),
            "Patient Clinical File")
        warning += warn
        total_error += error

        # CHECK: SEX
        warn, error = process_functions.check_col_and_values(
            clinicalDF, "SEX", sex_mapping['CODE'].tolist(),
            "Patient Clinical File", required=True)
        warning += warn
        total_error += error

        # CHECK: ETHNICITY
        warn, error = process_functions.check_col_and_values(
            clinicalDF, "ETHNICITY", ethnicity_mapping['CODE'].tolist(),
            "Patient Clinical File")
        warning += warn
        total_error += error

        return(total_error, warning)

    def _get_dataframe(self, filePathList):
        clinicalDf = pd.read_csv(filePathList[0], sep="\t", comment="#")
        if len(filePathList) > 1:
            otherClinicalDf = pd.read_csv(
                filePathList[1], sep="\t", comment="#")
            try:
                clinicalDf = clinicalDf.merge(otherClinicalDf, on="PATIENT_ID")
            except Exception:
                raise ValueError((
                    "If submitting separate patient and sample files, "
                    "they both must have the PATIENT_ID column"))
            # Must figure out which is sample and which is patient
            if "sample" in filePathList[0]:
                sample = clinicalDf
                patient = otherClinicalDf
            else:
                sample = otherClinicalDf
                patient = clinicalDf

            if not all(sample['PATIENT_ID'].isin(patient['PATIENT_ID'])):
                raise ValueError((
                    "Patient Clinical File: All samples must have associated "
                    "patient information"))

        return(clinicalDf)
