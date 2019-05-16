from __future__ import absolute_import
from genie import FileTypeFormat, process_functions
import os
import logging
import pandas as pd
import synapseclient
# import re
import datetime
logger = logging.getLogger(__name__)


def checkMapping(
        clinicalDF, colName, mapping, required=False, fileType="Patient"):
    """
    This function checks if the column exists then checks if the
    values in the column have the correct integer values

    :params clinicalDF          Patient/sample/flattened clinical file
    :params colName:            Expected column name
    :params mapping:            List of possible values

    :returns:                   A tuple warning, error
    """
    warning = ""
    error = ""
    haveColumn = process_functions.checkColExist(clinicalDF, colName)
    if not haveColumn:
        if required:
            error = "{}: clinical file must have {} column.\n".format(
                fileType, colName)
        else:
            warning = (
                "{}: clinical file doesn't have {} column. "
                "A blank column will be added\n".format(fileType, colName))
    else:
        if not all([i in mapping.tolist() for i in clinicalDF[colName]]):
            error = (
                "{}: Please double check your {} column.  "
                "This column must be these values {}or blank.\n".format(
                    fileType,
                    colName,
                    ", ".join(map(str, mapping)).replace(".0", "")))
    return(warning, error)


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


class clinical_individual(FileTypeFormat):

    _fileType = "clinical"

    _process_kwargs = [
        "newPath", "parentId", "databaseToSynIdMappingDf"]

    # VALIDATE FILE NAME
    def _validateFilename(self, filePath):
        assert os.path.basename(filePath) == "clinical_individual.csv", \
            "Individual clinical filename is not correct."

    # PROCESSING
    def uploadMissingData(self, df, col, dbSynId, 
                          stagingSynId, retractionSynId=None):
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
        clinical = clinical.fillna("")
        clinical['center'] = self.center

        return(clinical)

    def process_steps(
            self, filePath,
            databaseToSynIdMappingDf, newPath,
            parentId, oncotreeLink):
        patientSynId = databaseToSynIdMappingDf.Id[
            databaseToSynIdMappingDf['Database'] == "patient"][0]

        clinicalDf = pd.read_csv(filePath, sep="\t", comment="#")

        if "patient" in filePath.lower():
            clinicalTemplate = pd.DataFrame(columns=patientCols)
            patient = True
        elif "sample" in filePath.lower():
            clinicalTemplate = pd.DataFrame(columns=sampleCols)
            sample = True
        else:
            clinicalTemplate = pd.DataFrame(
                columns=set(patientCols + sampleCols))
            sample = True
            patient = True

        newClinicalDf = self._process(clinicalDf, clinicalTemplate)

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
                process_functions.get_oncotree_code_mappings(oncotreeLink)
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
    def _validate(self, clinicalDF, oncotreeLink):
        """
        This function validates the clinical file to make sure it adhere
        to the clinical SOP.

        Args:
            clinicalDF: Merged clinical file with patient and sample
                        information
            oncotreeLink: Link to oncotree

        Returns:
            Error message
        """
        total_error = ""
        warning = ""

        clinicalDF.columns = [col.upper() for col in clinicalDF.columns]
        clinicalDF = clinicalDF.fillna("")

        # oncotree_mapping = process_functions.get_oncotree_codes(oncotreeLink)
        # if oncotree_mapping.empty:
        oncotree_mapping = pd.DataFrame()
        oncotree_mapping_dict = \
            process_functions.get_oncotree_code_mappings(oncotreeLink)
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
                "Sample: clinical file must have SAMPLE_ID column.\n"
        else:
            if sum(clinicalDF[sampleId].duplicated()) > 0:
                total_error += (
                    "Sample: No duplicated SAMPLE_ID in the sample file "
                    "allowed.\nIf there are no duplicated SAMPLE_IDs, and "
                    "both sample and patient files are uploaded, then please "
                    "check to make sure no duplicated PATIENT_IDs exist in "
                    "the patient file.\n")

        # CHECK: PATIENT_ID
        patientId = "PATIENT_ID"
        # #CHECK: PATIENT_ID IN SAMPLE FILE
        havePatientColumn = \
            process_functions.checkColExist(clinicalDF, patientId)

        if not havePatientColumn:
            total_error += \
                "Patient: clinical file must have PATIENT_ID column.\n"

        # CHECK: within the sample file that the sample ids match
        # the patient ids
        if haveSampleColumn and havePatientColumn:
            if not all([patient in sample
                        for sample, patient in
                        zip(clinicalDF[sampleId], clinicalDF[patientId])]):

                total_error += (
                    "Sample: PATIENT_ID's much be contained in the "
                    "SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n")
            # #CHECK: All samples must have associated patient data
            # (GENIE requires patient data)
            if not all(clinicalDF[patientId] != ""):
                total_error += (
                    "Patient: All samples must have associated patient "
                    "information and no null patient ids allowed. These "
                    "samples are missing patient data: {}\n".format(
                        ", ".join(clinicalDF[sampleId][
                                  clinicalDF[patientId] == ""])))
            # CHECK: All patients should have associated sample data
            if not all(clinicalDF[sampleId] != ""):
                # ## MAKE WARNING FOR NOW###
                warning += (
                    "Sample: All patients must have associated sample "
                    "information. These patients are missing sample "
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
                clinicalDF[~clinicalDF[age].isin(['Unknown', ''])]

            age_seq_report_df[age] = \
                remove_greaterthan_lessthan_str(age_seq_report_df[age])

            if not all([process_functions.checkInt(i)
                        for i in age_seq_report_df[age]]):
                total_error += (
                    "Sample: Please double check your AGE_AT_SEQ_REPORT.  "
                    "It must be an integer or 'Unknown'.\n")
            else:
                age_seq_report_df[age] = age_seq_report_df[age].astype(int)
                median_age = pd.np.median(age_seq_report_df[age])
                if median_age < 100:
                    total_error += (
                        "Sample: Please double check your AGE_AT_SEQ_REPORT.  "
                        "You may be reporting this value in YEARS, "
                        "please report in DAYS.\n")
        else:
            total_error += \
                "Sample: clinical file must have AGE_AT_SEQ_REPORT column.\n"

        # CHECK: SAMPLE_TYPE
        haveColumn = process_functions.checkColExist(clinicalDF, "SAMPLE_TYPE")
        if haveColumn:
            if clinicalDF.SAMPLE_TYPE.dtype == int:
                if not all(clinicalDF['SAMPLE_TYPE'].isin(
                        sampleType_mapping['CODE'])):
                    total_error += (
                        "Sample: Please double check your SAMPLE_TYPE column. "
                        "This column must be {}.\n".format(
                            ", ".join(map(str, sampleType_mapping['CODE']))))
            else:
                total_error += (
                    "Sample: Please double check your SAMPLE_TYPE column. "
                    "No null values allowed.\n")

        else:
            total_error += \
                "Sample: clinical file must have SAMPLE_TYPE column.\n"

        # CHECK: SEQ_ASSAY_ID
        haveColumn = \
            process_functions.checkColExist(clinicalDF, "SEQ_ASSAY_ID")
        if haveColumn:
            if not all([i != "" for i in clinicalDF['SEQ_ASSAY_ID']]):
                total_error += (
                    "Sample: Please double check your SEQ_ASSAY_ID columns, "
                    "there are empty rows.\n")
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
                    "Sample: Please make sure your SEQ_ASSAY_IDs start with "
                    "your center abbreviation: {}.\n".format(
                        ", ".join(not_caps)))
        else:
            total_error += \
                "Sample: clinical file must have SEQ_ASSAY_ID column.\n"

        haveColumn = process_functions.checkColExist(clinicalDF, "SEQ_DATE")
        seq_date_error = (
            "Sample: SEQ_DATE must be one of five values- "
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
                total_error += \
                    "Sample: Samples without SEQ_DATEs will NOT be released.\n"
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
            total_error += "Sample: clinical file must SEQ_DATE column\n"

        # CHECK: BIRTH_YEAR
        birth_year = "BIRTH_YEAR"
        haveColumn = process_functions.checkColExist(clinicalDF, birth_year)
        if haveColumn:
            # Deal with HIPAA converted rows from DFCI
            # First for loop can't int(text) because there are
            # instances that have <YYYY
            # Remove '' for blank value support
            birth_year_df = \
                clinicalDF[~clinicalDF[birth_year].isin(['Unknown', ''])]
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
                    "Patient: Please double check your BIRTH_YEAR column, "
                    "it must be an integer in YYYY format > {year} or "
                    "'Unknown'.  Support for blank values will be deprecated "
                    "in 7...releases.\n".format(
                        year=datetime.datetime.utcnow().year))
        else:
            total_error += \
                "Patient: clinical file must have BIRTH_YEAR column.\n"

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
                    "Patient: Please double check your YEAR_DEATH column, "
                    "it must be an integer in YYYY format, "
                    "'Unknown', 'Not Applicable' or 'Not Collected'.\n")
        else:
            warning += (
                "Patient: Must have YEAR_DEATH column for "
                "7...release uploads.\n")

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
                    "Patient: Please double check your YEAR_CONTACT column, "
                    "it must be an integer in YYYY format, "
                    "'Unknown' or 'Not Collected'.\n")
        else:
            warning += (
                "Patient: Must have YEAR_CONTACT column for "
                "7...release uploads.\n")

        # INT CONTACT
        haveColumn = process_functions.checkColExist(clinicalDF, "INT_CONTACT")
        if haveColumn:
            if not all([
                    process_functions.checkInt(i)
                    for i in clinicalDF.INT_CONTACT if i not in
                    ['>32485', '<6570', 'Unknown', 'Not Collected']]):

                total_error += (
                    "Patient: Please double check your INT_CONTACT column, "
                    "it must be an integer, '>32485', '<6570', 'Unknown' "
                    "or 'Not Collected'.\n")
        else:
            warning += (
                "Patient: Must have INT_CONTACT column for "
                "7...release uploads.\n")

        # INT DOD
        haveColumn = process_functions.checkColExist(clinicalDF, "INT_DOD")
        if haveColumn:
            if not all([
                    process_functions.checkInt(i)
                    for i in clinicalDF.INT_DOD if i not in
                    ['>32485', '<6570', 'Unknown',
                     'Not Collected', 'Not Applicable']]):

                total_error += (
                    "Patient: Please double check your INT_DOD column, "
                    "it must be an integer, '>32485', '<6570', 'Unknown', "
                    "'Not Collected' or 'Not Applicable'.\n")
        else:
            warning += \
                "Patient: Must have INT_DOD column for 7...release uploads.\n"

        haveColumn = process_functions.checkColExist(clinicalDF, "DEAD")
        if haveColumn:
            # Need to have check_bool function
            if not all([
                    str(i).upper() in ['TRUE', 'FALSE']
                    for i in clinicalDF.DEAD if i not in
                    ['Unknown', 'Not Collected']]):
                total_error += (
                    "Patient: Please double check your DEAD column, "
                    "it must be True, False, 'Unknown' or 'Not Collected'.\n")
        else:
            warning += \
                "Patient: Must have DEAD column for 7...release uploads.\n"

        # CHECK: PRIMARY_RACE
        warn, error = checkMapping(
            clinicalDF, "PRIMARY_RACE", race_mapping['CODE'])
        warning += warn
        total_error += error

        # CHECK: SECONDARY_RACE
        warn, error = checkMapping(
            clinicalDF, "SECONDARY_RACE", race_mapping['CODE'])
        warning += warn
        total_error += error

        # CHECK: TERTIARY_RACE
        warn, error = checkMapping(
            clinicalDF, "TERTIARY_RACE", race_mapping['CODE'])
        warning += warn
        total_error += error

        # CHECK: SEX
        warn, error = checkMapping(
            clinicalDF, "SEX", sex_mapping['CODE'], required=True)
        warning += warn
        total_error += error

        # CHECK: ETHNICITY
        warn, error = checkMapping(
            clinicalDF, "ETHNICITY", ethnicity_mapping['CODE'])
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
                    "Patient: All samples must have associated "
                    "patient information"))

        return(clinicalDf)
