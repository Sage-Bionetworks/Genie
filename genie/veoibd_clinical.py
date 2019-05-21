from __future__ import absolute_import
from collections.abc import Sequence

from genie import FileTypeFormat, process_functions
import os
import logging
import pandas as pd
import synapseclient
# import re
import datetime

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

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

    _fileType = "veoibd_clinical"

    _process_kwargs = [
        "newPath", "parentId", "databaseToSynIdMappingDf"]
    
    _required_columns = ["individual_id", "age", "sex", "birth_country",
                         "ethnicity","family_hx_ibd","degree_one_with_ibd",
                         "degree_two_with_ibd","initial_dx","gi_site","eim",
                         "dx_perianal","dx_medication","comments"]

    _primary_key_columns = ["individual_id"]

    # VALIDATE FILE NAME
    def _validateFilename(self, filePath):

        if isinstance(filePath, Sequence):
            filePath = filePath[0]
        
        if os.path.basename(filePath) == "clinical_individual.csv":        
            logger.info("{} filename is validated.".format(self._fileType))
        else:
            logger.error("{} filename is not valid: {}.".format(self._fileType, filePath))
            raise ValueError("Individual clinical filename is not correct.")

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


    def process_steps(self, data, databaseToSynIdMappingDf, 
                      newPath, parentId):
        patientSynId = databaseToSynIdMappingDf.Id[
            databaseToSynIdMappingDf['Database'] == "patient"][0]
        
        patient = True

        data['center'] = self.center

        self.uploadMissingData(data, "individual_id", patientSynId, parentId)
        
        process_functions.updateData(syn=self.syn, databaseSynId=patientSynId, 
                                     newData=data, filterBy=self.center,
                                     filterByColumn="CENTER", col=patientCols,
                                     toDelete=True)
        
        data.to_csv(newPath, sep="\t", index=False)
        return(newPath)

    # VALIDATION
    def _validate(self, data):
        """
        This function validates the clinical file to make sure it adheres
        to the clinical SOP.

        Args:
            data: Pandas data frame with individual metadata

        Returns:
            Error and warning messages
        """
        total_error = ""
        warning = ""

        data = data.fillna("")

        # sampleType_mapping = \
        #     process_functions.getGenieMapping(self.syn, "syn7434273")

        # ethnicity_mapping = \
        #     process_functions.getGenieMapping(self.syn, "syn7434242")

        # race_mapping = \
        #     process_functions.getGenieMapping(self.syn, "syn7434236")

        # sex_mapping = \
        #     process_functions.getGenieMapping(self.syn, "syn7434222")

        # CHECK: SAMPLE_ID
        _hasColumnDict = dict()
        for column in self._required_columns:
            _hasColumnDict[column] = process_functions.checkColExist(data, 
                                                                     column)

            if not _hasColumnDict[column]:
                total_error += \
                    "File must have {} column.\n".format(column)

        # Check if the count of the primary key is not distinct
        res = data.groupby(self._primary_key_columns).size()
        if (res > 1).any():
            total_error += "Found duplicated {primaryKey}'s in the file.".format(primaryKey=self._primary_key_columns)

        # # CHECK: within the sample file that the sample ids match
        # # the patient ids
        # if haveSampleColumn and havePatientColumn:
        #     if not all([patient in sample
        #                 for sample, patient in
        #                 zip(clinicalDF[sampleId], clinicalDF[patientId])]):

        #         total_error += (
        #             "Sample: PATIENT_ID's much be contained in the "
        #             "SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n")
        #     # #CHECK: All samples must have associated patient data
        #     # (GENIE requires patient data)
        #     if not all(clinicalDF[patientId] != ""):
        #         total_error += (
        #             "Patient: All samples must have associated patient "
        #             "information and no null patient ids allowed. These "
        #             "samples are missing patient data: {}\n".format(
        #                 ", ".join(clinicalDF[sampleId][
        #                           clinicalDF[patientId] == ""])))
        #     # CHECK: All patients should have associated sample data
        #     if not all(clinicalDF[sampleId] != ""):
        #         # ## MAKE WARNING FOR NOW###
        #         warning += (
        #             "Sample: All patients must have associated sample "
        #             "information. These patients are missing sample "
        #             "data: {}\n".format(
        #                 ", ".join(clinicalDF[patientId][
        #                           clinicalDF[sampleId] == ""])))

        # # CHECK: AGE_AT_SEQ_REPORT
        # age = "AGE_AT_SEQ_REPORT"
        # haveColumn = process_functions.checkColExist(clinicalDF, age)
        # if haveColumn:
        #     # Deal with HIPAA converted rows from DFCI
        #     # First for loop can't int(text) because there
        #     # are instances that have <3435
        #     age_seq_report_df = \
        #         clinicalDF[~clinicalDF[age].isin(['Unknown', ''])]

        #     age_seq_report_df[age] = \
        #         remove_greaterthan_lessthan_str(age_seq_report_df[age])

        #     if not all([process_functions.checkInt(i)
        #                 for i in age_seq_report_df[age]]):
        #         total_error += (
        #             "Sample: Please double check your AGE_AT_SEQ_REPORT.  "
        #             "It must be an integer or 'Unknown'.\n")
        #     else:
        #         age_seq_report_df[age] = age_seq_report_df[age].astype(int)
        #         median_age = pd.np.median(age_seq_report_df[age])
        #         if median_age < 100:
        #             total_error += (
        #                 "Sample: Please double check your AGE_AT_SEQ_REPORT.  "
        #                 "You may be reporting this value in YEARS, "
        #                 "please report in DAYS.\n")
        # else:
        #     total_error += \
        #         "Sample: clinical file must have AGE_AT_SEQ_REPORT column.\n"

        # # CHECK: SAMPLE_TYPE
        # haveColumn = process_functions.checkColExist(clinicalDF, "SAMPLE_TYPE")
        # if haveColumn:
        #     if clinicalDF.SAMPLE_TYPE.dtype == int:
        #         if not all(clinicalDF['SAMPLE_TYPE'].isin(
        #                 sampleType_mapping['CODE'])):
        #             total_error += (
        #                 "Sample: Please double check your SAMPLE_TYPE column. "
        #                 "This column must be {}.\n".format(
        #                     ", ".join(map(str, sampleType_mapping['CODE']))))
        #     else:
        #         total_error += (
        #             "Sample: Please double check your SAMPLE_TYPE column. "
        #             "No null values allowed.\n")

        # else:
        #     total_error += \
        #         "Sample: clinical file must have SAMPLE_TYPE column.\n"

        # # CHECK: SEQ_ASSAY_ID
        # haveColumn = \
        #     process_functions.checkColExist(clinicalDF, "SEQ_ASSAY_ID")
        # if haveColumn:
        #     if not all([i != "" for i in clinicalDF['SEQ_ASSAY_ID']]):
        #         total_error += (
        #             "Sample: Please double check your SEQ_ASSAY_ID columns, "
        #             "there are empty rows.\n")
        #     # must remove empty seq assay ids first
        #     # Checking if seq assay ids start with the center name
        #     seqAssayIds = \
        #         clinicalDF.SEQ_ASSAY_ID[clinicalDF.SEQ_ASSAY_ID != ""]
        #     allSeqAssays = seqAssayIds.unique()
        #     notNormalized = []
        #     not_caps = []
        #     for seqassay in allSeqAssays:
        #         # SEQ Ids are all capitalized now, so no need to check
        #         # for differences in case
        #         if not seqassay.upper().startswith(self.center):
        #             not_caps.append(seqassay)
        #     if len(not_caps) > 0:
        #         total_error += (
        #             "Sample: Please make sure your SEQ_ASSAY_IDs start with "
        #             "your center abbreviation: {}.\n".format(
        #                 ", ".join(not_caps)))
        # else:
        #     total_error += \
        #         "Sample: clinical file must have SEQ_ASSAY_ID column.\n"

        # haveColumn = process_functions.checkColExist(clinicalDF, "SEQ_DATE")
        # seq_date_error = (
        #     "Sample: SEQ_DATE must be one of five values- "
        #     "For Jan-March: use Jan-YEAR. "
        #     "For Apr-June: use Apr-YEAR. "
        #     "For July-Sep: use Jul-YEAR. "
        #     "For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) "
        #     "For values that don't have SEQ_DATES that "
        #     "you want released use 'release'.\n")

        # if haveColumn:
        #     clinicalDF['SEQ_DATE'] = [
        #         i.title() for i in clinicalDF['SEQ_DATE'].astype(str)]

        #     seqDate = clinicalDF['SEQ_DATE'][
        #         clinicalDF['SEQ_DATE'] != 'Release']
        #     if sum(clinicalDF['SEQ_DATE'] == '') > 0:
        #         total_error += \
        #             "Sample: Samples without SEQ_DATEs will NOT be released.\n"
        #     try:
        #         if not seqDate.empty:
        #             dates = seqDate.apply(
        #                 lambda date: datetime.datetime.strptime(date, '%b-%Y'))
        #             # REMOVE JUN LATER
        #             if not all([i.startswith(("Jul", "Jan", "Oct", "Apr"))
        #                         for i in seqDate]):
        #                 total_error += seq_date_error
        #     except ValueError:
        #         total_error += seq_date_error
        # else:
        #     total_error += "Sample: clinical file must SEQ_DATE column\n"

        # # CHECK: BIRTH_YEAR
        # birth_year = "BIRTH_YEAR"
        # haveColumn = process_functions.checkColExist(clinicalDF, birth_year)
        # if haveColumn:
        #     # Deal with HIPAA converted rows from DFCI
        #     # First for loop can't int(text) because there are
        #     # instances that have <YYYY
        #     # Remove '' for blank value support
        #     birth_year_df = \
        #         clinicalDF[~clinicalDF[birth_year].isin(['Unknown', ''])]
        #     birth_year_df[birth_year] = \
        #         remove_greaterthan_lessthan_str(birth_year_df[birth_year])

        #     try:
        #         years = birth_year_df[birth_year].apply(
        #             lambda x: datetime.datetime.strptime(
        #                 str(int(x)), '%Y').year >
        #             datetime.datetime.utcnow().year)

        #         assert not years.any()
        #     except Exception:
        #         total_error += (
        #             "Patient: Please double check your BIRTH_YEAR column, "
        #             "it must be an integer in YYYY format > {year} or "
        #             "'Unknown'.  Support for blank values will be deprecated "
        #             "in 7...releases.\n".format(
        #                 year=datetime.datetime.utcnow().year))
        # else:
        #     total_error += \
        #         "Patient: clinical file must have BIRTH_YEAR column.\n"

        # # CHECK: VITAL_STATUS
        # # YEAR DEATH
        # haveColumn = process_functions.checkColExist(clinicalDF, "YEAR_DEATH")
        # if haveColumn:
        #     notNullYears = clinicalDF.YEAR_DEATH[~clinicalDF.YEAR_DEATH.isin(
        #         ['Unknown', 'Not Collected', 'Not Applicable'])]
        #     try:
        #         notNullYears.apply(
        #             lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
        #     except Exception:
        #         total_error += (
        #             "Patient: Please double check your YEAR_DEATH column, "
        #             "it must be an integer in YYYY format, "
        #             "'Unknown', 'Not Applicable' or 'Not Collected'.\n")
        # else:
        #     warning += (
        #         "Patient: Must have YEAR_DEATH column for "
        #         "7...release uploads.\n")

        # # YEAR CONTACT
        # haveColumn = process_functions.checkColExist(
        #     clinicalDF, "YEAR_CONTACT")
        # if haveColumn:
        #     notNullYears = clinicalDF.YEAR_CONTACT[
        #         ~clinicalDF.YEAR_CONTACT.isin(['Unknown', 'Not Collected'])]
        #     try:
        #         notNullYears.apply(
        #             lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
        #     except Exception:
        #         total_error += (
        #             "Patient: Please double check your YEAR_CONTACT column, "
        #             "it must be an integer in YYYY format, "
        #             "'Unknown' or 'Not Collected'.\n")
        # else:
        #     warning += (
        #         "Patient: Must have YEAR_CONTACT column for "
        #         "7...release uploads.\n")

        # # INT CONTACT
        # haveColumn = process_functions.checkColExist(clinicalDF, "INT_CONTACT")
        # if haveColumn:
        #     if not all([
        #             process_functions.checkInt(i)
        #             for i in clinicalDF.INT_CONTACT if i not in
        #             ['>32485', '<6570', 'Unknown', 'Not Collected']]):

        #         total_error += (
        #             "Patient: Please double check your INT_CONTACT column, "
        #             "it must be an integer, '>32485', '<6570', 'Unknown' "
        #             "or 'Not Collected'.\n")
        # else:
        #     warning += (
        #         "Patient: Must have INT_CONTACT column for "
        #         "7...release uploads.\n")

        # # INT DOD
        # haveColumn = process_functions.checkColExist(clinicalDF, "INT_DOD")
        # if haveColumn:
        #     if not all([
        #             process_functions.checkInt(i)
        #             for i in clinicalDF.INT_DOD if i not in
        #             ['>32485', '<6570', 'Unknown',
        #              'Not Collected', 'Not Applicable']]):

        #         total_error += (
        #             "Patient: Please double check your INT_DOD column, "
        #             "it must be an integer, '>32485', '<6570', 'Unknown', "
        #             "'Not Collected' or 'Not Applicable'.\n")
        # else:
        #     warning += \
        #         "Patient: Must have INT_DOD column for 7...release uploads.\n"

        # haveColumn = process_functions.checkColExist(clinicalDF, "DEAD")
        # if haveColumn:
        #     # Need to have check_bool function
        #     if not all([
        #             str(i).upper() in ['TRUE', 'FALSE']
        #             for i in clinicalDF.DEAD if i not in
        #             ['Unknown', 'Not Collected']]):
        #         total_error += (
        #             "Patient: Please double check your DEAD column, "
        #             "it must be True, False, 'Unknown' or 'Not Collected'.\n")
        # else:
        #     warning += \
        #         "Patient: Must have DEAD column for 7...release uploads.\n"

        # # CHECK: PRIMARY_RACE
        # warn, error = checkMapping(
        #     clinicalDF, "PRIMARY_RACE", race_mapping['CODE'])
        # warning += warn
        # total_error += error

        # # CHECK: SECONDARY_RACE
        # warn, error = checkMapping(
        #     clinicalDF, "SECONDARY_RACE", race_mapping['CODE'])
        # warning += warn
        # total_error += error

        # # CHECK: TERTIARY_RACE
        # warn, error = checkMapping(
        #     clinicalDF, "TERTIARY_RACE", race_mapping['CODE'])
        # warning += warn
        # total_error += error

        # # CHECK: SEX
        # warn, error = checkMapping(
        #     clinicalDF, "SEX", sex_mapping['CODE'], required=True)
        # warning += warn
        # total_error += error

        # # CHECK: ETHNICITY
        # warn, error = checkMapping(
        #     clinicalDF, "ETHNICITY", ethnicity_mapping['CODE'])
        # warning += warn
        # total_error += error

        return(total_error, warning)

    def _get_dataframe(self, filePathList):
        if isinstance(filePathList, Sequence):
            filePathList = filePathList[0]

        df = pd.read_csv(filePathList, comment="#")
        return(df)
