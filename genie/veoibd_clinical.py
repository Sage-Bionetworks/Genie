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

        return(total_error, warning)

    def _get_dataframe(self, filePathList):
        if isinstance(filePathList, Sequence):
            filePathList = filePathList[0]

        df = pd.read_csv(filePathList, comment="#")
        return(df)
