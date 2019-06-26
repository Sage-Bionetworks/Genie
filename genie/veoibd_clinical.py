from collections.abc import Sequence
import datetime
import logging
import os
import tempfile

import pandas as pd
import synapseclient

from .example_filetype_format import FileTypeFormat
from . import process_functions

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Clinical(FileTypeFormat):

    # This should match what is in the database mapping table
    _fileType = "clinical_filetype"

    _required_filename = "clinical_filetype.csv"

    _process_kwargs = [
        "newPath", "parentId", "databaseToSynIdMappingDf"]
    
    # This should match what is in the table that the data goes into
    _required_columns = ["col1", "col2", "center"]

    # This should match what the the primary key is set as an annotation
    # on the table the data goes into.
    _primary_key_columns = ["primary_key_col"]


    def _validateFilename(self, filePath):

        if isinstance(filePath, Sequence):
            filePath = filePath[0]
        
        if os.path.basename(filePath) == self._required_filename:
            logger.debug("{} filename is validated.".format(self._fileType))
        else:
            logger.debug("{} filename is not valid: {}.".format(self._fileType, filePath))
            raise AssertionError("{} filename ({}) is not correct. It should be {}".format(self._fileType,
                                                                                           os.path.basename(filePath),
                                                                                           self._required_filename))


    def _get_dataframe(self, filePathList):
        if isinstance(filePathList, Sequence):
            filePathList = filePathList[0]

        df = pd.read_csv(filePathList, comment="#")
        df['center'] = self.center

        return(df)


    def process_steps(self, data, databaseToSynIdMappingDf, 
                      newPath, parentId):
        patientSynId = databaseToSynIdMappingDf.Id[
            databaseToSynIdMappingDf['Database'] == self._fileType][0]
        
        process_functions.updateData(syn=self.syn, databaseSynId=patientSynId, 
                                     newData=data, filterBy=self.center,
                                     filterByColumn="center", col=self._required_columns,
                                     toDelete=True)
        
        data.to_csv(newPath, sep="\t", index=False)
        return(newPath)


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

        # CHECK: SAMPLE_ID
        _hasColumnDict = dict()
        logger.debug("My required columns are: {}".format(self._required_columns))
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

class ClinicalIndividual(Clinical):

    _fileType = "veoibd_clinical_individual"

    _required_filename = "clinical_individual.csv"
    
    _required_columns = ["individual_id", "age", "sex", "birth_country",
                         "ethnicity","family_hx_ibd","degree_one_with_ibd",
                         "degree_two_with_ibd","initial_dx","gi_site","eim",
                         "dx_perianal","dx_medication","comments", "center"]

    _primary_key_columns = ["individual_id"]


class ClinicalSample(Clinical):

    _fileType = "veoibd_clinical_sample"

    _required_filename = "clinical_sample.csv"
        
    _required_columns = ["sample_id", "individual_id", "assay_id", "center"]

    _primary_key_columns = ["sample_id"]
