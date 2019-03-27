from __future__ import absolute_import
from genie import FileTypeFormat, process_functions
import logging
import os
import pandas as pd
import synapseclient
import datetime
logger = logging.getLogger(__name__)


class sampleRetraction(FileTypeFormat):

    _fileType = "sampleRetraction"

    _process_kwargs = ["newPath", "databaseSynId","fileSynId"]

    def _get_dataframe(self, filePathList):
        '''
        This function by defaults assumes the filePathList is length of 1 
        and is a tsv file.  Could change depending on file type.
        '''
        filePath = filePathList[0]
        df = pd.read_csv(filePath,header=None)
        return(df)

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "%s.csv" % self._fileType

    def _process(self, deleteSamplesDf, modifiedOn):
        col = 'genieSampleId' if self._fileType == "sampleRetraction" else 'geniePatientId'
        deleteSamplesDf.rename(columns = {0:col}, inplace=True)
        samples = [process_functions.checkGenieId(sample, self.center) for sample in deleteSamplesDf[col]]
        deleteSamplesDf[col] = samples
        modifiedOn = synapseclient.utils.to_unix_epoch_time(datetime.datetime.strptime(modifiedOn, "%Y-%m-%dT%H:%M:%S"))
        deleteSamplesDf['retractionDate'] = modifiedOn
        deleteSamplesDf['center'] = self.center
        return(deleteSamplesDf)

    def process_steps(self, deleteSamples, fileSynId, databaseSynId, newPath):
        info = self.syn.get(fileSynId, downloadFile=False)
        deleteSamples = self._process(deleteSamples, info.modifiedOn.split(".")[0])
        process_functions.updateData(self.syn, databaseSynId, deleteSamples, self.center, filterByColumn="center", toDelete=True)
        return(newPath)
