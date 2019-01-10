from __future__ import absolute_import
from genie import example_filetype_format, process_functions
import logging
import os
import pandas as pd
import synapseclient
import datetime
logger = logging.getLogger(__name__)

class sampleRetraction(example_filetype_format.FileTypeFormat):

	_fileType = "sampleRetraction"

	_process_kwargs = ["newPath", "databaseSynId","fileSynId"]

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

	def process_steps(self, filePath, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		fileSynId = kwargs['fileSynId']
		databaseSynId = kwargs['databaseSynId']
		newPath = kwargs['newPath']
		info = self.syn.get(fileSynId, downloadFile=False)
		deleteSamples = pd.read_csv(filePath,header=None)
		deleteSamples = self._process(deleteSamples, info.modifiedOn.split(".")[0])
		process_functions.updateData(self.syn, databaseSynId, deleteSamples, self.center, filterByColumn="center", toDelete=True)
		return(newPath)

	def validate_steps(self, filePathList, **kwargs):
		filePath = filePathList[0]
		logger.info("VALIDATING %s" % os.path.basename(filePath))
		total_error = ""
		warning = ""
		#send email when there are updates to the GENIE retraction files
		self.syn.sendMessage([3324230,1968150],messageSubject="GENIE: New retraction for %s" % self.center, messageBody="Check GENIE files")
		logger.info("NO VALIDATION for retraction files")
		return(total_error, warning)	