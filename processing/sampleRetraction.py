import logging
import os
import pandas as pd
import process_functions
import synapseclient
import datetime
import example_filetype_format
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class sampleRetraction(example_filetype_format.FileTypeFormat):

	_fileType = "sampleRetraction"

	_process_kwargs = ["newPath", "databaseSynId","fileSynId"]

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "%s.csv" % self._fileType

	def process_steps(self, filePath, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		fileSynId = kwargs['fileSynId']
		databaseSynId = kwargs['databaseSynId']
		newPath = kwargs['newPath']
		info = self.syn.get(fileSynId, downloadFile=False)
		retractedSamples = self.syn.tableQuery("select * from %s where center = '%s'" % (databaseSynId,self.center))
		retractedSamplesDf = retractedSamples.asDataFrame()
		deleteSamples = pd.read_csv(filePath,header=None)
		col = 'genieSampleId' if self._fileType == "sampleRetraction" else 'geniePatientId'
		deleteSamples.rename(columns = {0:col}, inplace=True)
		samples = [process_functions.checkGenieId(sample, self.center) for sample in deleteSamples[col]]
		deleteSamples[col] = samples
		modifiedOn = synapseclient.utils.to_unix_epoch_time(datetime.datetime.strptime(info.modifiedOn.split(".")[0], "%Y-%m-%dT%H:%M:%S"))
		deleteSamples['retractionDate'] = modifiedOn
		deleteSamples['center'] = self.center
		process_functions.updateDatabase(self.syn, retractedSamplesDf, deleteSamples, databaseSynId, [col], toDelete=True)
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