import os
import logging
import synapseclient
import example_filetype_format
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class workflow(example_filetype_format.FileTypeFormat):

	_fileType = "md"

	_process_kwargs = ["databaseSynId"]

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]).startswith(self.center) and os.path.basename(filePath[0]).endswith(".md")

	def process_steps(self, filePath, *args, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		databaseSynId = kwargs['databaseSynId']
		self.syn.store(synapseclient.File(filePath, parent=databaseSynId))
		return(filePath)

	# def validate_steps(self, filePathList, **kwargs):
	# 	filePath = filePathList[0]
	# 	logger.info("VALIDATING %s" % os.path.basename(filePath))
	# 	total_error = ""
	# 	warning = ""
	# 	logger.info("NO VALIDATION for workflow files")
	# 	return(total_error, warning)