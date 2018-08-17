import os
import logging
import process_functions
import clinical
import pandas as pd
import example_filetype_format
import datetime
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class vitalStatus(example_filetype_format.FileTypeFormat):

	_fileType = "vitalStatus"

	## VALIDATING FILENAME
	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "vital_status.txt"
		

	def _validate(self, vitalStatusDf):
		total_error = ""
		warning = ""

		#PATIENT ID
		haveColumn = process_functions.checkColExist(vitalStatusDf, "PATIENT_ID")
		if haveColumn:
			if vitalStatusDf.PATIENT_ID.isnull().any():
				total_error += "Vital status file: Please double check your PATIENT_ID column. No null values allowed.\n"
		else:
			total_error += "Vital status file: Must have PATIENT_ID column.\n"

		#YEAR DEATH
		haveColumn = process_functions.checkColExist(vitalStatusDf, "YEAR_DEATH")
		if haveColumn:
				notNullYears = vitalStatusDf.YEAR_DEATH[~vitalStatusDf.YEAR_DEATH.isnull()]
				try:
					notNullYears.apply(lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
				except:
					total_error += "Vital status file: Please double check your YEAR_DEATH column, it must be an integer in YYYY format or NA/null/empty.\n"
		else:
			total_error += "Vital status file: Must have YEAR_DEATH column.\n"

		#YEAR CONTACT
		haveColumn = process_functions.checkColExist(vitalStatusDf, "YEAR_CONTACT")
		if haveColumn:
			try:
				notNullYears.apply(lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
			except:
				total_error += "Vital status file: Please double check your YEAR_CONTACT column, it must be an integer in YYYY format.\n"
		else:
			total_error += "Vital status file: Must have YEAR_CONTACT column.\n"

		#INT CONTACT
		haveColumn = process_functions.checkColExist(vitalStatusDf, "INT_CONTACT")
		if haveColumn:
			if not all([isinstance(i, int) for i in vitalStatusDf.INT_CONTACT]):
				total_error += "Vital status file: Please double check your INT_CONTACT column, it must be an integer.\n"
		else:
			total_error += "Vital status file: Must have INT_CONTACT column.\n"

		#INT DOD
		haveColumn = process_functions.checkColExist(vitalStatusDf, "INT_DOD")
		if haveColumn:
			if not all([process_functions.checkInt(i) for i in vitalStatusDf.INT_DOD if not pd.isnull(i)]):
				total_error += "Vital status file: Please double check your INT_DOD column, it must be an integer or NA/null/empty.\n"
		else:
			total_error += "Vital status file: Must have INT_DOD column.\n"

		haveColumn = process_functions.checkColExist(vitalStatusDf, "DEAD")
		if haveColumn:
			if not all([isinstance(i, bool) for i in vitalStatusDf.DEAD]):
				total_error += "Vital status file: Please double check your DEAD column, it must be a boolean value.\n"
		else:
			total_error += "Vital status file: Must have DEAD column.\n"

		return(total_error, warning)
	

	def _process(self, vitalStatusDf):
		#vitalStatus_mapping = process_functions.getGenieMapping(self.syn, "syn10888675")

		#noPhiCols = pd.Series(['PATIENT_ID','YEAR_DEATH','YEAR_CONTACT','INT_CONTACT','INT_DOD','DEAD'])

		#vitalStatusDf.VITAL_STATUS = [process_functions.getCODE(vitalStatus_mapping, status) for status in vitalStatusDf.VITAL_STATUS]
		vitalStatusDf.PATIENT_ID = [process_functions.checkGenieId(patient, self.center) for patient in vitalStatusDf.PATIENT_ID]
		vitalStatusDf['CENTER'] = self.center

		return(vitalStatusDf)

	# PROCESS
	def process_steps(self, filePath, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		databaseSynId = kwargs['databaseSynId']
		newPath = kwargs['newPath']
		vitalStatusDf = pd.read_csv(filePath, sep="\t", comment="#")
		vitalStatusDf = self._process(vitalStatusDf)
		cols = vitalStatusDf.columns
		process_functions.updateData(self.syn, databaseSynId, vitalStatusDf, self.center, cols)
		vitalStatusDf.to_csv(newPath, sep="\t",index=False)
		return(newPath)

	## VALIDATION
	# def validate_steps(self, filePathList, **kwargs):
	# 	logger.info("VALIDATING %s" % os.path.basename(filePathList[0]))
	# 	vitalStatusDf = pd.read_csv(filePathList[0], sep="\t")
	# 	total_error, warning = self.validate_helper(vitalStatusDf)
	# 	return(total_error, warning)
