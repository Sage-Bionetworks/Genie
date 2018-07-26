import os
import logging
import process_functions
import clinical
import pandas as pd
import example_filetype_format
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def validate_helper(filePath, vital_status_mapping):
	logger.info("VALIDATING %s" % os.path.basename(filePath))
	total_error = ""
	warning = ""
	noPhiCols = pd.Series(['PATIENT_ID','INT_CONTACT','INT_DOD','VITAL_STATUS'])
	phiCols = pd.Series(['PATIENT_ID','DATE_CONTACT','DOD','VITAL_STATUS'])
	#Null values allowed?
	vitalStatusDf = pd.read_csv(filePath, sep="\t")
	vitalStatusDf = vitalStatusDf.fillna('')
	noPhi = all(noPhiCols.isin(vitalStatusDf.columns))
	phi = all(phiCols.isin(vitalStatusDf.columns))
	if not noPhi and not phi:
		total_error += "Vital status file: Must be a tsv file and either have these columns: (%s) or these columns: (%s)\n" % (", ".join(noPhiCols), ", ".join(phiCols))
	else:
		#CHECK: VITAL_STATUS
		warn, error = clinical.checkMapping(vitalStatusDf,"VITAL_STATUS","VITAL_STATUS",vital_status_mapping['CODE'], required=True)
		warning += warn
		total_error += error
		if noPhi:
			#Remove < and > from all columns
			vitalStatusDf['INT_CONTACT'] = clinical.removeGreaterThanAndLessThanStr(vitalStatusDf['INT_CONTACT'])
			try:
				ints = [int(i) if i != "" else pd.np.nan for i in vitalStatusDf['INT_CONTACT'] ]
			except ValueError:
				total_error += "Vital status file: INT_CONTACT must be an integer or NA/null/empty.\n"                
			vitalStatusDf['INT_DOD'] = clinical.removeGreaterThanAndLessThanStr(vitalStatusDf['INT_DOD'])
			try:
				ints = [int(i) if i != "" else pd.np.nan for i in vitalStatusDf['INT_DOD'] ]
			except ValueError:
				 total_error += "Vital status file: INT_DOD must be an integer or NA/null/empty.\n" 
			   
		else:
			try:
				dates = vitalStatusDf.DATE_CONTACT.apply(lambda date: datetime.datetime.strptime(date, '%m%d%Y'))
			except ValueError as e:
				total_error += "Vital status file: DATE_CONTACT must be in this format mmddyyyy (ie, 03241980).\n"
			try:
				dates = vitalStatusDf.DOD.apply(lambda date: datetime.datetime.strptime(date, '%m%d%Y'))
			except ValueError as e:
				total_error += "Vital status file: DOD must be in this format mmddyyyy (ie, 03241980).\n"

	return(total_error, warning)

class vitalStatus(example_filetype_format.FileTypeFormat):

	_fileType = "vitalStatus"

	## VALIDATING FILENAME
	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "vital_status.txt"
		

	# PROCESS
	def process_steps(self, filePath, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		databaseSynId = kwargs['databaseSynId']
		newPath = kwargs['newPath']
		vitalStatusDf = pd.read_csv(filePath, sep="\t", comment="#")
		vitalStatus_mapping_ent = self.syn.tableQuery('SELECT * FROM syn10888675')
		vitalStatus_mapping = vitalStatus_mapping_ent.asDataFrame()
		#Fill NA's for the ones that are uncoded
		vitalStatus_mapping = vitalStatus_mapping.fillna("")
		vitalStatusDf.VITAL_STATUS = [process_functions.getCODE(vitalStatus_mapping, status) for status in vitalStatusDf.VITAL_STATUS]
		vitalStatusDf['CENTER'] = self.center
		cols = vitalStatusDf.columns
		process_functions.updateData(self.syn, databaseSynId, vitalStatusDf, self.center, cols)
		vitalStatusDf.to_csv(newPath, sep="\t",index=False)
		return(newPath)

	## VALIDATION
	def validate_steps(self, filePathList, **kwargs):
		vital_status_mapping = process_functions.getGenieMapping(self.syn, "syn10888675")
		total_error, warning = validate_helper(filePathList[0], vital_status_mapping)
		return(total_error, warning)

