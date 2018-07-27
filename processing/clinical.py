import process_functions
import os
import logging
import pandas as pd
import synapseclient
import example_filetype_format
import re
import datetime
#import time
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#CHECKS IF THE MAPPING IS CORRECT
def checkMapping(clinicalDF, colName, mapping, required=False, fileType = "Patient"):
	"""
	This function checks if the column exists then checks if the values in the column have the correct integer values
	
	:params clinicalDF          Patient/sample/flattened clinical file
	:params colName:        	Expected column name
	:params mapping:            List of possible values

	:returns:                   A tuple warning, error
	"""
	warning = ""
	error = ""
	haveColumn = process_functions.checkColExist(clinicalDF, colName)
	if not haveColumn:
		if required:
			error = "%s: clinical file must have %s column.\n" % (fileType, colName)
		else:
			warning = "%s: clinical file doesn't have %s column. A blank column will be added\n" % (fileType, colName)
	else:
		if not all([i in mapping.tolist() for i in clinicalDF[colName]]):
			error = "%s: Please double check your %s column.  This column must be these values %sor blank.\n" % (fileType, colName, ", ".join(map(str,mapping)).replace(".0",""))
	return(warning, error)

#In clinical file, there are redacted value such as >89 and <17.  These < and > signs must be removed
def removeGreaterThanAndLessThanStr(col):
	try:
		col = [text.replace(">","") if isinstance(text, str) else text for text in col]
		col = [int(text.replace("<","")) if isinstance(text, str) and text != "" else text for text in col]
	except ValueError:
		pass
	return(col)

class clinical(example_filetype_format.FileTypeFormat):

	_fileType = "clinical"

	#_process_kwargs = ["newPath", "patientSynId", "sampleSynId","parentId","retractedSampleSynId","retractedPatientSynId"]

	_process_kwargs = ["newPath", "patientSynId", "sampleSynId","parentId","oncotreeLink"]

	_validation_kwargs = ["oncotreeLink"]

	# VALIDATE FILE NAME
	def _validateFilename(self, filePath):
		if len(filePath) == 1:
			assert os.path.basename(filePath[0]) == "data_clinical_supp_%s.txt" % self.center
		else:
			required = pd.Series(["data_clinical_supp_sample_%s.txt" % self.center,"data_clinical_supp_patient_%s.txt" % self.center])
			assert all(required.isin([os.path.basename(i) for i in filePath]))

	# PROCESSING
	#Update clinical file with the correct mappings
	def update_clinical(self, x, sex_mapping,race_mapping,ethnicity_mapping,sample_type): 
		#PATIENT ID
		if x.get("PATIENT_ID") is not None:
			x['PATIENT_ID'] = process_functions.checkGenieId(x['PATIENT_ID'], self.center)
		# RACE
		if x.get('PRIMARY_RACE') is not None:
			x['PRIMARY_RACE'] = process_functions.getCODE(race_mapping, x['PRIMARY_RACE'])
		if x.get('SECONDARY_RACE') is not None:
			x['SECONDARY_RACE'] = process_functions.getCODE(race_mapping, x['SECONDARY_RACE'])
		if x.get('TERTIARY_RACE') is not None:
			x['TERTIARY_RACE'] = process_functions.getCODE(race_mapping, x['TERTIARY_RACE'])
		# ETHNICITY
		if x.get('ETHNICITY') is not None:
			x['ETHNICITY'] = process_functions.getCODE(ethnicity_mapping, x['ETHNICITY'])
		# BIRTH YEAR
		if x.get("BIRTH_YEAR") is not None:
			# BIRTH YEAR (Check if integer)
			if process_functions.checkInt(x['BIRTH_YEAR']):
				x['BIRTH_YEAR'] = int(x['BIRTH_YEAR'])
		# SEX
		if x.get("SEX") is not None:
			x['SEX'] = process_functions.getCODE(sex_mapping, x['SEX'])
		#TRIM EVERY COLUMN MAKE ALL DASHES 
		#SAMPLE ID
		if x.get('SAMPLE_ID') is not None:
			x['SAMPLE_ID'] = process_functions.checkGenieId(x['SAMPLE_ID'], self.center)
		#AGE AT SEQ REPORT
		if x.get('AGE_AT_SEQ_REPORT') is not None:
			if process_functions.checkInt(x['AGE_AT_SEQ_REPORT']):
				x['AGE_AT_SEQ_REPORT'] = int(x['AGE_AT_SEQ_REPORT'])

		#SEQ ASSAY ID
		if x.get('SEQ_ASSAY_ID') is not None:
			if not str(x['SEQ_ASSAY_ID']).startswith(self.center) and str(x['SEQ_ASSAY_ID']) != "":
				x['SEQ_ASSAY_ID'] = "%s-%s" % (self.center, str(x['SEQ_ASSAY_ID']))
			x['SEQ_ASSAY_ID'] = x['SEQ_ASSAY_ID'].replace('_','-')
			#standardize all SEQ_ASSAY_ID with uppercase
			x['SEQ_ASSAY_ID'] = x['SEQ_ASSAY_ID'].upper()
		
		#SAMPLE_TYPE
		if x.get('SAMPLE_TYPE') is not None:
			sampleType = x['SAMPLE_TYPE']
			x['SAMPLE_TYPE'] = process_functions.getCODE(sample_type, sampleType)
			#Trim spaces
			x['SAMPLE_TYPE_DETAILED'] = process_functions.getCODE(sample_type, sampleType, useDescription=True)

		if x.get('SEQ_DATE') is not None:
			x['SEQ_DATE'] = x['SEQ_DATE'].title()
			x['SEQ_YEAR'] = int(str(x['SEQ_DATE']).split("-")[1]) if str(x['SEQ_DATE']) != "Release" else pd.np.nan

		#TRIM EVERY COLUMN MAKE ALL DASHES 
		for i in x.keys():
			if isinstance(x[i],str):
				x[i] = x[i].strip(" ")
		return(x)

	def uploadMissingData(self, df, col, dbSynId, stagingSynId, retractionSynId=None):
		samples = "','".join(df[col])
		path = os.path.join(process_functions.SCRIPT_DIR, "%s_missing_%s.csv" % (self._fileType, col))
		missing = self.syn.tableQuery("select %s from %s where CENTER = '%s' and %s not in ('%s')" % (col, dbSynId, self.center, col, samples))
		missing.asDataFrame().to_csv(path, index=False)
		self.syn.store(synapseclient.File(path, parent=stagingSynId))
		os.remove(path)

	def preprocess(self, filePath, **kwargs):
		databaseToSynIdMappingDf = kwargs['databaseToSynIdMappingDf']
		patientSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "patient"][0]
		sampleSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "sample"][0]
		# retractedSampleSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "sampleRetraction"][0]
		# retractedPatientSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "patientRetraction"][0]
		#path = process_helper(syn, filePath, center, newPath, patientSynId, sampleSynId, centerStagingSynId, fileType)
		return({"patientSynId":patientSynId,"sampleSynId":sampleSynId})#"retractedSampleSynId":retractedSampleSynId,"retractedPatientSynId":retractedPatientSynId})


	def _process(self, clinical, clinicalTemplate):
		clinicalMerged = clinical.merge(clinicalTemplate,how='outer')
		clinicalMerged = clinicalMerged.drop(clinicalMerged.columns[~clinicalMerged.columns.isin(clinicalTemplate.columns)],1)

		ethnicity_mapping =process_functions.getGenieMapping(self.syn, "syn7434242")
		race_mapping = process_functions.getGenieMapping(self.syn, "syn7434236")
		sex_mapping = process_functions.getGenieMapping(self.syn, "syn7434222")
		sampleType_mapping = process_functions.getGenieMapping(self.syn, "syn7434273")
		sampleType_mapping = sampleType_mapping.append(pd.DataFrame([("","","")],columns=["CODE","CBIO_LABEL","DESCRIPTION"]))
		#Attach MSK to centers
		clinicalMerged = clinicalMerged.fillna("")
		clinicalRemapped = clinicalMerged.apply(lambda x: self.update_clinical(x, sex_mapping, race_mapping, ethnicity_mapping, sampleType_mapping),1)	
		clinicalRemapped['CENTER'] = self.center
		return(clinicalRemapped)

	def process_steps(self, filePath, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		patientSynId = kwargs['patientSynId']
		sampleSynId = kwargs['patientSynId']
		newPath = kwargs['newPath']
		patientSynId = kwargs['patientSynId']
		sampleSynId = kwargs['sampleSynId']
		centerStagingSynId = kwargs['parentId']
		oncotreeLink = kwargs['oncotreeLink']
		# retractedSampleSynId = kwargs['retractedSampleSynId']
		# retractedPatientSynId = kwargs['retractedPatientSynId']

		clinicalDf = pd.read_csv(filePath, sep="\t", comment="#")
		patient= False
		sample = False
		#These synapse ids for the clinical tier release scope is hardcoded because it never changes
		patientColsTable = self.syn.tableQuery('select fieldName from syn8545211 where patient is True')
		patientCols = patientColsTable.asDataFrame()['fieldName'].tolist()
		sampleColsTable = self.syn.tableQuery('select fieldName from syn8545211 where sample is True')
		sampleCols = sampleColsTable.asDataFrame()['fieldName'].tolist()

		if "patient" in filePath.lower():
			clinicalTemplate = pd.DataFrame(columns=patientCols)
			patient = True
		elif "sample" in filePath.lower():
			clinicalTemplate = pd.DataFrame(columns=sampleCols)
			sample = True
		else:
			clinicalTemplate = pd.DataFrame(columns=set(patientCols + sampleCols))
			sample = True 
			patient = True

		newClinicalDf = self._process(clinicalDf, clinicalTemplate)

		if patient:
			seqColumn = "BIRTH_YEAR"
			patientCols.append(seqColumn + "_NUMERICAL")
			newClinicalDf[seqColumn + "_NUMERICAL"] = [int(year) if process_functions.checkInt(year) else pd.np.nan for year in newClinicalDf[seqColumn]]
			patientClinical = newClinicalDf[patientCols].drop_duplicates("PATIENT_ID")
			self.uploadMissingData(patientClinical, "PATIENT_ID", patientSynId, centerStagingSynId)#retractedPatientSynId)
			process_functions.updateData(self.syn, patientSynId, patientClinical, self.center, patientCols, toDelete=True)
		if sample:
			seqColumn = "AGE_AT_SEQ_REPORT"
			sampleCols.extend([seqColumn + "_NUMERICAL"])
			newClinicalDf[seqColumn + "_NUMERICAL"] = [int(year) if process_functions.checkInt(year) else pd.np.nan for year in newClinicalDf[seqColumn]]
			if sum(newClinicalDf["SAMPLE_ID"].duplicated()) >0:
				logger.error("There are duplicated samples, and the duplicates are removed")
			sampleClinical = newClinicalDf[sampleCols].drop_duplicates("SAMPLE_ID")
			#Exclude all clinical samples with wrong oncotree codes
			oncotree_mapping = process_functions.get_oncotree_codes(oncotreeLink)
			if oncotree_mapping.empty:
				oncotree_mapping_dict = process_functions.get_oncotree_code_mappings(oncotreeLink)
				oncotree_mapping['ONCOTREE_CODE'] = oncotree_mapping_dict.keys()
			#Make oncotree codes uppercase (SpCC/SPCC)
			sampleClinical['ONCOTREE_CODE'] = sampleClinical['ONCOTREE_CODE'].astype(str).str.upper()
			sampleClinical = sampleClinical[sampleClinical['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])]
			self.uploadMissingData(sampleClinical, "SAMPLE_ID", sampleSynId, centerStagingSynId)#, retractedSampleSynId)
			process_functions.updateData(self.syn, sampleSynId, sampleClinical, self.center, sampleCols, toDelete=True)

		newClinicalDf.to_csv(newPath, sep="\t", index=False)
		return(newPath)

	# VALIDATION
	def validate_helper(self, clinicalDF, clinicalSampleDF, oncotreeLink):
		"""
		This function validates the clinical file to make sure it adhere to the clinical SOP.
		
		:params clinicalFilePath:              Flattened clinical file or patient clinical file
		:params clinicalSamplePath:            Sample clinical file if patient clinical file is provided

		:returns:                              Error message
		"""
		total_error = ""
		warning = ""

		clinicalDF.columns = [col.upper() for col in clinicalDF.columns]
		clinicalDF = clinicalDF.fillna("")

		clinicalSampleDF.columns = [col.upper() for col in clinicalSampleDF.columns]
		clinicalSampleDF = clinicalSampleDF.fillna("")
		oncotree_mapping = process_functions.get_oncotree_codes(oncotreeLink)
		if oncotree_mapping.empty:
			oncotree_mapping_dict = process_functions.get_oncotree_code_mappings(oncotreeLink)
			oncotree_mapping['ONCOTREE_CODE'] = oncotree_mapping_dict.keys()
		sampleType_mapping = process_functions.getGenieMapping(self.syn, "syn7434273")
		ethnicity_mapping =process_functions.getGenieMapping(self.syn, "syn7434242")
		race_mapping = process_functions.getGenieMapping(self.syn, "syn7434236")
		sex_mapping = process_functions.getGenieMapping(self.syn, "syn7434222")

		#CHECK: SAMPLE_ID
		sampleId = 'SAMPLE_ID'
		haveSampleColumn = process_functions.checkColExist(clinicalSampleDF, sampleId)
		if not haveSampleColumn:
			total_error += "Sample: clinical file must have SAMPLE_ID column.\n"
		else:
			if sum(clinicalSampleDF[sampleId] == "") > 0:
				total_error += "Sample: There can't be any blank values for SAMPLE_ID\n"
			if sum(clinicalSampleDF[sampleId].duplicated()) > 0:
				total_error += "Sample: There can't be any duplicated values for SAMPLE_ID\n"

		#CHECK: AGE_AT_SEQ_REPORT
		age = "AGE_AT_SEQ_REPORT"
		haveColumn = process_functions.checkColExist(clinicalSampleDF, age)
		if haveColumn:
			#Deal with HIPAA converted rows from DFCI
			#First for loop can't int(text) because there are instances that have <3435
			clinicalSampleDF[age] = removeGreaterThanAndLessThanStr(clinicalSampleDF[age]) 
			if not all([isinstance(i, (int,float)) or i == "" for i in clinicalSampleDF[age]]) or pd.np.median(clinicalSampleDF[age]) < 100:
				total_error += "Sample: Please double check your AGE_AT_SEQ_REPORT.  This is the interval in DAYS (integer) between the patient's date of birth and the date of the sequencing report that is associated with the sample.\n"
		else:
			total_error += "Sample: clinical file must have AGE_AT_SEQ_REPORT column.\n"

		#CHECK: ONCOTREE_CODE
		haveColumn = process_functions.checkColExist(clinicalSampleDF, "ONCOTREE_CODE")
		maleOncoCodes = ["TESTIS","PROSTATE","PENIS"]
		womenOncoCodes = ["CERVIX","VULVA","UTERUS","OVARY"]
		if haveColumn:
			#Make oncotree codes uppercase (SpCC/SPCC)
			clinicalSampleDF['ONCOTREE_CODE'] = clinicalSampleDF['ONCOTREE_CODE'].astype(str).str.upper()
			if not all(clinicalSampleDF['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])):
				unmapped_oncotrees = clinicalSampleDF['ONCOTREE_CODE'][~clinicalSampleDF['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])]
				total_error += "Sample: Please double check that all your ONCOTREE CODES exist in the mapping. You have %d samples that don't map. These are the codes that don't map: %s\n" % (len(unmapped_oncotrees),",".join(set(unmapped_oncotrees)))
			if process_functions.checkColExist(clinicalDF, "SEX") and 'oncotree_mapping_dict' in locals():
				wrongCodeSamples = []
				#This is to check if oncotree codes match the sex, returns list of samples that have conflicting codes and sex
				for code, patient, sample in zip(clinicalSampleDF['ONCOTREE_CODE'], clinicalSampleDF['PATIENT_ID'], clinicalSampleDF['SAMPLE_ID']):
					if oncotree_mapping_dict.get(code) is not None and sum(clinicalDF['PATIENT_ID'] == patient) > 0:
						primaryCode = oncotree_mapping_dict[code]['ONCOTREE_PRIMARY_NODE']
						sex = clinicalDF['SEX'][clinicalDF['PATIENT_ID'] == patient].values[0]
						if oncotree_mapping_dict[code]['ONCOTREE_PRIMARY_NODE'] in maleOncoCodes and sex != "Male":
							wrongCodeSamples.append(sample)
						if oncotree_mapping_dict[code]['ONCOTREE_PRIMARY_NODE'] in womenOncoCodes and sex != "Female":
							wrongCodeSamples.append(sample)
				if len(wrongCodeSamples) > 0:
					warning += "Sample: Some SAMPLE_IDs have conflicting SEX and ONCOTREE_CODES: %s\n" % ",".join(wrongCodeSamples)
		else:
			total_error += "Sample: clinical file must have ONCOTREE_CODE column.\n"
		
		#CHECK: SAMPLE_TYPE
		haveColumn = process_functions.checkColExist(clinicalSampleDF, "SAMPLE_TYPE")
		if haveColumn:
			if clinicalSampleDF.SAMPLE_TYPE.dtype == int:
				if not all(clinicalSampleDF['SAMPLE_TYPE'].isin(sampleType_mapping['CODE'])):
					total_error += "Sample: Please double check your SAMPLE_TYPE column. This column must be %s.\n" % ", ".join(map(str,sampleType_mapping['CODE']))
			else:
				total_error += "Sample: Please double check your SAMPLE_TYPE column. No null values allowed.\n"

		else:
			total_error += "Sample: clinical file must have SAMPLE_TYPE column.\n"

		#CHECK: SEQ_ASSAY_ID
		haveColumn = process_functions.checkColExist(clinicalSampleDF, "SEQ_ASSAY_ID")
		if haveColumn:
			if not all([i != "" for i in clinicalSampleDF['SEQ_ASSAY_ID']]):
				warning += "Sample: Please double check your SEQ_ASSAY_ID columns, there are empty rows.\n"
			#must remove empty seq assay ids first or else the case checking of values will fail
			#Checking for if there are seq assay ids like (SAGE-1, Sage-1)
			seqAssayIds = clinicalSampleDF.SEQ_ASSAY_ID[clinicalSampleDF.SEQ_ASSAY_ID != ""]
			allSeqAssays = seqAssayIds.unique()
			notNormalized = []
			for seqassay in allSeqAssays:
				checkIfNormalized = [bool(re.search(seqassay,seq,re.IGNORECASE)) for seq in allSeqAssays]
				if sum(checkIfNormalized) > 1:
					notNormalized.append(seqassay)
			if len(notNormalized) > 0:
				total_error += "Sample: Please normalize your SEQ_ASSAY_ID names.  You have these SEQ_ASSAY_IDs: %s.\n" % ", ".join(notNormalized)
		else:
			total_error += "Sample: clinical file must have SEQ_ASSAY_ID column.\n"

		haveColumn = process_functions.checkColExist(clinicalSampleDF, "SEQ_DATE")
		if haveColumn:
			clinicalSampleDF['SEQ_DATE'] = [i.title() for i in clinicalSampleDF['SEQ_DATE'].astype(str)]
			seqDate = clinicalSampleDF['SEQ_DATE'][clinicalSampleDF['SEQ_DATE'] != 'Release']
			if sum(clinicalSampleDF['SEQ_DATE'] == '') >0:
				total_error += "Sample: Samples without SEQ_DATEs will NOT be released.\n"
			try:
				if not seqDate.empty:
					dates = seqDate.apply(lambda date: datetime.datetime.strptime(date, '%b-%Y'))
					#REMOVE JUN LATER
					if not all([i.startswith(("Jul","Jan","Oct","Apr")) for i in seqDate]):
						total_error += "Sample: SEQ_DATE must be one of five values- For Jan-March: use Jan-YEAR. For Apr-June: use Apr-YEAR. For July-Sep: use Jul-YEAR. For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) For values that don't have SEQ_DATES that you want released use 'release'.\n"
			except ValueError as e:
				total_error += "Sample: SEQ_DATE must be one of five values- For Jan-March: use Jan-YEAR. For Apr-June: use Apr-YEAR. For July-Sep: use Jul-YEAR. For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) For values that don't have SEQ_DATES that you want released use 'release'.\n"
		else:
			total_error += "Sample: clinical file must SEQ_DATE column\n"
			#warning += "Sample: clinical file does not have SEQ_DATE column, ALL samples will be released\n"

		#CHECK: BIRTH_YEAR
		birth_year = "BIRTH_YEAR"
		haveColumn = process_functions.checkColExist(clinicalDF, birth_year)
		if haveColumn: 
			#Deal with HIPAA converted rows from DFCI
			#First for loop can't int(text) because there are instances that have <3435 
			clinicalDF[birth_year] = removeGreaterThanAndLessThanStr(clinicalDF[birth_year])
			if not all([isinstance(i, (int,float)) or i == "" for i in clinicalDF[birth_year]]):
				total_error += "Patient: Please double check your BIRTH_YEAR column.  This column must be integers or blank.\n"
		else:
			total_error += "Patient: clinical file must have BIRTH_YEAR column.\n"


		#CHECK: PATIENT_ID
		patientId = "PATIENT_ID"
		haveColumn = process_functions.checkColExist(clinicalDF, patientId)
		if not haveColumn:
			total_error += "Patient: clinical file must have PATIENT_ID column.\n"
		else:
			if sum(clinicalDF[patientId] == "") > 0:
				total_error += "Patient: There can't be any blank values for PATIENT_ID\n"

		#CHECK: PATIENT_ID IN SAMPLE FILE
		haveColumn = process_functions.checkColExist(clinicalSampleDF, patientId)
		if not haveColumn:
			total_error += "Sample: clinical file must have PATIENT_ID column.\n"
		else:
			if sum(clinicalSampleDF[patientId] == "") > 0:
				total_error += "Sample: There can't be any blank values for PATIENT_ID\n"

			#CHECK: within the sample file that the sample ids match the patient ids
			if haveSampleColumn:
				if not all([patient in sample for sample, patient in zip(clinicalSampleDF[sampleId], clinicalSampleDF[patientId])]):
					total_error += "Sample: PATIENT_ID's much be contained in the SAMPLE_ID's (ex. SAGE-1 <-> SAGE-1-2)\n"
			#Remove empty patient values from both files
			sample_patients = clinicalSampleDF[patientId][clinicalSampleDF[patientId] != ""]
			patient_patients = clinicalDF[patientId][clinicalDF[patientId] != ""]
			# #CHECK: All samples must have associated patient data (GENIE requires patient data)
			if not all(sample_patients.isin(patient_patients)):
				total_error += "Sample: All samples must have associated patient information. These samples are missing patient data: %s\n" % ", ".join(clinicalSampleDF[patientId][~clinicalSampleDF[patientId].isin(clinicalDF[patientId])])
			#CHECK: All patients must have associated sample data 
			if not all(patient_patients.isin(sample_patients)):
				### MAKE WARNING FOR NOW###
				warning += "Sample: All patients must have associated sample information. These patients are missing sample data: %s\n" % ", ".join(clinicalDF[patientId][~clinicalDF[patientId].isin(clinicalSampleDF[patientId])])

		#CHECK: PRIMARY_RACE
		warn, error = checkMapping(clinicalDF,"PRIMARY_RACE",race_mapping['CODE'])
		warning += warn
		total_error  += error 

		#CHECK: SECONDARY_RACE
		warn, error = checkMapping(clinicalDF,"SECONDARY_RACE",race_mapping['CODE'])
		warning += warn
		total_error  += error 

		#CHECK: TERTIARY_RACE
		warn, error = checkMapping(clinicalDF,"TERTIARY_RACE",race_mapping['CODE'])
		warning += warn
		total_error  += error 

		#CHECK: SEX
		warn, error = checkMapping(clinicalDF,"SEX",sex_mapping['CODE'], required=True)
		warning += warn
		total_error  += error 

		#CHECK: ETHNICITY
		warn, error = checkMapping(clinicalDF,"ETHNICITY",ethnicity_mapping['CODE'])
		warning += warn
		total_error  += error

		return(total_error, warning)

	def validate_steps(self, filePathList, **kwargs):

		logger.info("VALIDATING %s" % ", ".join([ os.path.basename(i) for i in filePathList]))
		oncotreeLink = kwargs['oncotreeLink']
		clinicalDF = pd.read_csv(filePathList[0],sep="\t",comment="#")

		if len(filePathList) > 1:
			otherClinicalDF = pd.read_csv(filePathList[1],sep="\t",comment="#")
			if "patient" in filePathList[0].lower():
				total_error, warning = self.validate_helper(clinicalDF, otherClinicalDF, oncotreeLink)
			else:
				total_error, warning = self.validate_helper(otherClinicalDF, clinicalDF, oncotreeLink)
		else:
			otherClinicalDF = clinicalDF.copy()
			total_error, warning = self.validate_helper(clinicalDF, otherClinicalDF, oncotreeLink)
		return(total_error, warning)