import os
import logging
import process_functions
import subprocess
import pandas as pd
import example_filetype_format
import synapseclient

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class maf(example_filetype_format.FileTypeFormat):

	_fileType = "maf"

	_process_kwargs = ["databaseToSynIdMappingDf","processing","path_to_GENIE",
					   "vcf2mafPath","veppath","vepdata",'validMAFs']

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "data_mutations_extended_%s.txt" % self.center

	#Format and storing maf file
	def formatMAF(self, mafDf):
		mafDf['Center'] = self.center
		mafDf['Tumor_Sample_Barcode'] = [process_functions.checkGenieId(i,self.center) for i in mafDf['Tumor_Sample_Barcode']]
		mafDf['Sequence_Source'] = pd.np.nan
		mafDf['Sequencer'] = pd.np.nan
		mafDf['Validation_Status'][mafDf['Validation_Status'].isin(["Unknown","unknown"])] = ''
		return(mafDf)

	def createFinalMaf(self, mafDf, filePath, maf=False):
		if not mafDf.empty:
			if os.stat(filePath).st_size == 0 or maf:
				mafSet = mafDf.to_csv(sep="\t",index=False)
			else:
				mafSet = mafDf.to_csv(sep="\t", index=False, header=None)
			writeOrAppend = "w" if maf else "a"
			with open(filePath, writeOrAppend) as maf:
				maf.write(mafSet)

	#There is a isNarrow option, but note that the number of rows of the maf file 
	#DOES NOT change in this function
	def storeProcessedMaf(self, filePath, mafSynId, centerMafSynId, isNarrow=False):
		logger.info('STORING %s' % filePath)
		database = self.syn.get(mafSynId)
		if isNarrow:
			self.syn.store(synapseclient.Table(database.id, filePath, separator="\t"))
		else:
			self.syn.store(synapseclient.File(filePath, parentId=centerMafSynId))
		return(filePath)
		
	def process_helper(self, filePath, path_to_GENIE, mafSynId, centerMafSynId,
					   vcf2mafPath, veppath, vepdata):
		logger.info('MAF2MAF %s' % filePath)
		fileName = "data_mutations_extended_%s_MAF.txt" % self.center
		newMafPath = os.path.join(path_to_GENIE,self.center,"staging",fileName)
		narrowMafPath = os.path.join(path_to_GENIE,self.center,"staging","data_mutations_extended_%s_MAF_narrow.txt" % self.center)
		narrowMafColumns = [col['name'] for col in self.syn.getTableColumns(mafSynId) if col['name'] != 'inBED']

		#Need to get rid of this center specific code
		if self.center == "GRCC":
			temp = pd.read_csv(filePath, sep="\t",comment="#")
			temp.columns = [col.upper() for col in temp.columns]
			#if temp.get('CHROMOSOME') is not None:
			temp = temp[temp['CHROMOSOME'] != "WT"]
			temp['CHROMOSOME'] = [i.replace("chr","") for i in temp['CHROMOSOME']]
			temp['TUMOR_SEQ_ALLELE2'] = temp['TUMOR_SEQ_ALLELE1']
			del temp['TUMOR_SEQ_ALLELE1']
			filePath = os.path.join(path_to_GENIE,"GRCCmaf.txt")
			temp.to_csv(filePath,sep="\t", index=False)

		tempdir = os.path.join(path_to_GENIE, self.center)
		commandCall = ["perl",os.path.join(vcf2mafPath,"maf2maf.pl"),
					   "--input-maf",filePath,
					   "--output-maf",newMafPath,
					   "--vep-fork", '8',
					   "--tmp-dir",tempdir,
					   '--vep-path', veppath,
					   '--vep-data', vepdata,
					   #'--ref-fasta','/root/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa',
					   "--custom-enst", os.path.join(vcf2mafPath,"data/isoform_overrides_uniprot")]
		maf = subprocess.call(commandCall) 

		process_functions.rmFiles(tempdir, recursive=False)
		open(narrowMafPath,"w").close()
		if os.path.exists(newMafPath):
			#This needs to switch to streaming at some point
			mafDf = pd.read_csv(newMafPath,sep="\t",comment="#")
			mafDf = self.formatMAF(mafDf)
			self.createFinalMaf(mafDf, newMafPath, maf=True)
			narrowMafDf = mafDf[narrowMafColumns]
			self.createFinalMaf(narrowMafDf, narrowMafPath, maf=True)
			#These functions have to be next to each other, because no modifications can happen 
			#Store Narrow MAF into db
			self.storeProcessedMaf(narrowMafPath, mafSynId, centerMafSynId, isNarrow=True)
			#Store MAF flat file into synapse
			self.storeProcessedMaf(newMafPath, mafSynId, centerMafSynId)
		else:
			logger.error('ERROR PROCESSING %s' % filePath)
			filePath = "NOTPROCESSED"
		return(filePath)

	def process_steps(self, filePath, **kwargs): 
		processing = kwargs['processing']
		mutationFiles = []
		if processing == self._fileType:
			databaseToSynIdMappingDf = kwargs['databaseToSynIdMappingDf']
			vcf2mafPath = kwargs['vcf2mafPath']
			veppath = kwargs['veppath']
			vepdata = kwargs['vepdata']
			validMAFs = kwargs['validMAFs']
			path_to_GENIE = kwargs['path_to_GENIE']
			mafProcessing = "mafSP" if self._fileType == "mafSP" else 'vcf2maf'
			mafSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == mafProcessing][0]
			centerMafSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "centerMaf"][0]
			for filePath in validMAFs:
				mafFilePath = self.process_helper(filePath, path_to_GENIE, mafSynId, centerMafSynId,
											 vcf2mafPath, veppath, vepdata)
				mutationFiles.append(mafFilePath)
			logger.info("UPDATED DATABASE WITH: %s" % ", ".join(mutationFiles))
		else:
			logger.info("Please run with `--process %s` parameter if you want to reannotate the %s files" % (self._fileType, self._fileType))
		return(mutationFiles)


	def validate_helper(self, mutationDF, SP=False):
		"""
		This function validates the clinical file to make sure it adhere to the clinical SOP.
		
		:params filePath:     Path to mutation file
		:returns:             Text with all the errors in the clinical file
		"""

		first_header = ['CHROMOSOME','HUGO_SYMBOL','TUMOR_SAMPLE_BARCODE']
		if SP:
			correct_column_headers = ['CHROMOSOME','START_POSITION','REFERENCE_ALLELE','TUMOR_SAMPLE_BARCODE'] #T_REF_COUNT + T_ALT_COUNT = T_DEPTH
		else:
			correct_column_headers = ['CHROMOSOME','START_POSITION','REFERENCE_ALLELE','TUMOR_SAMPLE_BARCODE','T_ALT_COUNT'] #T_REF_COUNT + T_ALT_COUNT = T_DEPTH
		optional_headers = ['T_REF_COUNT','N_DEPTH','N_REF_COUNT','N_ALT_COUNT']
		tumors = ['TUMOR_SEQ_ALLELE2','TUMOR_SEQ_ALLELE1']
		

		mutationDF.columns = [col.upper() for col in mutationDF.columns]

		total_error = ""
		warning = ""
		#CHECK: First column must be in the first_header list
		if mutationDF.columns[0] not in first_header:
			total_error += "First column header must be one of these: %s.\n" % ", ".join(first_header)
		
		if not process_functions.checkColExist(mutationDF, "T_DEPTH") and not SP:
			if not process_functions.checkColExist(mutationDF, "T_REF_COUNT"):
				total_error += "If you are missing T_DEPTH, you must have T_REF_COUNT!\n"

		#CHECK: Everything in correct_column_headers must be in mutation file
		if not all([process_functions.checkColExist(mutationDF, i) for i in correct_column_headers]):
			total_error += "Your mutation file must at least have these headers: %s.\n" % ",".join([i for i in correct_column_headers if i not in mutationDF.columns.values])
		
		#CHECK: Must have either Tumor_Seq_Allele1 or 2
		tumor_cols = [i for i in tumors if i in mutationDF.columns.values]
		if len(tumor_cols) == 0:
			total_error += "Your mutation file must also have at least one of these headers: %s.\n" % " or ".join(tumors)
		nullTumorSeqCount = 0
		for i in tumor_cols:
			#CHECK: The value "NA" can't be used as a placeholder
			if sum(mutationDF[i].fillna('') == "NA") > 0:
				warning += "Your %s column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n" % i
			#CHECK: There can't be any null values
			if sum(mutationDF[i].isnull()) > 0:
				nullTumorSeqCount += 1

		if nullTumorSeqCount == len(tumor_cols) and nullTumorSeqCount > 0:
			total_error += "Your mutation file must have at least one of these: %s and at least one of those columns cannot have any empty values.\n" % ", ".join(tumors)

		#CHECK: Mutation file would benefit from columns in optional_headers
		if not all([process_functions.checkColExist(mutationDF, i) for i in optional_headers]) and not SP:
			warning += "Your mutation file does not have the column headers that can give extra information to the processed mutation file: %s.\n" % ", ".join([i for i in optional_headers if i not in mutationDF.columns.values ])      

		if process_functions.checkColExist(mutationDF, "REFERENCE_ALLELE"):
			if sum(mutationDF['REFERENCE_ALLELE'] == "NA") > 0:
				warning += "Your REFERENCE_ALLELE column contains NA values, which cannot be placeholders for blank values.  Please put in empty strings for blank values.\n"
			#CHECK: mutation file must not have empty reference or variant alleles
			if sum(mutationDF['REFERENCE_ALLELE'].isnull()) > 0:
				total_error += "Your mutation file cannot have any empty REFERENCE_ALLELE values.\n"

		return(total_error, warning)

	def validate_steps(self, filePathList, **kwargs):
		logger.info("VALIDATING %s" % os.path.basename(filePathList[0]))
		mutationDF = pd.read_csv(filePathList[0],sep="\t",comment="#",na_values = ['-1.#IND', '1.#QNAN', '1.#IND', 
								 '-1.#QNAN', '#N/A N/A', '#N/A', 'N/A', '#NA', 'NULL', 'NaN', 
								 '-NaN', 'nan','-nan',''],keep_default_na=False)
		total_error, warning = self.validate_helper(mutationDF)
		return(total_error, warning)

