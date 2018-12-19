import pandas as pd
import logging
import os
logger = logging.getLogger(__name__)

class FileTypeFormat(object):

	_process_kwargs = ["newPath", "databaseSynId"]

	_fileType = "fileType"

	_validation_kwargs = []

	def __init__(self, syn, center, poolSize=1):
		self.syn = syn
		self.center = center
		#self.pool = multiprocessing.Pool(poolSize)
	
	def _get_dataframe(self, filePathList):
		'''
		This function by defaults assumes the filePathList is length of 1 
		and is a tsv file.  Could change depending on file type.
		'''
		filePath = filePathList[0]
		df = pd.read_csv(filePath,sep="\t",comment="#")
		return(df)

	def read_file(self, filePathList):
		'''
		Each file is to be read in for validation and processing.
		This is not to be changed in any functions.
		'''
		df = self._get_dataframe(filePathList)
		return(df)

	def _validateFilename(self, filePath):
		pass


	def validateFilename(self, filePath):
		self._validateFilename(filePath)
		return(self._fileType)



	def process_steps(self, df, **kwargs):
		pass

	def preprocess(self, filePath, **kwargs):
# - clinical
# - maf
# - vcf
		return(dict())


	def process(self, filePath, **kwargs):
		preprocess_args = self.preprocess(filePath, **kwargs)
		kwargs.update(preprocess_args)
		mykwargs = {}
		for required_parameter in self._process_kwargs:
			assert required_parameter in kwargs.keys(), "%s not in parameter list" % required_parameter
			mykwargs[required_parameter] = kwargs[required_parameter]
		logger.info('PROCESSING %s' % filePath)
		#If file type is vcf or maf file, processing requires a filepath
		if self._fileType not in ['vcf','maf','mafSP','md','clinical']:
			path_or_df = self.read_file([filePath])
		else:
			path_or_df = filePath
		path = self.process_steps(path_or_df, **mykwargs)
		return(path)

	def _validate(self, df, **kwargs):
		total_error =""
		warning = ""
		logger.info("NO VALIDATION for %s files" % self._fileType)
		return(total_error, warning)

	def validate(self, filePathList, **kwargs):
		mykwargs = {}
		for required_parameter in self._validation_kwargs:
			assert required_parameter in kwargs.keys(), "%s not in parameter list" % required_parameter
			mykwargs[required_parameter] = kwargs[required_parameter]
		logger.info("VALIDATING %s" % os.path.basename(",".join(filePathList)))
		df = self.read_file(filePathList)
		total_error, warning = self._validate(df, **mykwargs)
		return(total_error, warning)