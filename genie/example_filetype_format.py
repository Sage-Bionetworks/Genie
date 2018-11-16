
class FileTypeFormat(object):

	_process_kwargs = ["newPath", "databaseSynId"]

	_fileType = "fileType"

	_validation_kwargs = []

	def __init__(self, syn, center, poolSize=1):
		self.syn = syn
		self.center = center
		#self.pool = multiprocessing.Pool(poolSize)


	def _validateFilename(self, filePath):
		pass


	def validateFilename(self, filePath):
		self._validateFilename(filePath)
		return(self._fileType)



	def process_steps(self, filePath, *args, **kwargs):
		pass

	def preprocess(self, filePath, *args, **kwargs):
# - clinical
# - maf
# - vcf
		return(dict())


	def process(self, filePath, *args, **kwargs):
		
		preprocess_args = self.preprocess(filePath, **kwargs)
		kwargs.update(preprocess_args)
		mykwargs = {}
		for required_parameter in self._process_kwargs:
			assert required_parameter in kwargs.keys(), "%s not in parameter list" % required_parameter
			mykwargs[required_parameter] = kwargs[required_parameter]

		path = self.process_steps(filePath, **mykwargs)
		return(path)

	def validate_steps(self, filePathList, **kwargs):
		total_error = ""
		warning = ""
		return(total_error, warning)

	def validate(self, filePathList, **kwargs):
		mykwargs = {}
		for required_parameter in self._validation_kwargs:
			assert required_parameter in kwargs.keys(), "%s not in parameter list" % required_parameter
			mykwargs[required_parameter] = kwargs[required_parameter]

		total_error, warning = self.validate_steps(filePathList, **mykwargs)
		return(total_error, warning)