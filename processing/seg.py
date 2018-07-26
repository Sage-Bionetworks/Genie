import logging
import os
import pandas as pd
import process_functions
import example_filetype_format

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class seg(example_filetype_format.FileTypeFormat):

	_fileType = "seg"

	_process_kwargs = ["newPath", "databaseSynId"]

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "genie_data_cna_hg19_%s.seg" % self.center

	def process_steps(self, filePath, **kwargs):
		#For CBS files
		if kwargs.get("path") is not None:
			filePath = kwargs['path']
			newPath = filePath
		else:
			newPath = kwargs['newPath']
		logger.info('PROCESSING %s' % filePath)
		databaseSynId = kwargs['databaseSynId']
		seg = pd.read_csv(filePath, sep="\t")
		seg.columns = [col.upper() for col in seg.columns]
		newsamples = [process_functions.checkGenieId(i, self.center) for i in seg['ID']]
		seg['ID'] = newsamples
		seg = seg.drop_duplicates()
		seg = seg.rename(columns= {'LOC.START':'LOCSTART','LOC.END':'LOCEND','SEG.MEAN':'SEGMEAN','NUM.MARK':'NUMMARK'})
		seg['CENTER'] = self.center
		seg['LOCSTART'] = seg['LOCSTART'].astype(int)
		seg['LOCEND'] = seg['LOCEND'].astype(int)
		seg['NUMMARK'] = seg['NUMMARK'].astype(int)
		process_functions.updateData(self.syn, databaseSynId, seg, self.center, seg.columns)
		seg.to_csv(newPath,sep="\t",index=False)
		return(newPath)

	def validate_steps(self, filePathList, **kwargs):
		"""
		This function validates the SEG file to make sure it adhere to the genomic SOP.
		
		:params filePath:     Path to SEG/CBS files

		:returns:             Text with all the errors in the SEG file
		"""
		filePath = filePathList[0]
		logger.info("VALIDATING %s" % os.path.basename(filePath))
		total_error = ""
		warning = ""

		segDF = pd.read_csv(filePath,sep="\t",comment="#")
		segDF.columns = [col.upper() for col in segDF.columns]

		REQUIRED_HEADERS = ['ID','CHROM','LOC.START','LOC.END','NUM.MARK','SEG.MEAN']
		
		if not all(segDF.columns.isin(REQUIRED_HEADERS)):
			total_error += "Your seg file must at least have these headers: %s.\n" % ",".join([i for i in REQUIRED_HEADERS if i not in segDF.columns.values])

		return(total_error, warning)