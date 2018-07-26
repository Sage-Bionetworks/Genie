import os
import seg
import logging
import pandas as pd
import example_filetype_format
import process_functions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class cbs(seg.seg):

	_fileType = "cbs"

	_process_kwargs = ["validCBS", "path_to_GENIE", "databaseSynId", "path"]

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]).endswith("cbs") and os.path.basename(filePath[0]).startswith("GENIE-%s" % self.center)

	def preprocess(self, filePath, **kwargs):
		validCBS = kwargs['validCBS']
		path_to_GENIE = kwargs['path_to_GENIE']
		segDF = pd.DataFrame(columns = ['ID','CHROM','LOC.START','LOC.END','NUM.MARK','SEG.MEAN'])
		path = os.path.join(path_to_GENIE,"%s.seg" % self.center)
		for cbsFile in validCBS:
			newseg = pd.read_csv(cbsFile,sep="\t")
			newseg = newseg.rename(columns={newseg.columns[0]:'ID','Chr':'CHROM','Start':'LOC.START','End':'LOC.END','Probes':'NUM.MARK','Log2Ratio':'SEG.MEAN'})
			newseg['ID'] = os.path.basename(cbsFile).replace(".cbs","")
			segDF = segDF.append(newseg)
		noFloatsDf = process_functions.removeFloat(segDF)
		with open(path, "w") as segFile:
			segFile.write(noFloatsDf)
		return({"path":path})

	def validate_steps(self, filePathList, **kwargs):

		filePath = filePathList[0]
		logger.info("VALIDATING %s" % os.path.basename(filePath))
		total_error = ""
		warning = ""

		cbsDf = pd.read_csv(filePath,sep="\t",comment="#")
		REQUIRED_HEADERS = [cbsDf.columns[0], 'Chr','Start','End','Probes','Log2Ratio']
		if not all(cbsDf.columns.isin(REQUIRED_HEADERS)):
			total_error += "Your cbs file must at least have these headers: %s.\n" % ",".join([i for i in REQUIRED_HEADERS if i not in cbsDf.columns.values])

		return(total_error, warning)