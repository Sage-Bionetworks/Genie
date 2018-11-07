from __future__ import absolute_import
from genie import seg
from genie import example_filetype_format
from genie import process_functions

import os
import logging
import pandas as pd
#logging.basicConfig(level=logging.INFO)
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

	def _validate(self, cbsDf):
		total_error = ""
		warning = ""
		REQUIRED_HEADERS = pd.Series([cbsDf.columns[0], 'Chr','Start','End','Probes','Log2Ratio'])
		if not all(REQUIRED_HEADERS.isin(cbsDf.columns)):
			total_error += "Your cbs file must at least have these headers: %s.\n" % ",".join(REQUIRED_HEADERS[~REQUIRED_HEADERS.isin(cbsDf.columns)])
		else:
			intCols = ['Start','End','Probes']
			nonInts = [col for col in intCols if cbsDf[col].dtype != int]
			if len(nonInts) > 0:
				total_error += "cbs: Only integars allowed in these column(s): %s.\n" % ", ".join(nonInts)
			if not cbsDf['Log2Ratio'].dtype in [float, int]:
				total_error += "cbs: Only numerical values allowed in Log2Ratio.\n"

		checkNA = cbsDf.isna().apply(sum)
		nullCols = [ind for ind in checkNA.index if checkNA[ind] > 0]
		if len(nullCols) > 0:
			total_error += "cbs: No null or empty values allowed in column(s): %s.\n" % ", ".join(nullCols)
		return(total_error, warning)
		
	def validate_steps(self, filePathList, **kwargs):

		filePath = filePathList[0]
		logger.info("VALIDATING %s" % os.path.basename(filePath))


		cbsDf = pd.read_csv(filePath,sep="\t",comment="#")
		return(self._validate(cbsDf))