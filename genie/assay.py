from __future__ import absolute_import
from genie import example_filetype_format
import os
import logging
import synapseclient
import yaml
logger = logging.getLogger(__name__)


class Assayinfo(example_filetype_format.FileTypeFormat):
	'''
	Assay information file type
	'''
	_fileType = "assayinfo"

	_process_kwargs = ["newPath", "databaseSynId"]

	def _validateFilename(self, filepath_list):
		assert os.path.basename(filepath_list[0]) == "assay_information.txt"

	def process_steps(self, filePath, *args, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		databaseSynId = kwargs['databaseSynId']
		assay_info_df = self._get_dataframe(filePath)
		process_assay_info_df = self._process()

		process_functions.updateData(self.syn, sampleSynId, sampleClinical, self.center, col=sampleCols, toDelete=True)

		return(filePath)

	def _process(self, seg):
		seg.columns = [col.upper() for col in seg.columns]
		newsamples = [process_functions.checkGenieId(i, self.center) for i in seg['ID']]
		seg['ID'] = newsamples
		seg = seg.drop_duplicates()
		seg = seg.rename(columns= {'LOC.START':'LOCSTART','LOC.END':'LOCEND','SEG.MEAN':'SEGMEAN','NUM.MARK':'NUMMARK'})
		seg['CHROM'] = [str(chrom).replace("chr","") for chrom in seg['CHROM']]
		seg['CENTER'] = self.center
		seg['LOCSTART'] = seg['LOCSTART'].astype(int)
		seg['LOCEND'] = seg['LOCEND'].astype(int)
		seg['NUMMARK'] = seg['NUMMARK'].astype(int)
		return(seg)


	def _get_dataframe(self, filePathList):
		'''
		Return a dict
		'''
		filePath = filePathList[0]
		with open(filePath, 'r') as yamlfile:
			panel_info_dict = yaml.load(yamlfile)
		assay_info_df = pd.DataFrame(panel_info_dict)
		assay_info_df = assay_info_df.transpose()
		assay_info_df['SEQ_ASSAY_ID'] = assay_info_df.index
		assay_info_df.reset_index(drop=True, inplace=True)
		return(assay_info_df)

	# def process_steps(self, filePath, **kwargs):
	# 	logger.info('PROCESSING %s' % filePath)
	# 	newPath = kwargs['newPath']
	# 	databaseSynId = kwargs['databaseSynId']


	# 	seg = pd.read_csv(filePath, sep="\t")
	# 	seg = self._process(seg)
	# 	process_functions.updateData(self.syn, databaseSynId, seg, self.center, toDelete=True)
	# 	seg.to_csv(newPath,sep="\t",index=False)
	# 	return(newPath)

	def _validate(self, assay_dict):
		total_error = ""
		warning = ""
		segDF.columns = [col.upper() for col in segDF.columns]

		REQUIRED_HEADERS = pd.Series(['ID','CHROM','LOC.START','LOC.END','NUM.MARK','SEG.MEAN'])
		
		if not all(REQUIRED_HEADERS.isin(segDF.columns)):
			total_error += "Your seg file is missing these headers: %s.\n" % ", ".join(REQUIRED_HEADERS[~REQUIRED_HEADERS.isin(segDF.columns)])
		else:
			intCols = ['LOC.START','LOC.END','NUM.MARK']
			nonInts = [col for col in intCols if segDF[col].dtype != int]
			if len(nonInts) > 0:
				total_error += "Seg: Only integars allowed in these column(s): %s.\n" % ", ".join(sorted(nonInts))
			if not segDF['SEG.MEAN'].dtype in [float, int]:
				total_error += "Seg: Only numerical values allowed in SEG.MEAN.\n"

		checkNA = segDF.isna().apply(sum)
		nullCols = [ind for ind in checkNA.index if checkNA[ind] > 0]
		if len(nullCols) > 0:
			total_error += "Seg: No null or empty values allowed in column(s): %s.\n" % ", ".join(sorted(nullCols))

		return(total_error, warning)