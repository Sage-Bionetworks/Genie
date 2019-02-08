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
		assert os.path.basename(filepath_list[0]) == "assay_information.yaml"

	def process_steps(self, filePath, *args, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		databaseSynId = kwargs['databaseSynId']
		assay_info_df = self._get_dataframe(filePath)
		process_assay_info_df = self._process()
		process_functions.updateData(self.syn, databaseSynId, process_assay_info_df, self.center, filterByColumn="CENTER", toDelete=True)

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

	def _validate(self, assay_info_df):
		total_error = ""
		warning = ""
		#assay_info_df.columns = [col.upper() for col in assay_info_df.columns]

		required_headers = pd.Series(['SEQ_ASSAY_ID','is_paired_end','library_selection','library_strategy','platform','read_length','target_capture_kit','instrument_model','number_of_genes'])
		
		if not all(required_headers.isin(assay_info_df.columns)):
			total_error += "assay_information.yaml: Missing headers: %s.\n" % ", ".join(REQUIRED_HEADERS[~REQUIRED_HEADERS.isin(assay_info_df.columns)])
		else:
			all_seq_assays = assay_info_df.SEQ_ASSAY_ID.unique()
			all_seq_assays.apply(lambda assay: assay.startswith(self.center))


			assay_info_df.is_paired_end.isin([True, False]).all()
			assay_info_df.library_selection.isin(['Hybrid Selection','PCR','Affinity Enrichment','Poly-T Enrichment','Random','rRNA Depletion','miRNA Size Fractionation','Other']).all()
			assay_info_df.library_strategy.isin(['ATAC-Seq','Bisulfite-Seq','ChIP-Seq','miRNA-Seq','RNA-Seq','Targeted Sequencing','WGS','WXS']).all()
			assay_info_df.platform.isin(['Illumina','SOLiD','LS454','Ion Torrent','Complete Genomics','PacBio','Other']).all()
			assay_info_df.read_length.isin([True, False]).all()
			#assay_info_df.target_capture_kit.isin([True, False]).all()
			assay_info_df.instrument_model.isin([True, False]).all()
			assay_info_df.number_of_genes.isin([True, False]).all()




		return(total_error, warning)