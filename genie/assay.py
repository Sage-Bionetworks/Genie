from __future__ import absolute_import
from genie import example_filetype_format, process_functions
import os
import logging
import synapseclient
import yaml
import pandas as pd
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

	def _process(self, df):
		()
		#sequencing_center


	def _get_dataframe(self, filepath_list):
		'''
		Takes in yaml file, returns dataframe
		'''
		filepath = filepath_list[0]
		with open(filepath, 'r') as yamlfile:
			panel_info_dict = yaml.load(yamlfile)
		assay_info_df = pd.DataFrame(panel_info_dict)
		assay_info_df = assay_info_df.transpose()
		assay_info_df['SEQ_ASSAY_ID'] = assay_info_df.index
		assay_info_df.reset_index(drop=True, inplace=True)
		return(assay_info_df)

	def _validate(self, assay_info_df):
		total_error = ""
		warning = ""

		if process_functions.checkColExist(assay_info_df, "SEQ_ASSAY_ID"):
			all_seq_assays = assay_info_df.SEQ_ASSAY_ID.unique()
			if not all([assay.startswith(self.center) for assay in all_seq_assays]):
				total_error += "Assay_information.yaml: Please make sure your all your SEQ_ASSAY_IDs start with your center abbreviation.\n"
		else:
			total_error += "Assay_information.yaml: Must have SEQ_ASSAY_ID column.\n"

		warn, error = process_functions.check_col_and_values(assay_info_df, 'is_paired_end', [True, False], filename="Assay_information.yaml", required=True)
		warning += warn
		total_error  += error 
		warn, error = process_functions.check_col_and_values(assay_info_df, 'library_selection', ['Hybrid Selection','PCR','Affinity Enrichment','Poly-T Enrichment','Random','rRNA Depletion','miRNA Size Fractionation','Other'], filename="Assay_information.yaml", required=True)
		warning += warn
		total_error  += error 
		warn, error = process_functions.check_col_and_values(assay_info_df, 'library_strategy', ['ATAC-Seq','Bisulfite-Seq','ChIP-Seq','miRNA-Seq','RNA-Seq','Targeted Sequencing','WGS','WXS'], filename="Assay_information.yaml", required=True)
		warning += warn
		total_error  += error 
		warn, error = process_functions.check_col_and_values(assay_info_df, 'platform', ['Illumina','SOLiD','LS454','Ion Torrent','Complete Genomics','PacBio','Other'], filename="Assay_information.yaml", required=True)
		warning += warn
		total_error  += error 

		instrument_model = ['454 GS FLX Titanium','AB SOLiD 4','AB SOLiD 2','AB SOLiD 3','Complete Genomics','Illumina HiSeq X Ten',
				'Illumina HiSeq X Five','Illumina Genome Analyzer II','Illumina Genome Analyzer IIx','Illumina HiSeq 2000',
				'Illumina HiSeq 2500','Illumina HiSeq 4000','Illumina MiSeq','Illumina NextSeq','Ion Torrent PGM',
				'Ion Torrent Proton','PacBio RS','Other',None]
		warn, error = process_functions.check_col_and_values(assay_info_df, 'instrument_model', instrument_model, filename="Assay_information.yaml", required=True)
		warning += warn
		total_error  += error 

		variant_classes = ['Splice_Site','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins',
		'Nonstop_Mutation','Translation_Start_Site','In_Frame_Ins','In_Frame_Del','Missense_Mutation',
		'Intron','Splice_Region','Silent','RNA',"5'UTR","3'UTR",'IGR',"5'Flank","3'Flank"]

		warn, error = process_functions.check_col_and_values(assay_info_df, 'variant_consequences', variant_classes, filename="Assay_information.yaml")
		warning += warn
		total_error  += error 

		if not process_functions.checkColExist(assay_info_df, "target_capture_kit"):
			total_error += "Assay_information.yaml: Must have target_capture_kit column.\n"

		if process_functions.checkColExist(assay_info_df, "read_length"):
			if not all([process_functions.checkInt(i) for i in assay_info_df["read_length"]]):
				total_error += "Assay_information.yaml: Please double check your read_length.  It must be an integer or 'Unknown'.\n"
		else:
			total_error += "Assay_information.yaml: Must have read_length column.\n"

		if process_functions.checkColExist(assay_info_df, "number_of_genes"):
			if not all([process_functions.checkInt(i) for i in assay_info_df["number_of_genes"]]):
				total_error += "Assay_information.yaml: Please double check your number_of_genes.  It must be an integer or 'Unknown'.\n"
		else:
			total_error += "Assay_information.yaml: Must have number_of_genes column.\n"
		
		if process_functions.checkColExist(assay_info_df, "gene_padding"):
			if not all([process_functions.checkInt(i) for i in assay_info_df["gene_padding"] if i is not None and not pd.isnull(i)]):
				total_error += "Assay_information.yaml: Please double check your gene_padding.  It must be an integer or blank.\n"
		else:
			warning += "Assay_information.yaml: gene_padding is by default 10 if not specified.\n"

		return(total_error, warning)