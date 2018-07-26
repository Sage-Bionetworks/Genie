import subprocess
import os
import logging
import pandas as pd
import process_functions
import maf

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def contains_whitespace(x):
	"""
	Helper function for validateVCF.  No whitespace is allowed in VCF files

	:returns:     Sum of the the amount of whitespace in a string
	"""
	return(sum([" " in i for i in x if isinstance(i, str)]))


class vcf(maf.maf):

	_fileType = "vcf"

	_process_kwargs = ["validVCF", "processing", "path_to_GENIE", "databaseToSynIdMappingDf", 
					   "vcf2mafPath","veppath","vepdata"]

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]).startswith("GENIE-%s-" % self.center) and os.path.basename(filePath[0]).endswith(".vcf")
	
	def process_helper(self, vcffiles, path_to_GENIE, mafSynId, centerMafSynId,
					   vcf2mafPath,veppath, vepdata, 
					   reference="/home/tyu/reference/hg19/hg_19_all_chrs.fasta"):
		logger.info('VCF2MAF %s' % self.center)
		centerInputFolder = os.path.join(path_to_GENIE, self.center,"input")
		centerStagingFolder = os.path.join(path_to_GENIE,self.center,"staging")
		mafFiles = []
		for path in vcffiles:
			vcfName = os.path.basename(path)
			logger.info(vcfName)
			newVCFPath = os.path.join(centerInputFolder, vcfName)
			#remove chr from each row
			command = ["sed", "'s/^chr//'", path, ">", newVCFPath]
			subprocess.call(" ".join(command), shell=True)
			#Empty spaces must be replaced with a period
			command = ["sed", '-i', "'s/\t\t/\t.\t/g'", newVCFPath]
			subprocess.call(" ".join(command), shell=True)
			#All INFO/HGVS values have a whitespace, which is not allowed in VCF specs. Replace that with a comma
			command = ['sed', '-i', "'s/ p\./,p./'", newVCFPath]
			subprocess.call(" ".join(command), shell=True)
			#Strips out windows indentations \r
			command = ['dos2unix',newVCFPath]
			subprocess.call(command)
			vcfCols = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
			with open(newVCFPath,"r") as f:
				for line in f:
					if line.startswith("#CHROM"):
						cols = line

			cols = cols.replace("\n","")
			cols = cols.replace("\r","")
			cols = cols.split("\t")

			samples = [i for i in cols if i not in vcfCols]

			tumorName = vcfName.replace(".vcf","")

			if len(samples) == 1:
				tumor = samples[0]
				normal = "NORMAL"
				### If the tumor name isn't TUMOR, set the sample id to be the tumor name
				if tumor != "TUMOR":
					tumorName = tumor
			elif len(samples) == 2:
				#Tumor is always first, normal is second
				tumor = samples[0]
				normal = samples[1]
				tumorName = tumor   
			else:
				tumor = "TUMOR"
				normal = "NORMAL"
			newMAFPath = newVCFPath + ".maf"
			if os.path.isfile(newMAFPath):
				mafFiles.append(newMAFPath)
			else:
				command = ['perl', os.path.join(vcf2mafPath,'vcf2maf.pl'),
						   '--input-vcf', newVCFPath, 
						   '--output-maf', newMAFPath,
						   '--vep-path', veppath,
						   '--vep-data', vepdata,
						   '--vep-forks', '8',
						   '--tumor-id',tumorName,
						   '--normal-id',normal,
						   '--vcf-tumor-id',tumor,
						   #'--ref-fasta','/root/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa',
						   '--custom-enst', os.path.join(vcf2mafPath, 'data/isoform_overrides_uniprot')]
				subprocess.call(command)
				if (os.path.isfile(newMAFPath)):
					mafFiles.append(newMAFPath)
		
		logger.info("MERGING MAFS")
		#maf = pd.DataFrame()
		newMafPath = os.path.join(centerStagingFolder,"data_mutations_extended_%s_VCF.txt" % self.center)
		narrowMafPath = os.path.join(centerStagingFolder,"data_mutations_extended_%s_VCF_narrow.txt" % self.center)
		open(newMafPath,"w").close()
		open(narrowMafPath,"w").close()
		narrowMafColumns = [col['name'] for col in self.syn.getTableColumns(mafSynId) if col['name'] != 'inBED']

		for mafFile in mafFiles:
			mafDf = pd.read_csv(mafFile,sep="\t",comment="#")
			mafDf = self.formatMAF(mafDf)
			self.createFinalMaf(mafDf, newMafPath)
			narrowMafDf = mafDf[narrowMafColumns]
			self.createFinalMaf(narrowMafDf, narrowMafPath)

		if len(mafFiles) > 0:
			#These functions have to be next to each other, because no modifications can happen 
			#Store Narrow MAF into db
			self.storeProcessedMaf(narrowMafPath, mafSynId, centerMafSynId, isNarrow=True)
			#Store MAF flat file into synapse
			self.storeProcessedMaf(newMafPath, mafSynId, centerMafSynId)
		return(newMafPath)

	def process_steps(self, filePath, **kwargs): 
		processing = kwargs['processing']
		mutationFiles = []
		if processing == self._fileType:
			databaseToSynIdMappingDf = kwargs['databaseToSynIdMappingDf']
			vcf2mafPath = kwargs['vcf2mafPath']
			veppath = kwargs['veppath']
			vepdata = kwargs['vepdata']
			validVCF = kwargs['validVCF']
			path_to_GENIE = kwargs['path_to_GENIE']
			mafProcessing = "mafSP" if self._fileType == "mafSP" else 'vcf2maf'
			mafSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == mafProcessing][0]
			centerMafSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "centerMaf"][0]
			logger.info(validVCF)
			vcfFilePath = self.process_helper(validVCF, path_to_GENIE, mafSynId,centerMafSynId,
										 vcf2mafPath, veppath, vepdata)
			mutationFiles = [vcfFilePath]
			logger.info("UPDATED DATABASE WITH: %s" % ", ".join(mutationFiles))
		else:
			logger.info("Please run with `--process %s` parameter if you want to reannotate the %s files" % (self._fileType, self._fileType))
		return(mutationFiles)

	# Resolve missing read counts
	def _validate(self, vcf):
		REQUIRED_HEADERS = pd.Series(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])

		total_error = ""
		warning = ""
		if not all(REQUIRED_HEADERS.isin(vcf.columns)):
			total_error += "Your vcf file must have these headers: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.\n"

		if len(vcf.columns) > 8:
			if "FORMAT" not in vcf.columns:
				total_error += "Your vcf file must have FORMAT header if genotype columns exist.\n"
	   
		#Require that they report variants mapped to either GRCh37 or hg19 without 
		#the chr-prefix. variants on chrM are not supported
		haveColumn = process_functions.checkColExist(vcf, "#CHROM")
		if haveColumn:
			nochr = ["chr" in i for i in vcf['#CHROM'] if isinstance(i, str)]
			if sum(nochr) > 0:
				warning += "Your vcf file should not have the chr prefix in front of chromosomes.\n"
			if sum(vcf['#CHROM'].isin(["chrM"])) > 0:
				total_error += "Your vcf file must not have variants on chrM.\n"

		#No white spaces
		temp = vcf.apply(lambda x: contains_whitespace(x), axis=1)
		if sum(temp) >0:
			warning += "Your vcf file should not have any white spaces in any of the columns.\n"
		#I can also recommend a `bcftools query` command that will parse a VCF in a detailed way, 
		#and output with warnings or errors if the format is not adhered too
		return(total_error, warning)

	# Resolve missing read counts
	def validate_steps(self, filePathList, **kwargs):
		"""
		This function validates the VCF file to make sure it adhere to the genomic SOP.
		
		:params filePath:     Path to VCF file
		:returns:             Text with all the errors in the VCF file
		"""  
		filePath = filePathList[0]
		logger.info("VALIDATING %s" % os.path.basename(filePath))
		#FORMAT is optional
		headers = None
		with open(filePath,"r") as foo:
			for i in foo:
				if i.startswith("#CHROM"):
					headers = i.replace("\n","").replace("\r","").split("\t")
		if headers is not None:
			vcf = pd.read_csv(filePath, sep="\t",comment="#",header=None,names=headers)
		else:
			raise ValueError("Your vcf must start with the header #CHROM")

		total_error, warning = self._validate(vcf)
		return(total_error, warning)