#!/usr/local/bin/python
import pandas as pd 
import argparse
import synapseclient
from synapseclient import File, Link
import os
import process_functions as process
import logging
import math
import datetime
import time
import re
import subprocess
import synapseutils as synu
import create_case_lists
import dashboard_table_updater

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GENIE_RELEASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),"GENIE_release")
CASE_LIST_PATH = os.path.join(GENIE_RELEASE_DIR,'case_lists')
CNA_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR,"data_CNA_%s.txt")
SAMPLE_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR,'data_clinical_supp_sample_%s.txt')
PATIENT_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR,'data_clinical_supp_patient_%s.txt')
MUTATIONS_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR,'data_mutations_extended_%s.txt')
FUSIONS_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR,'data_fusions_%s.txt')
SEG_CENTER_PATH = os.path.join(GENIE_RELEASE_DIR,'genie_data_cna_hg19_%s.seg')
BED_DIFFS_SEQASSAY_PATH = os.path.join(GENIE_RELEASE_DIR,'diff_%s.csv')

def findCaseListId(syn, parentId):
	releaseEnts = synu.walk(syn, parentId)
	releaseFolders = next(releaseEnts)
	if len(releaseFolders[1]) == 0:
		caselistId = syn.store(synapseclient.Folder(name="case_lists", parent=parentId)).id
	else:
		caselistId = releaseFolders[1][0][1]
	return(caselistId)

def storeFile(syn, filePath, genieVersion="database", name=None, parent=None, fileFormat = None, cBioFileFormat = None, staging=False, caseLists=False, centerStaging=False):
	#Take care of clinical case (Clinical files go elsewhere)
	# ANONYMIZE_CENTER = syn.tableQuery('SELECT * FROM syn10170510')
	# ANONYMIZE_CENTER_DF = ANONYMIZE_CENTER.asDataFrame()

	if staging:
		parent = "syn9689654" if caseLists else "syn7551278"
	logger.info("STORING FILE: %s" % os.path.basename(filePath))
	# if not centerStaging:
	# 	process.center_anon(filePath, ANONYMIZE_CENTER_DF)
	if name is None:
		name = os.path.basename(filePath)
	temp = File(filePath, name=name, parent = parent, versionComment=genieVersion)
	if fileFormat is not None:
		temp.fileFormat = fileFormat
	if cBioFileFormat is not None:
		temp.cBioFileFormat = cBioFileFormat
	temp = syn.store(temp)
	# if not centerStaging:
	# 	process.center_convert_back(filePath, ANONYMIZE_CENTER_DF)
	return(temp)

#Remove PHI data
def reAnnotatePHI(mergedClinical):
	PHIcutoff = 365*89
	pedsCutoff = 365*18
	mergedClinical['AGE_AT_SEQ_REPORT'] = mergedClinical['AGE_AT_SEQ_REPORT'].fillna('')
	mergedClinical['BIRTH_YEAR'] = mergedClinical['BIRTH_YEAR'].fillna('')
	AGE_AT_SEQ_REPORT =[pedsCutoff - 1 if "<" in str(age) else age for age in mergedClinical['AGE_AT_SEQ_REPORT']]
	AGE_AT_SEQ_REPORT =[PHIcutoff + 1 if ">" in str(age) else age for age in AGE_AT_SEQ_REPORT]
	toRedact  = [int(float(age)) > PHIcutoff if age != '' else False for age in AGE_AT_SEQ_REPORT]
	toRedactPeds  = [int(float(age)) < pedsCutoff if age != '' else False for age in AGE_AT_SEQ_REPORT]
	logger.info("Redacting >89")
	logger.info(mergedClinical[toRedact])
	logger.info("Redacting <18")
	logger.info(mergedClinical[toRedactPeds])
	#mergedClinical['BIRTH_YEAR'][toRedact] = "<1926"
	mergedClinical['BIRTH_YEAR'][toRedact] = "cannotReleaseHIPAA"
	mergedClinical['AGE_AT_SEQ_REPORT'][toRedact] = ">32485"
	#mergedClinical['BIRTH_YEAR'][toRedactPeds] = ">1998"
	mergedClinical['BIRTH_YEAR'][toRedactPeds] = "withheld"
	mergedClinical['AGE_AT_SEQ_REPORT'][toRedactPeds] = "<6570"
	#BIRTH_YEAR =[">89y" if str(year) in [">89","<1926"] else year for year in mergedClinical['BIRTH_YEAR']]
	#BIRTH_YEAR =["<18y" if str(year) in ["<18",">1998"] else year for year in BIRTH_YEAR]
	#mergedClinical['BIRTH_YEAR'] = BIRTH_YEAR
	mergedClinical['INT_CONTACT'] = mergedClinical['INT_CONTACT'].fillna('')
	mergedClinical['INT_DOD'] = mergedClinical['INT_DOD'].fillna('')
	INT_CONTACT =[pedsCutoff - 1 if "<" in str(age) else age for age in mergedClinical['INT_CONTACT']]
	INT_CONTACT =[PHIcutoff + 1 if ">" in str(age) else age for age in INT_CONTACT]
	toRedact  = [int(float(age)) > PHIcutoff if age != '' else False for age in INT_CONTACT]
	toRedactPeds  = [int(float(age)) < pedsCutoff if age != '' else False for age in INT_CONTACT]
	mergedClinical['INT_CONTACT'][toRedact] = ">32485"
	mergedClinical['INT_CONTACT'][toRedactPeds] = "<6570"
	mergedClinical['BIRTH_YEAR'][toRedact] = "cannotReleaseHIPAA"
	mergedClinical['BIRTH_YEAR'][toRedactPeds] = "withheld"

	INT_DOD =[pedsCutoff - 1 if "<" in str(age) else age for age in mergedClinical['INT_DOD']]
	INT_DOD =[PHIcutoff + 1 if ">" in str(age) else age for age in INT_DOD]
	toRedact  = [int(float(age)) > PHIcutoff if age != '' else False for age in INT_DOD]
	toRedactPeds  = [int(float(age)) < pedsCutoff if age != '' else False for age in INT_DOD]
	mergedClinical['INT_DOD'][toRedact] = ">32485"
	mergedClinical['INT_DOD'][toRedactPeds] = "<6570"
	mergedClinical['BIRTH_YEAR'][toRedact] = "cannotReleaseHIPAA"
	mergedClinical['BIRTH_YEAR'][toRedactPeds] = "withheld"

	#Redact DFCI inputed fields
	mergedClinical['BIRTH_YEAR'] = ["cannotReleaseHIPAA" if ">" in str(year) else year for year in mergedClinical['BIRTH_YEAR']]
	mergedClinical['BIRTH_YEAR'] = ["withheld" if "<" in str(year) else year for year in mergedClinical['BIRTH_YEAR']]

	return(mergedClinical)

#remove string float
def removeStringFloat(string):
	string = string.replace(".0\t","\t")
	string = string.replace(".0\n","\n")
	return(string)

#remove decimal for integers due to pandas
def removePandasDfFloat(df):
	text = df.to_csv(sep="\t",index=False)
	text = removeStringFloat(text)
	return(text)

#Configure each maf row
def configureMafRow(rowArray, headers, keepSamples, remove_variants):
	chrom = str(rowArray[headers.index('Chromosome')])
	start = str(rowArray[headers.index('Start_Position')])
	end = str(rowArray[headers.index('End_Position')])
	ref = str(rowArray[headers.index('Reference_Allele')])
	seq = str(rowArray[headers.index('Tumor_Seq_Allele2')])
	sampleId = str(rowArray[headers.index('Tumor_Sample_Barcode')])
	variant = chrom +' '+ start+ ' '+end +' '+ref + ' '+ seq+ ' ' + sampleId
	if pd.Series(sampleId).isin(keepSamples).any() and not pd.Series(variant).isin(remove_variants).any():
		fillnas = ['t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count']
		for i in fillnas:
			#mutationsDf[i] = mutationsDf[i].fillna("NA")
			#mutationsDf[i] = ["NA" if str(each) == "." else each for each in  mutationsDf[i]]
			value = rowArray[headers.index(i)]
			rowArray[headers.index(i)] = "" if str(value) == "." else value

		nDepth = rowArray[headers.index("n_depth")]
		rowArray[headers.index("Match_Norm_Seq_Allele2")] = '' if str(nDepth) in ["NA","0.0"] else nDepth
		rowArray[headers.index("Match_Norm_Seq_Allele1")] = '' if str(nDepth) in ["NA","0.0"] else nDepth
		#rowArray.pop(headers.index('inBED'))
		newRow = "\t".join(rowArray)
		newRow += "\n"
		newRow = removeStringFloat(newRow)
		return(newRow)
	else:
		return(None)

#Run MAF in BED script, filter data and update MAFinBED database
def runMAFinBED(syn, CENTER_MAPPING_DF, databaseSynIdMappingDf, test=False, genieVersion="test"):
	MAFinBED_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../analyses/genomicData/MAFinBED.R')
	command = ['Rscript', MAFinBED_script, str(test)]
	subprocess.check_call(command)

	mutationSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "vcf2maf"][0]
	removedVariants = syn.tableQuery("select Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Center from %s where inBED is False and Center in ('%s')" % (mutationSynId,"','".join(CENTER_MAPPING_DF.center)))
	removedVariantsDf = removedVariants.asDataFrame()
	removedVariantsDf['removeVariants'] = removedVariantsDf['Chromosome'].astype(str) + ' ' + removedVariantsDf['Start_Position'].astype(str) + ' ' + removedVariantsDf['End_Position'].astype(str) + ' ' + removedVariantsDf['Reference_Allele'].astype(str) + ' ' + removedVariantsDf['Tumor_Seq_Allele2'].astype(str) + ' ' + removedVariantsDf['Tumor_Sample_Barcode'].astype(str)
	#Store filtered vairants
	for center in removedVariantsDf['Center'].unique():
		center_mutation = removedVariantsDf[removedVariantsDf['Center'] == center]

		#mafText = removePandasDfFloat(center_mutation)
		center_mutation.to_csv("mafinbed_filtered_variants.csv",index=False)

		storeFile(syn, "mafinbed_filtered_variants.csv", parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True, genieVersion=genieVersion)
		os.unlink("mafinbed_filtered_variants.csv")
	return(removedVariantsDf['removeVariants'])

def seq_date_filter(clinicalDf, processingDate, consortiumReleaseCutOff):
	removeSeqDateSamples = process.seqDateFilter(clinicalDf, processingDate, consortiumReleaseCutOff)
	return(removeSeqDateSamples)

def mutation_in_cis_filter(syn, skipMutationsInCis, test, variant_filtering_synId, CENTER_MAPPING_DF, genieVersion):
	if not skipMutationsInCis:
		mergeCheck_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../analyses/mergeFlag/mergeCheck.R')
		command = ['Rscript', mergeCheck_script, str(test)]
		subprocess.check_call(command)
		# Store each centers mutations in cis to their staging folder
		mergeCheck = syn.tableQuery("select * from %s where Center in ('%s')" % (variant_filtering_synId,"','".join(CENTER_MAPPING_DF.center)))
		mergeCheckDf = mergeCheck.asDataFrame()
		for center in mergeCheckDf.Center.unique():
			if not pd.isnull(center):
				stagingSynId = CENTER_MAPPING_DF.stagingSynId[CENTER_MAPPING_DF['center'] == center]
				mergeCheckDf[mergeCheckDf['Center'] == center].to_csv("mutationsInCis_filtered_samples.csv",index=False)
				storeFile(syn, "mutationsInCis_filtered_samples.csv", parent = stagingSynId[0], centerStaging=True, genieVersion=genieVersion)
				os.unlink("mutationsInCis_filtered_samples.csv")
	variant_filtering = syn.tableQuery("SELECT Tumor_Sample_Barcode FROM %s where Flag = 'TOSS' and Tumor_Sample_Barcode is not null" % variant_filtering_synId)

	filtered_samples = variant_filtering.asDataFrame()
	# #Alex script #1 removed patients
	remove_samples = filtered_samples['Tumor_Sample_Barcode'].drop_duplicates()
	return(remove_samples)

def seq_assay_id_filter(clinicalDf):
	remove_seqassayId = clinicalDf['SEQ_ASSAY_ID'].value_counts()[clinicalDf['SEQ_ASSAY_ID'].value_counts() < 50]
	clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'].isin(remove_seqassayId.keys().values)]
	return(clinicalDf.SAMPLE_ID)


#def stagingToCbio(syn, releaseSynId, releaseDate, filtering=False, release=False, current_release_staging=False):	
def stagingToCbio(syn, processingDate, genieVersion, CENTER_MAPPING_DF, databaseSynIdMappingDf, oncotree_url=None, consortiumReleaseCutOff=183, filtering=False, current_release_staging=False, skipMutationsInCis=False, test=False):	
	CNA_PATH = os.path.join(GENIE_RELEASE_DIR,"data_CNA_%s.txt" % genieVersion)
	CLINCICAL_PATH = os.path.join(GENIE_RELEASE_DIR,'data_clinical_%s.txt' % genieVersion)
	CLINCICAL_SAMPLE_PATH = os.path.join(GENIE_RELEASE_DIR,'data_clinical_sample_%s.txt' % genieVersion)
	CLINCICAL_PATIENT_PATH = os.path.join(GENIE_RELEASE_DIR,'data_clinical_patient_%s.txt' % genieVersion)
	DATA_GENE_PANEL_PATH = os.path.join(GENIE_RELEASE_DIR,'data_gene_matrix_%s.txt' % genieVersion)
	MUTATIONS_PATH = os.path.join(GENIE_RELEASE_DIR,'data_mutations_extended_%s.txt' % genieVersion)
	FUSIONS_PATH = os.path.join(GENIE_RELEASE_DIR,'data_fusions_%s.txt' % genieVersion)
	SEG_PATH = os.path.join(GENIE_RELEASE_DIR,'genie_private_data_cna_hg19_%s.seg' % genieVersion)
	COMBINED_BED_PATH = os.path.join(GENIE_RELEASE_DIR,'genie_combined_%s.bed' % genieVersion)
	consortiumReleaseSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "consortium"][0]
	variant_filtering_synId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "mutationsInCis"][0]

	if not os.path.exists(GENIE_RELEASE_DIR):
		os.mkdir(GENIE_RELEASE_DIR)

	### ADD CHECKS TO CODE BEFORE UPLOAD. Throw error if things don't go through
	logger.info("MERING, FILTERING, STORING CLINICAL FILES")
	### STORING CLINICAL FILES INTO CBIOPORTAL
	patientSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "patient"][0]
	sampleSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "sample"][0]
	
	#Using center mapping df to gate centers in release fileStage
	patient = syn.tableQuery("SELECT * FROM %s where CENTER in ('%s')" % (patientSynId,"','".join(CENTER_MAPPING_DF.center)))
	sample = syn.tableQuery("SELECT * FROM %s where CENTER in ('%s')" % (sampleSynId,"','".join(CENTER_MAPPING_DF.center)))
	patientDf = patient.asDataFrame()
	sampleDf = sample.asDataFrame()
	del sampleDf['AGE_AT_SEQ_REPORT_NUMERICAL']
	del sampleDf['CENTER']
	del patientDf['BIRTH_YEAR_NUMERICAL']
	#Clinical release scope filter
	#If private -> Don't release to public
	clinicalReleaseScope = syn.tableQuery("SELECT * FROM syn8545211 where releaseScope <> 'private'")
	clinicalReleaseScopeDf = clinicalReleaseScope.asDataFrame()

	patientCols = clinicalReleaseScopeDf['fieldName'][clinicalReleaseScopeDf['level'] == "patient"].tolist()
	sampleCols = clinicalReleaseScopeDf['fieldName'][clinicalReleaseScopeDf['level'] == "sample"].tolist()

	totalSample = ['PATIENT_ID']
	totalSample.extend(sampleCols)
	sampleCols = totalSample
	clinicalDf = sampleDf.merge(patientDf, on="PATIENT_ID",how="outer")
	#Remove patients without any sample or patient ids
	clinicalDf = clinicalDf[~clinicalDf['SAMPLE_ID'].isnull()]
	clinicalDf = clinicalDf[~clinicalDf['PATIENT_ID'].isnull()]

	#add in vital status
	logger.info("MERGE VITAL STATUS")
	vitalStatusSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "vitalStatus"][0]
	vitalStatus = syn.tableQuery('select * from %s' % vitalStatusSynId)
	vitalStatusDf = vitalStatus.asDataFrame()
	del vitalStatusDf['CENTER']
	#Make sure to grab only the patients that exist in the clinical database
	vitalStatusDf = vitalStatusDf[vitalStatusDf.PATIENT_ID.isin(clinicalDf.PATIENT_ID)]
	clinicalDf = clinicalDf.merge(vitalStatusDf, on = "PATIENT_ID",how="outer")

	#########FILTERING#########
	logger.info("REMOVING PHI")
	clinicalDf = reAnnotatePHI(clinicalDf)
	logger.info("MAF IN BED FILTER")
	remove_mafInBed_variants = runMAFinBED(syn, CENTER_MAPPING_DF, databaseSynIdMappingDf, test, genieVersion=genieVersion)
	logger.info("MUTATION IN CIS FILTER")
	remove_mutationInCis_samples = mutation_in_cis_filter(syn, skipMutationsInCis, test, variant_filtering_synId, CENTER_MAPPING_DF, genieVersion=genieVersion)
	logger.info("SEQ DATE FILTER")
	remove_seqDate_samples = seq_date_filter(clinicalDf, processingDate, consortiumReleaseCutOff)
	#Only certain samples are removed for the files that go into staging center folder
	removeForCenterConsortiumSamples = remove_mutationInCis_samples
	#Most filteres are applied for the files that go into the merged consortium release
	removeForMergedConsortiumSamples = set(remove_seqDate_samples)#set(remove_seqAssayId_samples)#.union(set(remove_seqDate_samples))
	removeForMergedConsortiumSamples = set(removeForMergedConsortiumSamples).union(set(removeForCenterConsortiumSamples))

	logger.info("ADD CANCER TYPES")
	#This removes support for both oncotree urls (only support json)
	oncotreeDict = process.get_oncotree_code_mappings(oncotree_url)
	#Add in unknown key which maps to UNKNOWN everything
	oncotree_mapping_dict['UNKNOWN']= {'CANCER_TYPE': 'UNKNOWN',
			  'CANCER_TYPE_DETAILED': 'UNKNOWN',
			  'ONCOTREE_PRIMARY_NODE': 'UNKNOWN',
			  'ONCOTREE_SECONDARY_NODE': 'UNKNOWN'}
	clinicalDf['CANCER_TYPE'] = [oncotreeDict[code.upper()].get("CANCER_TYPE",float('nan')) for code in clinicalDf['ONCOTREE_CODE']]
	clinicalDf['CANCER_TYPE_DETAILED'] = [oncotreeDict[code.upper()].get("CANCER_TYPE_DETAILED",float('nan')) for code in clinicalDf['ONCOTREE_CODE']]
	clinicalDf['ONCOTREE_PRIMARY_NODE'] = [oncotreeDict[code.upper()].get("ONCOTREE_PRIMARY_NODE",float('nan')) for code in clinicalDf['ONCOTREE_CODE']]
	clinicalDf['ONCOTREE_SECONDARY_NODE'] = [oncotreeDict[code.upper()].get("ONCOTREE_SECONDARY_NODE",float('nan')) for code in clinicalDf['ONCOTREE_CODE']]

	#CANCER TYPES are added which is why the clinical file is written out.
	#clinicalDf.to_csv(CLINCICAL_PATH, sep="\t", index=False)
	#add_cancerType_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../analyses/clinicalData/oncotree_code_converter.py')
	#command = ['python',add_cancerType_script,'-o',oncotree_url,'-c',CLINCICAL_PATH]
	#subprocess.check_call(command)
	#clinicalDf = pd.read_csv(CLINCICAL_PATH, sep="\t", comment="#")

	#All cancer types that are null should have null oncotree codes
	clinicalDf['ONCOTREE_CODE'][clinicalDf['CANCER_TYPE'].isnull()] = float('nan')
	# Suggest using AGE_AT_SEQ_REPORT_DAYS instead so that the descriptions can match
	clinicalDf['AGE_AT_SEQ_REPORT_DAYS'] = clinicalDf['AGE_AT_SEQ_REPORT']
	clinicalDf['AGE_AT_SEQ_REPORT'] = [int(math.floor(int(float(i))/365.25)) if process.checkInt(i) else i for i in clinicalDf['AGE_AT_SEQ_REPORT']]
	clinicalDf['AGE_AT_SEQ_REPORT'][clinicalDf['AGE_AT_SEQ_REPORT'] == ">32485"] = ">89"
	clinicalDf['AGE_AT_SEQ_REPORT'][clinicalDf['AGE_AT_SEQ_REPORT'] == "<6570"] = "<18"

	############################################################
	#CENTER SPECIFIC CODE FOR RIGHT NOW (REMOVE UHN-555-V1)
	############################################################
	clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'] != "UHN-555-V1"]
	#clinicalDf = clinicalDf[clinicalDf['CENTER'] != "WAKE"]
	#clinicalDf = clinicalDf[clinicalDf['CENTER'] != "CRUK"]
	############################################################
	############################################################

	clinicalDf.drop_duplicates("SAMPLE_ID",inplace=True)

	#samples must be removed after reading in the clinical file again

	clinicalDfStaging = clinicalDf[~clinicalDf['SAMPLE_ID'].isin(removeForCenterConsortiumSamples)]
	if not current_release_staging:
		for center in CENTER_MAPPING_DF.center:
			center_clinical =  clinicalDfStaging[clinicalDfStaging['CENTER'] == center]
			center_sample = center_clinical[sampleCols].drop_duplicates('SAMPLE_ID')
			center_patient = center_clinical[patientCols].drop_duplicates('PATIENT_ID')
			center_sample.to_csv(SAMPLE_CENTER_PATH % center, sep="\t",index=False)
			center_patient.to_csv(PATIENT_CENTER_PATH % center, sep="\t",index=False)
			storeFile(syn, SAMPLE_CENTER_PATH % center, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)
			storeFile(syn, PATIENT_CENTER_PATH % center, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)
	
	clinicalDf = clinicalDf[~clinicalDf['SAMPLE_ID'].isin(removeForMergedConsortiumSamples)]
	#This must happen here because the seq assay filter must happen after all the other filters
	# logger.info("SEQ ASSAY FILTER")
	# remove_seqAssayId_samples = seq_assay_id_filter(clinicalDf)
	#removeForMergedConsortiumSamples = removeForMergedConsortiumSamples.union(set(remove_seqAssayId_samples))
	# clinicalDf = clinicalDf[~clinicalDf['SAMPLE_ID'].isin(remove_seqAssayId_samples)]

	keepForCenterConsortiumSamples = clinicalDfStaging.SAMPLE_ID
	keepForMergedConsortiumSamples = clinicalDf.SAMPLE_ID
	#This mapping table is the GENIE clinical code to description mapping to generate the headers of the clinical file
	mapping_table = syn.tableQuery('SELECT * FROM syn9621600')
	mapping = mapping_table.asDataFrame()
	
	process.addClinicalHeaders(clinicalDf, mapping, patientCols, sampleCols, CLINCICAL_SAMPLE_PATH, CLINCICAL_PATIENT_PATH)
	storeFile(syn, CLINCICAL_SAMPLE_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_clinical_sample.txt", staging=current_release_staging)
	storeFile(syn, CLINCICAL_PATIENT_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_clinical_patient.txt", staging=current_release_staging)
	clinicalDf.to_csv(CLINCICAL_PATH, sep="\t", index=False)
	storeFile(syn, CLINCICAL_PATH, parent= consortiumReleaseSynId, name="data_clinical.txt", staging=current_release_staging)


	logger.info("CREATING DATA GENE PANEL FILE")
	#Samples have already been removed
	data_gene_panel = pd.DataFrame(columns = ["SAMPLE_ID","SEQ_ASSAY_ID"])
	data_gene_panel = data_gene_panel.append(clinicalDf[['SAMPLE_ID', 'SEQ_ASSAY_ID']])
	data_gene_panel = data_gene_panel.rename(columns={"SEQ_ASSAY_ID":"mutations"})
	data_gene_panel = data_gene_panel[data_gene_panel['SAMPLE_ID'] != ""]
	data_gene_panel.drop_duplicates("SAMPLE_ID",inplace=True)
	# Gene panel file is written below CNA, because of the "cna" column

	samples = data_gene_panel['SAMPLE_ID']
	sequenced_samples = "#sequenced_samples: " + " ".join(samples)

	logger.info("FILTERING, STORING MUTATION FILES")
	centerMafFileViewSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "centerMafView"][0]
	centerMafSynIds = syn.tableQuery("select id from %s " % centerMafFileViewSynId + "where name like '%mutation%'")
	centerMafSynIdsDf = centerMafSynIds.asDataFrame()
	with open(MUTATIONS_PATH, 'w') as f:
		f.write(sequenced_samples + "\n") 
	for index, mafSynId in enumerate(centerMafSynIdsDf.id):
		mafEnt = syn.get(mafSynId)
		logger.info(mafEnt.path)
		with open(mafEnt.path,"r") as mafFile:
			header = mafFile.readline()
			headers = header.replace("\n","").split("\t")
			if index == 0:
				with open(MUTATIONS_PATH, 'a') as f:
					f.write(header)
				#Create maf file per center for their staging directory
				for center in clinicalDf['CENTER'].unique():
					with open(MUTATIONS_CENTER_PATH % center, 'w') as f:
						f.write(header)
		# with open(mafEnt.path,"r") as newMafFile:
		# 	newMafFile.readline()
			center = mafEnt.path.split("_")[3]
			#Make sure to only write the centers that release = True
			if center in CENTER_MAPPING_DF.center.tolist():
				for row in mafFile:
					rowArray = row.replace("\n","").split("\t")
					center = rowArray[headers.index('Center')]
					newMergedRow = configureMafRow(rowArray, headers, keepForMergedConsortiumSamples, remove_mafInBed_variants)
					if newMergedRow is not None:
						with open(MUTATIONS_PATH, 'a') as f:
							f.write(newMergedRow)
					newCenterRow = configureMafRow(rowArray, headers, keepForCenterConsortiumSamples, remove_mafInBed_variants)
					if newCenterRow is not None:
						with open(MUTATIONS_CENTER_PATH % center, 'a') as f:
							f.write(newCenterRow)
	storeFile(syn, MUTATIONS_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_mutations_extended.txt", staging=current_release_staging)
	if not current_release_staging:
		for center in clinicalDf['CENTER'].unique():
			storeFile(syn, MUTATIONS_CENTER_PATH % center, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)

	#Only need to upload these files once
	#if filtering:
		#Normalize gene panels
	logger.info("STORING GENE PANELS FILES")
	fileviewSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "fileview"][0]
	genePanels = syn.tableQuery("select id from %s where cBioFileFormat = 'genePanel' and fileStage = 'staging'" % fileviewSynId)
	genePanelDf = genePanels.asDataFrame()
	genePanelEntities = []
	for synId in genePanelDf['id']:
		genePanel = syn.get(synId)
		panelNames = set(data_gene_panel['mutations'])
		genePanelName = os.path.basename(genePanel.path)
		newGenePanelPath = os.path.join(GENIE_RELEASE_DIR, genePanelName.replace(".txt","_%s.txt" % genieVersion))
		if genePanelName.replace(".txt","").replace("data_gene_panel_","") in panelNames:
			os.rename(genePanel.path, newGenePanelPath)
			genePanelEntities.append(storeFile(syn, newGenePanelPath, parent= consortiumReleaseSynId, genieVersion=genieVersion, name=genePanelName, cBioFileFormat = "genePanel", staging=current_release_staging))
	
	## CNA
	logger.info("MERING, FILTERING, STORING CNA FILES")
	centerCNASynIds = syn.tableQuery("select id from %s " % centerMafFileViewSynId + "where name like 'data_CNA%'")
	centerCNASynIdsDf = centerCNASynIds.asDataFrame()
	#Grab all unique symbols and form cnaTemplate
	allSymbols = set()

	for cnaSynId in centerCNASynIdsDf.id:
		cnaEnt = syn.get(cnaSynId)
		with open(cnaEnt.path,"r") as cnaFile:
			#Read first line first to get all the samples
			samples = cnaFile.readline()
			#Get all hugo symbols
			allSymbols = allSymbols.union(set(line.split("\t")[0] for line in cnaFile))
	cnaTemplate = pd.DataFrame({"Hugo_Symbol":list(allSymbols)})
	cnaTemplate.sort_values("Hugo_Symbol",inplace=True)
	cnaTemplate.to_csv(CNA_PATH, sep="\t",index=False)
	#Loop through to create finalized CNA file
	withCenterHugoSymbol = pd.Series("Hugo_Symbol")
	withCenterHugoSymbol = withCenterHugoSymbol.append(pd.Series(keepForCenterConsortiumSamples))

	withMergedHugoSymbol = pd.Series("Hugo_Symbol")
	withMergedHugoSymbol = withMergedHugoSymbol.append(pd.Series(keepForMergedConsortiumSamples))

	cnaSamples = []

	for cnaSynId in centerCNASynIdsDf.id:
		cnaEnt = syn.get(cnaSynId)
		center = cnaEnt.name.replace("data_CNA_","").replace(".txt","")
		logger.info(cnaEnt.path)
		if center in CENTER_MAPPING_DF.center.tolist():
			centerCNA = pd.read_csv(cnaEnt.path,sep="\t")
			merged = cnaTemplate.merge(centerCNA, on="Hugo_Symbol", how="outer")
			merged.sort_values("Hugo_Symbol",inplace=True)

			merged = merged[merged.columns[merged.columns.isin(withCenterHugoSymbol)]]

			cnaText = removePandasDfFloat(merged)
			#Replace blank with NA's
			cnaText = cnaText.replace("\t\t","\tNA\t").replace("\t\t","\tNA\t").replace('\t\n',"\tNA\n")

			#Store center CNA file in staging dir
			with open(CNA_CENTER_PATH % center, "w") as cnaFile:
				cnaFile.write(cnaText)
			storeFile(syn, CNA_CENTER_PATH % center, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)
			#This is to remove more samples for the final cna file
			merged = merged[merged.columns[merged.columns.isin(withMergedHugoSymbol)]]

			cnaText = removePandasDfFloat(merged)
			cnaText = cnaText.replace("\t\t","\tNA\t").replace("\t\t","\tNA\t").replace('\t\n',"\tNA\n")

			with open(CNA_CENTER_PATH % center, "w") as cnaFile:
				cnaFile.write(cnaText)
			#Join CNA file
			cnaSamples.extend(merged.columns[1:].tolist())
			joinCommand = ["join",CNA_PATH, CNA_CENTER_PATH % center]
			output = subprocess.check_output(joinCommand)
			with open(CNA_PATH, "w") as cnaFile:
				cnaFile.write(output.decode("utf-8").replace(" ","\t"))

	storeFile(syn, CNA_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_CNA.txt", staging=current_release_staging)

	# cnaText = removePandasDfFloat(mergedCNA)
	# with open(CNA_PATH, "w") as cnaFile:
	# 	cnaFile.write(cnaText)
	# storeFile(syn, CNA_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_CNA.txt", staging=current_release_staging)

	# mergedCNA = pd.DataFrame(columns=["Hugo_Symbol"],index=[])
	# cnaSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "cna"][0]

	# for center in CENTER_MAPPING_DF.center:
	# 	logger.info("MERGING %s CNA" % center)
	# 	cna = syn.tableQuery("SELECT TUMOR_SAMPLE_BARCODE,CNAData FROM %s where CENTER = '%s'" % (cnaSynId,center))
	# 	cnaDf = cna.asDataFrame()
	# 	if not cnaDf.empty:
	# 		allSymbols = set()
	# 		for data in cnaDf['CNAData']:
	# 			cnadata = data.split("\n")
	# 			symbols = cnadata[0].split(",")
	# 			allSymbols.update(set(symbols))
	# 		center_cna = pd.DataFrame(index = allSymbols)
	# 		if len(symbols) == len(allSymbols):
	# 			center_cna = center_cna.ix[symbols]
	# 		for sample,data in zip(cnaDf['TUMOR_SAMPLE_BARCODE'],cnaDf['CNAData']):
	# 			cnadata = data.split("\n")
	# 			symbols = cnadata[0].split(",")
	# 			cnavalues = cnadata[1].split(",")
	# 			indexMethod = True
	# 			if len(center_cna.index) == len(symbols):
	# 				if all(center_cna.index == symbols):
	# 					center_cna[sample] = cnavalues
	# 					indexMethod = False
	# 			if indexMethod:
	# 				center_cna[sample] = 0
	# 				center_cna[sample].ix[symbols] = cnavalues
	# 		#missing = center_cna.columns[~center_cna.columns.isin(samples)]
	# 		center_cna = center_cna[center_cna.columns[center_cna.columns.isin(keepForCenterConsortiumSamples)]]
	# 		center_cna['Hugo_Symbol'] = center_cna.index
	# 		center_cna = center_cna.fillna('NA')
	# 		cols = center_cna.columns.tolist()
	# 		cols = cols[-1:] + cols[:-1]
	# 		center_cna = center_cna[cols]
	# 		cnaText = removePandasDfFloat(center_cna)
	# 		with open(CNA_CENTER_PATH % center, "w") as cnaFile:
	# 			cnaFile.write(cnaText)
	# 		if not current_release_staging:
	# 			storeFile(syn,CNA_CENTER_PATH % center, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)
	# 		mergedCNA = mergedCNA.merge(center_cna, on='Hugo_Symbol', how="outer")

	# mergedCNA = mergedCNA[mergedCNA.columns[mergedCNA.columns.isin(keepForMergedConsortiumSamples.append(pd.Series("Hugo_Symbol")))]]
	# mergedCNA = mergedCNA.fillna('NA')
	# cnaText = removePandasDfFloat(mergedCNA)
	# with open(CNA_PATH, "w") as cnaFile:
	# 	cnaFile.write(cnaText)
	# storeFile(syn, CNA_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_CNA.txt", staging=current_release_staging)

	# Add in CNA column into gene panel file
	cnaSeqIds = data_gene_panel['mutations'][data_gene_panel['SAMPLE_ID'].isin(cnaSamples)].unique()
	data_gene_panel['cna'] = data_gene_panel['mutations']
	data_gene_panel['cna'][~data_gene_panel['cna'].isin(cnaSeqIds)] = "NA"
	data_gene_panel.to_csv(DATA_GENE_PANEL_PATH,sep="\t",index=False)
	storeFile(syn, DATA_GENE_PANEL_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_gene_matrix.txt", staging=current_release_staging)

	### FUSIONS
	logger.info("MERING, FILTERING, STORING FUSION FILES")
	fusionSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "fusions"][0]

	Fusions = syn.tableQuery('SELECT HUGO_SYMBOL,ENTREZ_GENE_ID,CENTER,TUMOR_SAMPLE_BARCODE,FUSION,DNA_SUPPORT,RNA_SUPPORT,METHOD,FRAME FROM %s' % fusionSynId)
	FusionsDf = Fusions.asDataFrame()
	#FusionsDf['ENTREZ_GENE_ID'] = FusionsDf['ENTREZ_GENE_ID'].fillna(0).apply(int)
	FusionsDf['ENTREZ_GENE_ID'][FusionsDf['ENTREZ_GENE_ID'] == 0] = pd.np.nan
	# missing = FusionsDf[~FusionsDf['TUMOR_SAMPLE_BARCODE'].isin(samples)]['TUMOR_SAMPLE_BARCODE']
	# missing.drop_duplicates().to_csv("fusion_missing.txt",sep="\t",index=False)
	if not current_release_staging:
		FusionsStagingDf = FusionsDf[FusionsDf['TUMOR_SAMPLE_BARCODE'].isin(keepForCenterConsortiumSamples)]
		for center in CENTER_MAPPING_DF.center:
			center_fusion = FusionsStagingDf[FusionsStagingDf['CENTER'] == center]
			if not center_fusion.empty:
				center_fusion.to_csv(FUSIONS_CENTER_PATH % center, sep="\t",index=False)
				storeFile(syn, FUSIONS_CENTER_PATH % center, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)
		
	FusionsDf = FusionsDf[FusionsDf['TUMOR_SAMPLE_BARCODE'].isin(keepForMergedConsortiumSamples)]
	FusionsDf = FusionsDf.rename(columns= {'HUGO_SYMBOL':'Hugo_Symbol','ENTREZ_GENE_ID':'Entrez_Gene_Id','CENTER':'Center','TUMOR_SAMPLE_BARCODE':'Tumor_Sample_Barcode','FUSION':'Fusion',
										   'DNA_SUPPORT':'DNA_support','RNA_SUPPORT':'RNA_support','METHOD':'Method','FRAME':'Frame'})
	#Remove duplicated Fusions
	FusionsDf = FusionsDf[~FusionsDf[['Hugo_Symbol','Tumor_Sample_Barcode','Fusion']].duplicated()]
	#FusionsDf.to_csv(FUSIONS_PATH, sep="\t", index=False)
	fusionText = removePandasDfFloat(FusionsDf)
	with open(FUSIONS_PATH, "w") as fusionFile:
		fusionFile.write(fusionText)
	storeFile(syn, FUSIONS_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="data_fusions.txt", staging=current_release_staging)

	### SEG
	logger.info("MERING, FILTERING, STORING SEG FILES")
	segSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "seg"][0]
	seg = syn.tableQuery('SELECT ID,CHROM,LOCSTART,LOCEND,NUMMARK,SEGMEAN,CENTER FROM %s' % segSynId)
	segDf = seg.asDataFrame()
	segDf = segDf.rename(columns= {'CHROM':'chrom','LOCSTART':'loc.start','LOCEND':'loc.end','SEGMEAN':'seg.mean','NUMMARK':'num.mark'})
	# missing = segDf[~segDf['ID'].isin(samples)]['ID']
	# missing.drop_duplicates().to_csv("seg_missing.txt",sep="\t",index=False)
	if not current_release_staging:
		segStagingDf = segDf[segDf['ID'].isin(keepForCenterConsortiumSamples)]
		for center in CENTER_MAPPING_DF.center:
			center_seg = segStagingDf[segStagingDf['CENTER'] == center]
			if not center_seg.empty:
				del center_seg['CENTER']
				segText = removePandasDfFloat(center_seg)
				with open(SEG_CENTER_PATH % center, "w") as segFile:
					segFile.write(segText)
				storeFile(syn, SEG_CENTER_PATH % center, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)
	del segDf['CENTER']
	segDf = segDf[segDf['ID'].isin(keepForMergedConsortiumSamples)]
	segText = removePandasDfFloat(segDf)
	with open(SEG_PATH, "w") as segFile:
		segFile.write(segText)
	storeFile(syn, SEG_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="genie_private_data_cna_hg19.seg", staging=current_release_staging)

	#BED
	logger.info("STORING COMBINED BED FILE")
	bedSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == "bed"][0]
	bed = syn.tableQuery("SELECT Chromosome,Start_Position,End_Position,Hugo_Symbol,ID,SEQ_ASSAY_ID,Feature_Type,includeInPanel FROM %s where CENTER in ('%s')" % (bedSynId,"','".join(CENTER_MAPPING_DF.center)))
	bedDf = bed.asDataFrame()
	if not current_release_staging:
		for seqAssay in bedDf['SEQ_ASSAY_ID'].unique():
			bedSeqDf = bedDf[bedDf['SEQ_ASSAY_ID'] == seqAssay]
			center = seqAssay.split("-")[0]
			bedSeqDf = bedSeqDf[bedSeqDf['Hugo_Symbol'] != bedSeqDf['ID']]
			if not bedSeqDf.empty:
				bedSeqDf.to_csv(BED_DIFFS_SEQASSAY_PATH % seqAssay, index=False)
				storeFile(syn, BED_DIFFS_SEQASSAY_PATH % seqAssay, genieVersion=genieVersion, parent = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0], centerStaging=True)
	#This clinicalDf is already filtered through most of the filters
	bedDf = bedDf[bedDf['SEQ_ASSAY_ID'].isin(clinicalDf.SEQ_ASSAY_ID)]
	bedDf.to_csv(COMBINED_BED_PATH, sep="\t",index=False)
	storeFile(syn, COMBINED_BED_PATH, parent= consortiumReleaseSynId, genieVersion=genieVersion, name="genie_combined.bed", staging=current_release_staging)

	return(genePanelEntities)

def command_reviseMetadataFiles(syn, args, databaseSynIdMappingDf):
	reviseMetadataFiles(syn, args.staging, databaseSynIdMappingDf, args.genieVersion)

def reviseMetadataFiles(syn, staging, databaseSynIdMappingDf, genieVersion=None):
	# if staging:
	# 	consortiumId = "syn7551278"
	# 	#allFiles = syn.getChildren('syn7551278')
	# else:
	consortiumId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'consortium'].values[0]
	allFiles = syn.getChildren(consortiumId)
	metadataEnts = [syn.get(i['id'], downloadLocation=GENIE_RELEASE_DIR, ifcollision="overwrite.local") for i in allFiles if 'meta' in i['name']]
	for metaEnt in metadataEnts:
		with open(metaEnt.path, "r+") as meta:
			metaText = meta.read()
			if "meta_study" not in metaEnt.path:
				version = ''
			else:
				version = re.search(".+GENIE.+v(.+)", metaText).group(1)
			#Fix this line
			genieVersion = version if genieVersion is None else genieVersion
			dataFileVersion = re.search(".+data_(.+)[.]txt",metaText)
			if dataFileVersion is None:
				dataFileVersion = re.search(".+data_(.+)[.]seg",metaText)
			if dataFileVersion is not None:
				dataFileVersion = dataFileVersion.group(1).split("_")[-1]

			if version != genieVersion:
				metaText = metaText.replace("AACR Project GENIE Cohort v%s" % version,"AACR Project GENIE Cohort v%s" % genieVersion)
				metaText = metaText.replace("GENIE v%s" % version,"GENIE v%s" % genieVersion)
				if dataFileVersion is not None:
					metaText = metaText.replace(dataFileVersion, genieVersion)
					metaText = metaText.replace(dataFileVersion, genieVersion)
				meta.seek(0)
				meta.write(metaText)
				meta.truncate()
				storeFile(syn, metaEnt.path, parent=consortiumId, genieVersion=genieVersion, staging=staging)

def createLinkVersion(syn,genieVersion, caseListEntities, genePanelEntities, databaseSynIdMappingDf):
	versioning = genieVersion.split(".")
	main = versioning[0]
	#second = ".".join(versioning[1:])
	releaseSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'release'].values[0]
	releases = synu.walk(syn, releaseSynId)
	mainReleaseFolders = next(releases)[1]
	releaseFolderSynId = [synId for folderName, synId in mainReleaseFolders if folderName == "Release %s" % main] 
	if len(releaseFolderSynId) > 0:
		secondRelease = synu.walk(syn, releaseFolderSynId[0])
		secondReleaseFolders = next(secondRelease)[1]
		secondReleaseFolderSynIdList = [synId for folderName, synId in secondReleaseFolders if folderName == genieVersion] 
		if len(secondReleaseFolderSynIdList) > 0:
			secondReleaseFolderSynId = secondReleaseFolderSynIdList[0]
		else:
			secondReleaseFolderSynId = syn.store(synapseclient.Folder(genieVersion, parent = releaseFolderSynId[0]))
	else:
		mainReleaseFolderId = syn.store(synapseclient.Folder("Release %s" % main, parent = releaseSynId)).id
		secondReleaseFolderSynId = syn.store(synapseclient.Folder(genieVersion, parent = mainReleaseFolderId)).id

	caselistId=findCaseListId(syn, secondReleaseFolderSynId)
	consortiumSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'consortium'].values[0]
	consortiumRelease = syn.getChildren(consortiumSynId)
	#[syn.store(synapseclient.Link(ents['id'], parent=secondReleaseFolderSynId, targetVersion=ents['versionNumber'])) for ents in publicRelease if ents['type'] != "org.sagebionetworks.repo.model.Folder" and ents['name'] != "data_clinical.txt"]
	#data_clinical.txt MUST be pulled in because the clinical file is needed in the consortium_to_public.py
	[syn.store(synapseclient.Link(ents['id'], parent=secondReleaseFolderSynId, targetVersion=ents['versionNumber'])) for ents in consortiumRelease if ents['type'] != "org.sagebionetworks.repo.model.Folder" and not ents['name'].startswith("data_gene_panel")]
	releaseFiles = syn.getChildren(secondReleaseFolderSynId)
	clinEnt = [ents['id'] for ents in releaseFiles if ents['name'] == "data_clinical.txt"][0]
	#Set private permission for the data_clinical.txt Link
	syn.setPermissions(clinEnt, principalId=3346558, accessType=[])
	syn.setPermissions(clinEnt, principalId=3326313, accessType=[])
	#caselistEnts = syn.getChildren("syn9689663")
	[syn.store(synapseclient.Link(ents.id, parent=caselistId, targetVersion=ents.versionNumber)) for ents in caseListEntities]
	#Store gene panels
	[syn.store(synapseclient.Link(ents.id, parent=secondReleaseFolderSynId, targetVersion=ents.versionNumber)) for ents in genePanelEntities]

def main():
	parser = argparse.ArgumentParser(description='Staging to CBioOutput')
	parser.add_argument("processingDate", type=str, metavar="Jan-2017",
						help="The processing date of GENIE in Month-Year format (ie. Apr-2017)")
	parser.add_argument("cbioportalPath", type=str, metavar="/path/to/cbioportal",
						help="Make sure you clone the cbioportal github: git clone https://github.com/cBioPortal/cbioportal.git")
	parser.add_argument("genieVersion", type=str, metavar="1.0.1",
						help="GENIE release version")
	parser.add_argument("--oncotreeLink", type=str,
						help="Link to oncotree code")
	parser.add_argument("--consortiumReleaseCutOff", type=int, metavar=184, default=184,
						help="Consortium release cut off time in days")
	parser.add_argument("--staging", action='store_true',
						help="Store into staging folder")
	parser.add_argument("--skipMutationsInCis", action='store_true',
						help="Skip running mutation in cis script")
	parser.add_argument("--test", action='store_true',
						help="Run test")
	parser.add_argument("--pemFile", type=str, 
						help="Path to PEM file (genie.pem)")
	parser.add_argument("--debug", action='store_true', 
						help="Synapse debug feature")
	args = parser.parse_args()
	syn = process.synLogin(args.pemFile, debug=args.debug)

	assert not (args.test and args.staging), "You can only specify --test or --staging, not both"

	if args.test:
		databaseSynIdMappingId = 'syn11600968'
		args.genieVersion = "TESTING"
	elif args.staging:
		args.skipMutationsInCis = True
		databaseSynIdMappingId = 'syn12094210'	
	else:
		databaseSynIdMappingId = 'syn10967259'
	#Database/folder syn id mapping
	databaseSynIdMapping = syn.tableQuery('select * from %s' % databaseSynIdMappingId)
	databaseSynIdMappingDf = databaseSynIdMapping.asDataFrame()
	consortiumSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'consortium'].values[0]

	if args.oncotreeLink is None:
		oncoLink = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
		oncoLinkEnt = syn.get(oncoLink)
		args.oncotreeLink = oncoLinkEnt.externalURL

	#Check if you can connect to oncotree link, if not then don't run validation / processing
	process.checkUrl(args.oncotreeLink)

	#get syn id of case list folder in consortium release
	caseListSynId = findCaseListId(syn, consortiumSynId)

	if not args.test and not args.staging:
		processTrackerSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'processTracker'].values[0]
		processTracker = syn.tableQuery("SELECT timeStartProcessing FROM %s where center = 'SAGE' and processingType = 'dbToStage'" % processTrackerSynId)
		processTrackerDf = processTracker.asDataFrame()
		processTrackerDf['timeStartProcessing'][0] = str(int(time.time()*1000))
		syn.store(synapseclient.Table(processTrackerSynId,processTrackerDf))

	syn.table_query_timeout = 50000
	centerMappingSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'centerMapping'].values[0]
	#Only release files where release is true
	CENTER_MAPPING = syn.tableQuery('SELECT * FROM %s where release is true' % centerMappingSynId)
	CENTER_MAPPING_DF = CENTER_MAPPING.asDataFrame()
	processingDate = datetime.datetime.strptime(args.processingDate, '%b-%Y')
#######
	cbioValidatorPath = os.path.join(args.cbioportalPath,"core/src/main/scripts/importer/validateData.py")
	assert os.path.exists(cbioValidatorPath), "Please specify correct cbioportalPath"

	logger.info("STAGING TO CONSORTIUM")
	genePanelEntities = stagingToCbio(syn, processingDate, args.genieVersion, CENTER_MAPPING_DF, databaseSynIdMappingDf, oncotree_url=args.oncotreeLink, consortiumReleaseCutOff= args.consortiumReleaseCutOff, current_release_staging=args.staging, skipMutationsInCis=args.skipMutationsInCis, test=args.test)
	
	#No need to run twice anymore
	#stagingToCbio(syn, processingDate, args.genieVersion, CENTER_MAPPING_DF, databaseSynIdMappingDf, oncotree_url=args.oncotreeLink, consortiumReleaseCutOff= args.consortiumReleaseCutOff, filtering=True, current_release_staging=args.staging, test=args.test)
	#Create case lists files
	logger.info("CREATE CASE LIST FILES")
	#Remove old caselists first
	if not os.path.exists(CASE_LIST_PATH):
		os.mkdir(CASE_LIST_PATH)
	caselists = os.listdir(CASE_LIST_PATH)
	[os.remove(os.path.join(CASE_LIST_PATH,caselist)) for caselist in caselists]
	CLINICAL_PATH = os.path.join(GENIE_RELEASE_DIR,'data_clinical_%s.txt' % args.genieVersion)
	GENE_MATRIX_PATH = os.path.join(GENIE_RELEASE_DIR,"data_gene_matrix_%s.txt" % args.genieVersion)
	create_case_lists.create_case_lists(CLINICAL_PATH, GENE_MATRIX_PATH, CASE_LIST_PATH, "genie_private")
	caseListFiles = os.listdir(CASE_LIST_PATH)
	caseListEntities = []
	for casePath in caseListFiles:
		casePath = os.path.join(CASE_LIST_PATH, casePath)
		caseListEntities.append(storeFile(syn, casePath, parent=caseListSynId, staging=args.staging, caseLists=True, genieVersion=args.genieVersion))

	logger.info("REMOVING UNNECESSARY FILES")
	genie_files = os.listdir(GENIE_RELEASE_DIR)
	#deletePatterns = ('data_clinical_supp_patient_','data_clinical_supp_sample_','data_CNA_','data_mutations_extended_','data_fusions_','genie_private_data_cna_hg19_')
	#[os.remove(os.path.join(GENIE_RELEASE_DIR,genieFile)) for genieFile in genie_files if genieFile.startswith(deletePatterns)]
	[os.remove(os.path.join(GENIE_RELEASE_DIR,genieFile)) for genieFile in genie_files if args.genieVersion not in genieFile and "meta" not in genieFile and "case_lists" not in genieFile]
	os.remove(CLINICAL_PATH)
#######	
	logger.info("REVISE METADATA FILES")
	command_reviseMetadataFiles(syn, args, databaseSynIdMappingDf)
	logger.info("CBIO VALIDATION")
	#Must be exit 0 because the validator sometimes fails, but we still want to capture the output	
	command = [cbioValidatorPath,'-s',GENIE_RELEASE_DIR,'-n','; exit 0']
	cbioOutput = subprocess.check_output(" ".join(command), shell=True)
	logger.info(cbioOutput.decode("utf-8"))
	if not args.test and not args.staging:
		with open("cbioValidatorLogsConsortium_%s.txt" % args.genieVersion, "w") as cbioLog:
			cbioLog.write(cbioOutput.decode("utf-8"))
		syn.store(File("cbioValidatorLogsConsortium_%s.txt" % args.genieVersion, parentId = "syn10155804"))
		os.remove("cbioValidatorLogsConsortium_%s.txt" % args.genieVersion)
	logger.info("REMOVING OLD FILES")

	process.rmFiles(CASE_LIST_PATH)
	if os.path.exists('%s/genie_private_meta_cna_hg19_seg.txt' % GENIE_RELEASE_DIR):
		os.unlink('%s/genie_private_meta_cna_hg19_seg.txt' % GENIE_RELEASE_DIR)
	logger.info("CREATING LINK VERSION")
	createLinkVersion(syn, args.genieVersion, caseListEntities, genePanelEntities, databaseSynIdMappingDf)

	if not args.test and not args.staging:
		processTracker = syn.tableQuery("SELECT timeEndProcessing FROM %s where center = 'SAGE' and processingType = 'dbToStage'" % processTrackerSynId)
		processTrackerDf = processTracker.asDataFrame()
		processTrackerDf['timeEndProcessing'][0] = str(int(time.time()*1000))
		syn.store(synapseclient.Table(processTrackerSynId,processTrackerDf))
	logger.info("COMPLETED DATABASE TO STAGING")

	if not args.test:
		logger.info("DASHBOARD UPDATE")
		dashboard_table_updater.run_dashboard(syn, databaseSynIdMappingDf, args.genieVersion, staging=args.staging)
		dashboard_markdown_html_commands = ['Rscript', os.path.join(os.path.dirname(os.path.abspath(__file__)),'dashboard_markdown_generator.R'), args.genieVersion]
		if args.staging:
			dashboard_markdown_html_commands.append('--staging')
		subprocess.check_call(dashboard_markdown_html_commands)
		logger.info("DASHBOARD UPDATE COMPLETE")

if __name__ == "__main__":
	main()

