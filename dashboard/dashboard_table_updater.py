from genie import process_functions as process
import sys
import os
import synapseclient
import synapseutils as synu
import pandas as pd
import argparse
import re
import datetime

def getCenterDataCompletion(center, df):
	centerDf = df[df['CENTER'] == center]
	total = len(centerDf)
	centerData = pd.DataFrame()
	for col in ['AGE_AT_SEQ_REPORT','BIRTH_YEAR', 'SEQ_ASSAY_ID']:
		if centerDf.get(col) is not None:
			completeness = float(sum(~centerDf[col].isnull())) / int(total)
			returned = pd.DataFrame([[col, center, total, completeness]])
			centerData = centerData.append(returned)
	for col in ['SAMPLE_TYPE','PRIMARY_RACE', 'SECONDARY_RACE', 'TERTIARY_RACE', 'SEX', 'ETHNICITY','VITAL_STATUS']:
		if centerDf.get(col) is not None:

			completeness = float(sum(centerDf[col] != "Unknown")) / int(total)
			returned = pd.DataFrame([[col, center, total, completeness]])
			centerData = centerData.append(returned)
	return(centerData)


def update_samples_in_release_table(syn, files, release, samples_in_release_synid):
	'''
	Convenience function that updates the sample in release table
	This tracks the samples of each release.  1 means it exists, and 0 means it doesn't

	params: 
		syn: synapse object
		files: file mapping generated from file mapping function
		release:  GENIE release number (ie. 5.3-consortium)
		samples_in_release_synid: Synapse Id of 'samples in release' Table
	'''
	clinical_ent = syn.get(files['clinical'],followLink=True)
	clinicaldf = pd.read_csv(clinical_ent.path,sep="\t",comment="#")
	cols = [i['name'] for i in list(syn.getTableColumns(samples_in_release_synid))]

	if release not in cols:
		schema = syn.get(samples_in_release_synid)
		new_column = syn.store(synapseclient.Column(name=release, columnType='INTEGER', defaultValue=0))
		schema.addColumn(new_column)
		schema = syn.store(schema)
	#Columns of samples in release
	samples_per_release = syn.tableQuery('SELECT SAMPLE_ID, "%s" FROM %s' % (release, samples_in_release_synid))
	samples_per_releasedf = samples_per_release.asDataFrame()
	new_samples = clinicaldf[['SAMPLE_ID']][~clinicaldf.SAMPLE_ID.isin(samples_per_releasedf.SAMPLE_ID)]
	new_samples[release] = 1
	old_samples = clinicaldf[['SAMPLE_ID']][clinicaldf.SAMPLE_ID.isin(samples_per_releasedf.SAMPLE_ID)]
	old_samples[release] = 1
	samples_in_releasedf = new_samples.append(old_samples)
	print(samples_in_releasedf)
	#process.updateDatabase(syn, samples_per_releasedf, samples_in_releasedf, samples_in_release_synid, ["SAMPLE_ID"])


def update_cumalative_sample_table(syn, files, release, cumulativeSampleCountSynId):
	### CONSORTIUM RELEASE SAMPLE COUNT TABLE UPDATE
	sampleCountPerRound = syn.tableQuery('SELECT * FROM %s' % cumulativeSampleCountSynId)
	sampleCountPerRoundDf = sampleCountPerRound.asDataFrame()

	clinical = syn.get(files['clinical'],followLink=True)
	clinicalDf = pd.read_csv(clinical.path,sep="\t",comment="#")
	clinicalDf.columns = [i.upper() for i in clinicalDf.columns]
	if clinicalDf.get("CENTER") is None:
		clinicalDf['CENTER'] = [sample.split("-")[1] for sample in clinicalDf.SAMPLE_ID]
	clinCounts = clinicalDf['CENTER'].value_counts()
	clinCounts['Total'] = sum(clinCounts)
	clinCounts.name = "Clinical"

	if release not in cols:
		schema = syn.get(samplesInReleaseSynId)
		new_column = syn.store(synapseclient.Column(name=release, columnType='INTEGER', defaultValue=0))
		schema.addColumn(new_column)
		schema = syn.store(schema)
	#Columns of samples in release
	print(samplesInReleaseSynId)
	samplesPerRelease = syn.tableQuery('SELECT SAMPLE_ID, "%s" FROM %s' % (release, samplesInReleaseSynId))
	samplesPerReleaseDf = samplesPerRelease.asDataFrame()
	newSamples = clinicalDf[['SAMPLE_ID']][~clinicalDf.SAMPLE_ID.isin(samplesPerReleaseDf.SAMPLE_ID)]
	newSamples[release] = 1
	oldSamples = clinicalDf[['SAMPLE_ID']][clinicalDf.SAMPLE_ID.isin(samplesPerReleaseDf.SAMPLE_ID)]
	oldSamples[release] = 1
	samplesInReleaseDf = newSamples.append(oldSamples)

	#process.updateDatabase(syn, samplesPerReleaseDf, samplesInReleaseDf, samplesInReleaseSynId, ["SAMPLE_ID"])

	fusion = syn.get(files['fusion'],followLink=True)
	fusionDf = pd.read_csv(fusion.path, sep="\t",comment="#")
	fusionDf.columns = [i.upper() for i in fusionDf.columns]

	fusionCounts = fusionDf['CENTER'][~fusionDf['TUMOR_SAMPLE_BARCODE'].duplicated()].value_counts()
	fusionCounts['Total'] = sum(fusionCounts)

	cna = syn.get(files['cna'],followLink=True)
	cnaDf = pd.read_csv(cna.path, sep="\t", comment="#")
	cnaCounts = pd.Series([i.split("-")[1] for i in cnaDf.columns[1:]]).value_counts()
	cnaCounts['Total'] = sum(cnaCounts)

	seg = syn.get(files['seg'],followLink=True)
	segDf = pd.read_csv(seg.path, sep="\t", comment="#")
	segDf.columns = [i.upper() for i in segDf.columns]

	segDf['CENTER'] = [i.split("-")[1] for i in segDf['ID']]
	segCounts = segDf['CENTER'][~segDf['ID'].duplicated()].value_counts()
	segCounts['Total'] = sum(segCounts)

	total = pd.DataFrame(clinCounts)
	total['Fusions'] = fusionCounts
	total['CNV'] = cnaCounts
	total['Mutation'] = clinCounts
	total['SEG'] = segCounts
	total = total.fillna(0)
	total = total.applymap(int)
	total['Center'] =  total.index

	total['Release'] = release
	print(total)
	#process.updateDatabase(syn, sampleCountPerRoundDf, total, cumulativeSampleCountSynId, ["Center", "Release"])

# def getFileDict(filenames):
# 	files=dict()
# 	for filename, synId in filenames:
# 		if not filename.startswith("meta"):
# 			if filename.startswith("data_clinical_sample"):
# 				files['clinical'] = synId
# 			elif filename.endswith("fusions.txt"):
# 				files['fusion'] = synId
# 			elif filename.endswith("CNA.txt"):
# 				files['cna'] = synId
# 			elif filename.endswith(".seg"):
# 				files['seg'] = synId
# 	return(files)

def get_file_mapping(syn, release_folder_synid):
	"""
	Get file mapping between important files needed for dashboard and 
	their synapse ids

	params:
		syn:  synapse object
		release_folder_synid: synapse id of release

	"""
	files = syn.getChildren(release_folder_synid)
	file_mapping = dict()
	for metadata in files:
		filename = metadata['name']
		synid = metadata['id']
		if not filename.startswith("meta"):
			if filename.startswith("data_clinical_sample"):
				file_mapping['clinical'] = synid
			elif filename.endswith("fusions.txt"):
				file_mapping['fusion'] = synid
			elif filename.endswith("CNA.txt"):
				file_mapping['cna'] = synid
			elif filename.endswith(".seg"):
				file_mapping['seg'] = synid
	return(file_mapping)

def update_release_numbers(syn, database_mappingdf, release = None):
	#Update release table with current release or all releases
	samples_in_release_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'samplesInRelease'].values[0]
	cumulative_sample_count_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'cumulativeSampleCount'].values[0]

	release_folder_fileview_synid = "syn17019650"
	release_folder = syn.tableQuery("select id,name from %s" % release_folder_fileview_synid + " where name not like 'Release%' and name <> 'case_lists'")
	release_folderdf = release_folder.asDataFrame()

	for release_synid, release_name in zip(release_folderdf.id, release_folderdf.name):
		file_mapping = get_file_mapping(syn, release_synid)
		if release is None:
			update_samples_in_release_table(syn, file_mapping, release_name, samples_in_release_synid)
			update_cumalative_sample_table(syn, file_mapping, release_name, cols, cumulative_sample_count_synid, samples_in_release_synid)
		elif release == release_name:
			update_samples_in_release_table(syn, file_mapping, release_name, samples_in_release_synid)
			update_cumalative_sample_table(syn, file_mapping, release_name, cols, cumulative_sample_count_synid, samples_in_release_synid)
		else:
			pass

# def update_database_numbers():
# 	#Update release table with releases
# 	releaseSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'release'].values[0]
# 	samplesInReleaseSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'samplesInRelease'].values[0]
# 	cumulativeSampleCountSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'cumulativeSampleCount'].values[0]
# 	oncotreeSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'oncotree'].values[0]
# 	primaryCodeSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'primaryCode'].values[0]
# 	sampleDiffCountSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'sampleDiffCount'].values[0]
# 	dataCompletionSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'dataCompletion'].values[0]

# 	# ## Database 
# 	sampleCountPerRound = syn.tableQuery('SELECT * FROM %s' % cumulativeSampleCountSynId)
# 	sampleCountPerRoundDf = sampleCountPerRound.asDataFrame()
# 	clinical = syn.tableQuery('select * from syn7517674')
# 	clinicalDf = clinical.asDataFrame()
# 	clinCounts = clinicalDf['CENTER'].value_counts()
# 	clinCounts['Total'] = sum(clinCounts)
# 	clinCounts.name = "Clinical"

# 	fusion = syn.tableQuery('select * from syn7893268')
# 	fusionDf = fusion.asDataFrame()
# 	fusionCounts = fusionDf['CENTER'][~fusionDf['TUMOR_SAMPLE_BARCODE'].duplicated()].value_counts()
# 	fusionCounts['Total'] = sum(fusionCounts)

# 	cna = syn.tableQuery('select * from syn8356049')
# 	cnaDf = cna.asDataFrame()
# 	cnaCounts = cnaDf['CENTER'].value_counts()
# 	cnaCounts['Total'] = sum(cnaCounts)

# 	seg = syn.tableQuery('select * from syn7893341')
# 	segDf = seg.asDataFrame()
# 	segCounts = segDf['CENTER'][~segDf['ID'].duplicated()].value_counts()
# 	segCounts['Total'] = sum(segCounts)

# 	total = pd.DataFrame(clinCounts)
# 	total['Fusions'] = fusionCounts
# 	total['CNV'] = cnaCounts
# 	total['Mutation'] = clinCounts
# 	total['SEG'] = segCounts
# 	total = total.fillna(0)
# 	total = total.applymap(int)
# 	total['Center'] =  total.index
# 	total['Release'] = "Database"

# 	process.updateDatabase(syn, sampleCountPerRoundDf, total, cumulativeSampleCountSynId, ["Center", "Release"])

# 	# ### DISTRIBUTION OF ONCOTREE CODE TABLE UPDATE
# 	oncotreeCodeDistributionDf = pd.DataFrame(columns=set(clinicalDf['CENTER']), index=set(clinicalDf['ONCOTREE_CODE']))
# 	for center in oncotreeCodeDistributionDf.columns:
# 		onc_counts = clinicalDf['ONCOTREE_CODE'][clinicalDf['CENTER'] == center].value_counts()
# 		oncotreeCodeDistributionDf[center] = onc_counts
# 	oncotreeCodeDistributionDf = oncotreeCodeDistributionDf.fillna(0)
# 	oncotreeCodeDistributionDf = oncotreeCodeDistributionDf.applymap(int)
# 	oncotreeCodeDistributionDf['Total'] = oncotreeCodeDistributionDf.apply(sum, axis=1)
# 	oncotreeCodeDistributionDf['Oncotree_Code'] = oncotreeCodeDistributionDf.index

# 	oncotreeDatabase = syn.tableQuery('SELECT %s FROM %s' % ("Oncotree_Code," + ",".join(clinicalDf['CENTER'].unique()) + ",Total", oncotreeSynId))
# 	oncotreeDatabaseDf = oncotreeDatabase.asDataFrame()
# 	process.updateDatabase(syn, oncotreeDatabaseDf, oncotreeCodeDistributionDf, oncotreeSynId, ["Oncotree_Code"],toDelete=True)
	
# 	# ### DISTRIBUTION OF PRIMARY CODE TABLE UPDATE

# 	ONCOTREE_MAP = pd.read_csv(oncotreeLink,sep="\t")
# 	if not ONCOTREE_MAP.empty:
# 		levels = [col for col in ONCOTREE_MAP.columns if "level_" in col]
# 		oncotreeDict = {}
# 		for level in levels:
# 			oncotreeDict[level] = pd.Series([re.sub(".+[(](.+)[)]","\\1",code) if not pd.isnull(code) else '' for code in ONCOTREE_MAP[level]])

# 		primary = oncotreeDict.pop('level_1')
# 		#Need to optimize this at some point
# 		clinicalDf['PRIMARY_CODES'] = clinicalDf.ONCOTREE_CODE.apply(lambda code: process.getPrimary(code, oncotreeDict, primary))
# 	else:
# 		ONCOTREE_MAP = process.get_oncotree_code_mappings(oncotreeLink)
# 		clinicalDf['PRIMARY_CODES'] = [ONCOTREE_MAP[i.upper()]['ONCOTREE_PRIMARY_NODE'] for i in clinicalDf.ONCOTREE_CODE]

# 	# ### DISTRIBUTION OF PRIMARY ONCOTREE CODE TABLE UPDATE
# 	oncotreeCodeDistributionDf = pd.DataFrame(columns=set(clinicalDf['CENTER']), index=set(clinicalDf['PRIMARY_CODES']))
# 	for center in oncotreeCodeDistributionDf.columns:
# 		onc_counts = clinicalDf['PRIMARY_CODES'][clinicalDf['CENTER'] == center].value_counts()
# 		oncotreeCodeDistributionDf[center] = onc_counts
# 	oncotreeCodeDistributionDf = oncotreeCodeDistributionDf.fillna(0)
# 	oncotreeCodeDistributionDf = oncotreeCodeDistributionDf.applymap(int)
# 	oncotreeCodeDistributionDf['Total'] = oncotreeCodeDistributionDf.apply(sum, axis=1)
# 	oncotreeCodeDistributionDf['Oncotree_Code'] = oncotreeCodeDistributionDf.index

# 	oncotreeDatabase = syn.tableQuery('SELECT %s FROM %s' % ("Oncotree_Code," + ",".join(clinicalDf['CENTER'].unique()) + ",Total",primaryCodeSynId))
# 	oncotreeDatabaseDf = oncotreeDatabase.asDataFrame()
# 	process.updateDatabase(syn, oncotreeDatabaseDf, oncotreeCodeDistributionDf, primaryCodeSynId, ["Oncotree_Code"], toDelete=True)
	
# 	#UPDATE DIFF TABLE
# 	sampleCountPerRound = syn.tableQuery('SELECT * FROM %s' % cumulativeSampleCountSynId)
# 	sampleCountPerRoundDf = sampleCountPerRound.asDataFrame()
# 	releases = list(set(sampleCountPerRoundDf['Release'][~sampleCountPerRoundDf['Release'].isin(['Database'])]))
# 	releases.sort()
# 	noTotalSampleCountPerRoundDf = sampleCountPerRoundDf[sampleCountPerRoundDf['Center']!='Total']
# 	difftable = noTotalSampleCountPerRoundDf[noTotalSampleCountPerRoundDf['Release'] == releases[0]]
# 	for index, rel in enumerate(releases[1:]):
# 		prior = noTotalSampleCountPerRoundDf[noTotalSampleCountPerRoundDf['Release'] == releases[index]]
# 		current = noTotalSampleCountPerRoundDf[noTotalSampleCountPerRoundDf['Release'] == rel]
# 		prior.index=prior['Center']
# 		current.index=current['Center']
# 		del prior['Center']
# 		del prior['Release']
# 		del current['Center']
# 		del current['Release']
# 		currDiff = current - prior
# 		currDiff['Center'] = currDiff.index
# 		currDiff['Release'] = rel
# 		difftable = difftable.append(currDiff)
	

# 	difftableDatabase = syn.tableQuery('SELECT * FROM %s' % sampleDiffCountSynId)
# 	difftableDatabaseDf = difftableDatabase.asDataFrame()
# 	difftableDatabaseDf = difftableDatabaseDf.fillna(0)
# 	difftable[['Clinical','Mutation','CNV','SEG','Fusions']] = difftable[['Clinical','Mutation','CNV','SEG','Fusions']].fillna(0).applymap(int)
# 	process.updateDatabase(syn, difftableDatabaseDf, difftable, sampleDiffCountSynId, ["Center","Release"])

# 	sample = syn.tableQuery('select * from syn7517674')
# 	sampleDf = sample.asDataFrame()
# 	patient = syn.tableQuery('select * from syn7517669')
# 	patientDf = patient.asDataFrame()

# 	vitalStatus = syn.tableQuery('select * from syn11559910')
# 	vitalStatusDf = vitalStatus.asDataFrame()
# 	del vitalStatusDf['CENTER']
# 	vitalStatusDf = vitalStatusDf[vitalStatusDf['PATIENT_ID'].isin(patientDf.PATIENT_ID)]
# 	patientDf = patientDf.merge(vitalStatusDf, on = "PATIENT_ID",how="outer")

# 	dataCompleteness = pd.DataFrame()
# 	centerInfos = sampleDf.CENTER.drop_duplicates().apply(lambda center: getCenterDataCompletion(center, sampleDf))
# 	for centerInfo in centerInfos:
# 		dataCompleteness = dataCompleteness.append(centerInfo)

# 	centerInfos = patientDf.CENTER.drop_duplicates().apply(lambda center: getCenterDataCompletion(center, patientDf))

# 	for centerInfo in centerInfos:
# 		dataCompleteness = dataCompleteness.append(centerInfo)
# 	dataCompletenessDatabase = syn.tableQuery('select * from %s' % dataCompletionSynId) 
# 	dataCompletenessDatabaseDf = dataCompletenessDatabase.asDataFrame()
# 	dataCompleteness.columns = dataCompletenessDatabaseDf.columns
# 	process.updateDatabase(syn, dataCompletenessDatabaseDf, dataCompleteness, dataCompletionSynId, ["FIELD","CENTER"])

# 	#Updates to query and date dashboard was updated
# 	now = datetime.datetime.now()
# 	markdown = ["_Updated %s/%s/%s_\n\n" % (now.month,now.day,now.year),
# 			   "##Count of Clinical Samples\n",
# 			   "${synapsetable?query=SELECT Center%2C Clinical%2C Release FROM " + cumulativeSampleCountSynId+ "}\n\n",
# 			   "\n\n##Oncotree Codes\n\n",
# 			   "${synapsetable?query=SELECT Oncotree%5FCode%2C " + "%2C ".join(clinicalDf['CENTER'].unique()) + "%2C Total FROM " +  primaryCodeSynId + " ORDER BY Total DESC&limit=15}\n\n"]
# 	wikiPage = syn.getWiki("syn3380222",235803)
# 	wikiPage.markdown = "".join(markdown)
# 	syn.store(wikiPage)


def main():
	parser = argparse.ArgumentParser(description='Update dashboard tables')
	parser.add_argument('--release', help = "GENIE release number (ie. 5.3-consortium)", default=None)
	parser.add_argument("--oncotree_link", type=str, help="Link to oncotree code")
	parser.add_argument("--pem_file", type=str, help="Path to PEM file (genie.pem)")
	parser.add_argument("--staging", action='store_true', help = "Using staging directory files")
	args = parser.parse_args()
	syn = process.synLogin(args)

	if args.staging:
		#Database to Synapse Id mapping Table
		database_mapping_synid = 'syn12094210'	
	else:
		database_mapping_synid = 'syn10967259'	
	
	database_mapping = syn.tableQuery('select * from %s' % database_mapping_synid)
	database_mappingdf = database_mapping.asDataFrame()
	if args.oncotree_link is None:
		oncotree_link = database_mappingdf['Id'][database_mappingdf['Database'] == 'oncotreeLink'].values[0]
		oncotree_link_ent = syn.get(oncotree_link)
		args.oncotree_link = oncotree_link_ent.externalURL
	

	update_release_numbers(syn, database_mappingdf, release = args.release)
	
	#update_tables(syn, args.oncotree_link, database_mappingdf, args.release)

if __name__ == "__main__":
	main()

