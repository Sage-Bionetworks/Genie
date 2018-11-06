import genie
import os
import synapseclient
import pandas as pd
import argparse
import datetime

def get_center_data_completion(center, df):
	'''
	Get center data completion.  Calulates the percentile of 
	how complete a clinical data element is:
	Number of not blank/Unknown/NA divded by total number of patients or samples

	params:
		center: GENIE center
		df: sample or patient dataframe
	'''
	centerdf = df[df['CENTER'] == center]
	total = len(centerdf)
	center_data = pd.DataFrame()
	for col in centerdf:
		if not col.endswith("_NUMERICAL") and col not in ['CENTER', 'PATIENT_ID','SAMPLE_ID','SAMPLE_TYPE_DETAILED']:
			not_missing = [pd.isnull(value) for value in centerdf[col]]
			completeness = float(sum(not_missing)) / int(total)
			returned = pd.DataFrame([[col, center, total, completeness]])
			center_data = center_data.append(returned)
	return(center_data)

def update_samples_in_release_table(syn, file_mapping, release, samples_in_release_synid):
	'''
	Convenience function that updates the sample in release table
	This tracks the samples of each release.  1 means it exists, and 0 means it doesn't

	params: 
		syn: synapse object
		file_mapping: file mapping generated from file mapping function
		release:  GENIE release number (ie. 5.3-consortium)
		samples_in_release_synid: Synapse Id of 'samples in release' Table
	'''
	clinical_ent = syn.get(file_mapping['clinical'],followLink=True)
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
	genie.process_functions.updateDatabase(syn, samples_per_releasedf, samples_in_releasedf, samples_in_release_synid, ["SAMPLE_ID"])


def update_cumulative_sample_table(syn, file_mapping, release, cumulative_sample_count_synid):
	'''
	Consortium release sample count table update function
	This gets the cumulative sample count of each file type in each release
	
	params:
		syn: synapse object
		file_mapping: file mapping generated from file mapping function
		release:  GENIE release number (ie. 5.3-consortium)
		cumulative_sample_count_synid: Synapse Id of 'Cumulative sample count' Table
	'''

	sample_count_per_round = syn.tableQuery('SELECT * FROM %s' % cumulative_sample_count_synid)
	sample_count_per_rounddf = sample_count_per_round.asDataFrame()

	clinical_ent = syn.get(file_mapping['clinical'],followLink=True)
	clinicaldf = pd.read_csv(clinical_ent.path,sep="\t",comment="#")
	clinicaldf.columns = [i.upper() for i in clinicaldf.columns]
	if clinicaldf.get("CENTER") is None:
		clinicaldf['CENTER'] = [sample.split("-")[1] for sample in clinicaldf.SAMPLE_ID]
	clinical_counts = clinicaldf['CENTER'].value_counts()
	clinical_counts['Total'] = sum(clinical_counts)
	clinical_counts.name = "Clinical"

	fusion_ent = syn.get(file_mapping['fusion'],followLink=True)
	fusiondf = pd.read_csv(fusion_ent.path, sep="\t",comment="#")
	fusiondf.columns = [i.upper() for i in fusiondf.columns]

	fusion_counts = fusiondf['CENTER'][~fusiondf['TUMOR_SAMPLE_BARCODE'].duplicated()].value_counts()
	fusion_counts['Total'] = sum(fusion_counts)

	cna_ent = syn.get(file_mapping['cna'],followLink=True)
	cnadf = pd.read_csv(cna_ent.path, sep="\t", comment="#")
	cna_counts = pd.Series([i.split("-")[1] for i in cnadf.columns[1:]]).value_counts()
	cna_counts['Total'] = sum(cna_counts)

	seg_ent= syn.get(file_mapping['seg'],followLink=True)
	segdf = pd.read_csv(seg_ent.path, sep="\t", comment="#")
	segdf.columns = [i.upper() for i in segdf.columns]

	segdf['CENTER'] = [i.split("-")[1] for i in segdf['ID']]
	seg_counts = segdf['CENTER'][~segdf['ID'].duplicated()].value_counts()
	seg_counts['Total'] = sum(seg_counts)

	total_counts = pd.DataFrame(clinical_counts)
	total_counts['Fusions'] = fusion_counts
	total_counts['CNV'] = cna_counts
	total_counts['Mutation'] = clinical_counts
	total_counts['SEG'] = seg_counts
	total_counts = total_counts.fillna(0)
	total_counts = total_counts.applymap(int)
	total_counts['Center'] =  total_counts.index

	total_counts['Release'] = release
	genie.process_functions.updateDatabase(syn, sample_count_per_rounddf, total_counts, cumulative_sample_count_synid, ["Center", "Release"])


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
	'''
	Function that updates all release dashboard numbers or
	specific release number

	params:
		syn: synapse object
		database_mappingdf: mapping between synapse ids and database
		release: GENIE release (ie. 5.3-consortium).  Defaults to None
	'''
	#Update release table with current release or all releases
	samples_in_release_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'samplesInRelease'].values[0]
	cumulative_sample_count_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'cumulativeSampleCount'].values[0]

	release_folder_fileview_synid = "syn17019650"
	release_folder = syn.tableQuery("select id,name from %s" % release_folder_fileview_synid + " where name not like 'Release%' and name <> 'case_lists'")
	release_folderdf = release_folder.asDataFrame()

	for release_synid, release_name in zip(release_folderdf.id, release_folderdf.name):
		file_mapping = get_file_mapping(syn, release_synid)
		#If release is specified, only process on that, otherwise process for all
		if release is None or release == release_name:
			update_samples_in_release_table(syn, file_mapping, release_name, samples_in_release_synid)
			update_cumulative_sample_table(syn, file_mapping, release_name, cumulative_sample_count_synid)
		else:
			pass

def update_database_numbers(syn, database_mappingdf):
	'''
	Updates database cumulative numbers
	
	params:
		syn: synapse object
		database_mappingdf: mapping between synapse ids and database
	'''
	cumulative_sample_count_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'cumulativeSampleCount'].values[0]
	# ## Database 
	database_count = syn.tableQuery("SELECT * FROM %s where Release = 'Database'" % cumulative_sample_count_synid)
	database_countdf = database_count.asDataFrame()
	clinical = syn.tableQuery('select * from syn7517674')
	clinicaldf = clinical.asDataFrame()
	clinincal_counts = clinicaldf['CENTER'].value_counts()
	clinincal_counts['Total'] = sum(clinincal_counts)
	clinincal_counts.name = "Clinical"

	fusion = syn.tableQuery('select * from syn7893268')
	fusiondf = fusion.asDataFrame()
	fusion_counts = fusiondf['CENTER'][~fusiondf['TUMOR_SAMPLE_BARCODE'].duplicated()].value_counts()
	fusion_counts['Total'] = sum(fusion_counts)

	center_flat_files = syn.getChildren("syn12278118")
	cna_file_paths = [syn.get(file['id']).path for file in center_flat_files if file['name'].startswith("data_CNA")]
	cna_numbers = {}
	for cna_file in cna_file_paths:
		center = os.path.basename(cna_file).replace(".txt","").split("_")[2]
		with open(cna_file,'r') as cna:
			header = cna.readline()
			samples = header.split("\t")
			# Minus one because of Hugo_Symbol
			cna_numbers[center] = len(samples) - 1
	cna_counts = pd.Series(cna_numbers)
	cna_counts['Total'] = sum(cna_counts)

	seg = syn.tableQuery('select * from syn7893341')
	segdf = seg.asDataFrame()
	seg_counts = segdf['CENTER'][~segdf['ID'].duplicated()].value_counts()
	seg_counts['Total'] = sum(seg_counts)

	db_counts = pd.DataFrame(clinincal_counts)
	db_counts['Fusions'] = fusion_counts
	db_counts['CNV'] = cna_counts
	db_counts['Mutation'] = clinincal_counts
	db_counts['SEG'] = seg_counts
	db_counts = db_counts.fillna(0)
	db_counts = db_counts.applymap(int)
	db_counts['Center'] =  db_counts.index
	db_counts['Release'] = "Database"
	genie.process_functions.updateDatabase(syn, database_countdf, db_counts, cumulative_sample_count_synid, ["Center", "Release"])


def update_oncotree_code_tables(syn, database_mappingdf):
	'''
	Function that updates database statistics of oncotree codes and primary onocotree codes

	params:
		syn: synapse object
		database_mappingdf: mapping between synapse ids and database
	'''
	oncotree_distribution_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'oncotree'].values[0]
	clinical = syn.tableQuery('select * from syn7517674')
	clinicaldf = clinical.asDataFrame()

	#### DISTRIBUTION OF ONCOTREE CODE TABLE UPDATE
	oncotree_code_distributiondf = pd.DataFrame(columns=set(clinicaldf['CENTER']), index=set(clinicaldf['ONCOTREE_CODE']))
	for center in oncotree_code_distributiondf.columns:
		onc_counts = clinicaldf['ONCOTREE_CODE'][clinicaldf['CENTER'] == center].value_counts()
		oncotree_code_distributiondf[center] = onc_counts
	oncotree_code_distributiondf = oncotree_code_distributiondf.fillna(0)
	oncotree_code_distributiondf = oncotree_code_distributiondf.applymap(int)
	oncotree_code_distributiondf['Total'] = oncotree_code_distributiondf.apply(sum, axis=1)
	oncotree_code_distributiondf['Oncotree_Code'] = oncotree_code_distributiondf.index

	oncotree_distribution_db = syn.tableQuery('SELECT %s FROM %s' % ("Oncotree_Code," + ",".join(clinicaldf['CENTER'].unique()) + ",Total", oncotree_distribution_synid))
	oncotree_distribution_dbdf = oncotree_distribution_db.asDataFrame()
	genie.process_functions.updateDatabase(syn, oncotree_distribution_dbdf, oncotree_code_distributiondf, oncotree_distribution_synid, ["Oncotree_Code"],toDelete=True)
	
	#### DISTRIBUTION OF PRIMARY CODE TABLE UPDATE
	oncotree_link_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'oncotreeLink'].values[0]
	primary_code_synId = database_mappingdf['Id'][database_mappingdf['Database'] == 'primaryCode'].values[0]
	#Can also use most up to date oncotree code, because these tables are updated
	#from the database
	oncotree_link_ent = syn.get(oncotree_link_synid)
	oncotree_link = oncotree_link_ent.externalURL
	oncotree_mapping = genie.process_functions.get_oncotree_code_mappings(oncotree_link)
	clinicaldf['PRIMARY_CODES'] = [oncotree_mapping[i.upper()]['ONCOTREE_PRIMARY_NODE'] for i in clinicaldf.ONCOTREE_CODE]

	# ### DISTRIBUTION OF PRIMARY ONCOTREE CODE TABLE UPDATE
	primary_code_distributiondf = pd.DataFrame(columns=set(clinicaldf['CENTER']), index=set(clinicaldf['PRIMARY_CODES']))
	for center in primary_code_distributiondf.columns:
		onc_counts = clinicaldf['PRIMARY_CODES'][clinicaldf['CENTER'] == center].value_counts()
		primary_code_distributiondf[center] = onc_counts
	primary_code_distributiondf = primary_code_distributiondf.fillna(0)
	primary_code_distributiondf = primary_code_distributiondf.applymap(int)
	primary_code_distributiondf['Total'] = primary_code_distributiondf.apply(sum, axis=1)
	primary_code_distributiondf['Oncotree_Code'] = primary_code_distributiondf.index

	primary_code_dist_db = syn.tableQuery('SELECT %s FROM %s' % ("Oncotree_Code," + ",".join(clinicaldf['CENTER'].unique()) + ",Total",primary_code_synId))
	primary_code_dist_dbdf = primary_code_dist_db.asDataFrame()
	genie.process_functions.updateDatabase(syn, primary_code_dist_dbdf, primary_code_distributiondf, primary_code_synId, ["Oncotree_Code"], toDelete=True)
	

def update_sample_difference_table(syn, database_mappingdf):
	'''
	Function that updates sample difference table between
	consortium releases

	params:
		syn: synapse object
		database_mappingdf: mapping between synapse ids and database
	'''
	cumulative_sample_count_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'cumulativeSampleCount'].values[0]
	sample_diff_count_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'sampleDiffCount'].values[0]

	#UPDATE DIFF TABLE
	sample_count_per_round = syn.tableQuery("SELECT * FROM %s where Center <> 'Total' and Release <> 'Database'" % cumulative_sample_count_synid)
	sample_count_per_rounddf = sample_count_per_round.asDataFrame()
	releases = list(sample_count_per_rounddf['Release'].unique())
	#sort the releases and remove public releases
	releases.sort()
	consortium_releases = [release for release in releases if "public" not in release and ".0." not in release]
	diff_between_releasesdf = sample_count_per_rounddf[sample_count_per_rounddf['Release'] == consortium_releases[0]]
	for index, release_name in enumerate(consortium_releases[1:]):
		prior_release = sample_count_per_rounddf[sample_count_per_rounddf['Release'] == consortium_releases[index]]
		current_release = sample_count_per_rounddf[sample_count_per_rounddf['Release'] == release_name]
		prior_release.index=prior_release['Center']
		current_release.index=current_release['Center']

		del prior_release['Center']
		del prior_release['Release']
		del current_release['Center']
		del current_release['Release']
		#Append new rows of centers that are new and just added to the releases
		new_centers = current_release.index[~current_release.index.isin(prior_release.index)]
		if not new_centers.empty:
			prior_release = prior_release.append(pd.DataFrame(index=new_centers))
			prior_release = prior_release.fillna(0)
		difference = current_release - prior_release
		difference['Center'] = difference.index
		difference['Release'] = release_name
		diff_between_releasesdf = diff_between_releasesdf.append(difference)
	
	difftable_db = syn.tableQuery('SELECT * FROM %s' % sample_diff_count_synid)
	difftable_dbdf = difftable_db.asDataFrame()
	difftable_dbdf = difftable_dbdf.fillna(0)
	diff_between_releasesdf[['Clinical','Mutation','CNV','SEG','Fusions']] = diff_between_releasesdf[['Clinical','Mutation','CNV','SEG','Fusions']].fillna(0).applymap(int)
	genie.process_functions.updateDatabase(syn, difftable_dbdf, diff_between_releasesdf, sample_diff_count_synid, ["Center","Release"], toDelete=True)


def update_data_completeness_table(syn, database_mappingdf):
	'''
	Function that updates the data completeness of the database
	
	params:
		syn: synapse object
		database_mappingdf: mapping between synapse ids and database
	'''
	data_completion_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'dataCompletion'].values[0]

	sample = syn.tableQuery('select * from syn7517674')
	sampledf = sample.asDataFrame()
	patient = syn.tableQuery('select * from syn7517669')
	patientdf = patient.asDataFrame()

	vital_status = syn.tableQuery('select * from syn11559910')
	vital_statusdf = vital_status.asDataFrame()
	del vital_statusdf['CENTER']
	vital_statusdf = vital_statusdf[vital_statusdf['PATIENT_ID'].isin(patientdf.PATIENT_ID)]
	patientdf = patientdf.merge(vital_statusdf, on = "PATIENT_ID",how="outer")

	data_completenessdf = pd.DataFrame()
	center_infos = sampledf.CENTER.drop_duplicates().apply(lambda center: get_center_data_completion(center, sampledf))
	for center_info in center_infos:
		data_completenessdf = data_completenessdf.append(center_info)

	center_infos = patientdf.CENTER.drop_duplicates().apply(lambda center: get_center_data_completion(center, patientdf))
	for center_info in center_infos:
		data_completenessdf = data_completenessdf.append(center_info)

	data_completeness_db = syn.tableQuery('select * from %s' % data_completion_synid) 
	data_completeness_dbdf = data_completeness_db.asDataFrame()
	data_completenessdf.columns = data_completeness_dbdf.columns
	genie.process_functions.updateDatabase(syn, data_completeness_dbdf, data_completenessdf, data_completion_synid, ["FIELD","CENTER"], toDelete=True)


def update_wiki(syn, database_mappingdf):
	'''
	Updates the GENIE project dashboard wiki timestamp

	params:
		syn: synapse object
		database_mappingdf: mapping between synapse ids and database

	'''
	#Updates to query and date dashboard was updated
	cumulative_sample_count_synid = database_mappingdf['Id'][database_mappingdf['Database'] == 'cumulativeSampleCount'].values[0]
	primary_code_synId = database_mappingdf['Id'][database_mappingdf['Database'] == 'primaryCode'].values[0]

	centers = syn.tableQuery('select distinct(CENTER) as CENTER from syn7517674')
	centersdf = centers.asDataFrame()
	now = datetime.datetime.now()
	markdown = ["_Updated %s/%s/%s_\n\n" % (now.month,now.day,now.year),
			   "##Count of Clinical Samples\n",
			   "${synapsetable?query=SELECT Center%2C Clinical%2C Release FROM " + cumulative_sample_count_synid+ "}\n\n",
			   "\n\n##Primary Oncotree Codes\n\n",
			   "${synapsetable?query=SELECT Oncotree%5FCode%2C " + "%2C ".join(centersdf['CENTER'].unique()) + "%2C Total FROM " +  primary_code_synId + " ORDER BY Total DESC&limit=15}\n\n"]
	wikiPage = syn.getWiki("syn3380222",235803)
	wikiPage.markdown = "".join(markdown)
	syn.store(wikiPage)

def run_dashboard(syn, database_mappingdf, release):
	'''
	Function that runs the dashboard scripts

	params:
		syn: synapse object
		database_mappingdf: mapping between synapse ids and database
		release: GENIE release (ie. 5.3-consortium)

	'''
	update_release_numbers(syn, database_mappingdf, release = release)
	update_database_numbers(syn, database_mappingdf)
	update_oncotree_code_tables(syn, database_mappingdf)
	update_sample_difference_table(syn, database_mappingdf)
	update_data_completeness_table(syn, database_mappingdf)
	update_wiki(syn,database_mappingdf)

def main():
	parser = argparse.ArgumentParser(description='Update dashboard tables')
	parser.add_argument('--release', help = "GENIE release number (ie. 5.3-consortium)", default=None)
	parser.add_argument("--pem_file", type=str, help="Path to PEM file (genie.pem)")
	parser.add_argument("--staging", action='store_true', help = "Using staging directory files")
	args = parser.parse_args()
	syn = genie.process_functions.synLogin(args)
	if args.staging:
		#Database to Synapse Id mapping Table
		database_mapping_synid = 'syn12094210'	
	else:
		database_mapping_synid = 'syn10967259'	
	
	database_mapping = syn.tableQuery('select * from %s' % database_mapping_synid)
	database_mappingdf = database_mapping.asDataFrame()
	run_dashboard(syn, database_mappingdf, args.release)


if __name__ == "__main__":
	main()

