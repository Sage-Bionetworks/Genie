import pytest
import os
import synapseclient
import synapseutils as synu
import mock
import pandas as pd
from genie import input_to_database
from genie import validate

syn = mock.create_autospec(synapseclient.Synapse) 
sample_clinical_synid = 'syn2222'
patient_clincal_synid = 'syn11111'
vcf1synid = 'syn6666'
vcf2synid = 'syn8888'
first = ([('inputs',"syn12345")], [('vcfs','syn33333')], [('data_clinical_supp_sample_SAGE.txt',sample_clinical_synid), ('data_clinical_supp_patient_SAGE.txt',patient_clincal_synid)])
second = ([('vcfs',"syn33333")], [], [('GENIE-SAGE-000-1111.vcf',vcf1synid), ('GENIE-SAGE-111-2222.vcf',vcf2synid)])
center = "SAGE"
oncotreeurl = "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"

def test_samename_rename_file():
	filename = synapseclient.utils.make_bogus_data_file()
	entity = synapseclient.Entity(path = filename, name=os.path.basename(filename))
	expectedpath = filename
	with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get:
		path = input_to_database.rename_file(syn, "syn12345")
		assert path == expectedpath
		assert os.stat(path) == os.stat(expectedpath)
		patch_syn_get.assert_called_once_with("syn12345")

	os.remove(filename)

def test_diffname_rename_file():
	filename = synapseclient.utils.make_bogus_data_file()
	entity = synapseclient.Entity(path = filename, name="testname")
	expectedpath = os.path.join(os.path.dirname(filename), "testname")
	with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get:
		path = input_to_database.rename_file(syn, "syn12345")
		assert path == expectedpath
		assert os.stat(path) == os.stat(expectedpath)
		patch_syn_get.assert_called_once_with("syn12345")
	os.remove(filename)
	os.remove(path)


def walk_return():
	yield first
	yield second

def test_main_get_center_input_files():
	filename = synapseclient.utils.make_bogus_data_file()
	entity = synapseclient.Entity(path = filename, name=os.path.basename(filename))
	expected_center_file_list = [([sample_clinical_synid,patient_clincal_synid],[filename,filename])]
	with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get,\
		 mock.patch.object(synu, "walk", return_value=walk_return()) as patch_synu_walk:
		center_file_list = input_to_database.get_center_input_files(syn, "syn12345", center, process="main")
		assert center_file_list == expected_center_file_list
		patch_synu_walk.assert_called_once_with(syn, 'syn12345')
	os.remove(filename)

def test_vcf_get_center_input_files():
	filename = synapseclient.utils.make_bogus_data_file()
	entity = synapseclient.Entity(path = filename, name=os.path.basename(filename))
	expected_center_file_list = [([sample_clinical_synid,patient_clincal_synid],[filename,filename]),([vcf1synid],[filename]),([vcf2synid],[filename])]
	with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get,\
		 mock.patch.object(synu, "walk", return_value=walk_return()) as patch_synu_walk:
		center_file_list = input_to_database.get_center_input_files(syn, "syn12345", center, process="vcf")
		assert center_file_list == expected_center_file_list
		patch_synu_walk.assert_called_once_with(syn, 'syn12345')
	os.remove(filename)

def walk_return_empty():
	yield ([], [], [])

def test_empty_get_center_input_files():
	filename = synapseclient.utils.make_bogus_data_file()
	entity = synapseclient.Entity(path = filename, name=os.path.basename(filename))
	expected_center_file_list = [([sample_clinical_synid,patient_clincal_synid],[filename,filename]),([vcf1synid],[filename]),([vcf2synid],[filename])]
	with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get,\
		 mock.patch.object(synu, "walk", return_value=walk_return_empty()) as patch_synu_walk:
		center_file_list = input_to_database.get_center_input_files(syn, "syn12345", center, process="vcf")
		assert center_file_list == []
	os.remove(filename)


@pytest.fixture( params=[
	# tuple with (input, expectedOutput)
	(["data_CNA_SAGE.txt"], "cna"),
	(["data_clinical_supp_SAGE.txt"], "clinical"),
	(["data_clinical_supp_sample_SAGE.txt","data_clinical_supp_patient_SAGE.txt"], "clinical")
	])
def filename_fileformat_map(request):
	return request.param

def test_perfect_get_filetype(filename_fileformat_map):
	(filepath_list, fileformat) = filename_fileformat_map
	assert input_to_database.get_filetype(syn, filepath_list, center) == fileformat

def test_wrongfilename_get_filetype():
	assert input_to_database.get_filetype(syn, ['wrong.txt'], center) is None

def test_unvalidatedinput_check_existing_file_status():
	validation_statusdf = pd.DataFrame(columns=['id'],dtype=str)
	error_trackerdf = pd.DataFrame(columns=['id'],dtype=str)
	entity = synapseclient.Entity(id='syn1234')
	entities = [entity]
	input_filenames = ['first.txt']
	file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	assert file_status['to_validate']
	assert file_status['status_list'] == []
	assert file_status['error_list'] == []

def test_valid_check_existing_file_status():
	validation_statusdf = pd.DataFrame({'id':['syn1234','syn2345'],'status':['VALID','INVALID'],'md5':['3333','44444'],'name':['first.txt','second.txt']})
	error_trackerdf = pd.DataFrame({'id':['syn2345'],'errors':['Invalid file format']})
	entity = synapseclient.Entity(id='syn1234',md5='3333')
	entities = [entity]
	input_filenames = ['first.txt']
	file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	assert not file_status['to_validate']
	assert file_status['status_list'] == ['VALID']
	assert file_status['error_list'] == []	

def test_invalid_check_existing_file_status():
	validation_statusdf = pd.DataFrame({'id':['syn1234','syn2345'],'status':['VALID','INVALID'],'md5':['3333','44444'],'name':['first.txt','second.txt']})
	error_trackerdf = pd.DataFrame({'id':['syn2345'],'errors':['Invalid file format']})
	entity = synapseclient.Entity(id='syn2345',md5='44444')
	entities = [entity]
	input_filenames = ['second.txt']
	file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	assert not file_status['to_validate']
	assert file_status['status_list'] == ['INVALID']
	assert file_status['error_list'] == ['Invalid file format']	

def test_nostorederrors_check_existing_file_status():
	validation_statusdf = pd.DataFrame({'id':['syn1234','syn2345'],'status':['VALID','INVALID'],'md5':['3333','44444'],'name':['first.txt','second.txt']})
	error_trackerdf = pd.DataFrame(columns=['id'],dtype=str)
	entity = synapseclient.Entity(id='syn2345',md5='44444')
	entities = [entity]
	input_filenames = ['second.txt']
	file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	assert file_status['to_validate']
	assert file_status['status_list'] == ['INVALID']
	assert file_status['error_list'] == []	

	
def test_diffmd5validate_check_existing_file_status():
	validation_statusdf = pd.DataFrame({'id':['syn1234','syn2345'],'status':['VALID','INVALID'],'md5':['3333','44444'],'name':['first.txt','second.txt']})
	error_trackerdf = pd.DataFrame(columns=['id'],dtype=str)
	entity = synapseclient.Entity(id='syn1234',md5='44444')
	entities = [entity]
	input_filenames = ['first.txt']
	file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	assert file_status['to_validate']
	assert file_status['status_list'] == ['VALID']
	assert file_status['error_list'] == []	

def test_diffnametovalidate_check_existing_file_status():
	validation_statusdf = pd.DataFrame({'id':['syn1234','syn2345'],'status':['VALID','INVALID'],'md5':['3333','44444'],'name':['first.txt','second.txt']})
	error_trackerdf = pd.DataFrame(columns=['id'],dtype=str)
	entity = synapseclient.Entity(id='syn1234',md5='3333')
	entities = [entity]
	input_filenames = ['second.txt']
	file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	assert file_status['to_validate']
	assert file_status['status_list'] == ['VALID']
	assert file_status['error_list'] == []	

def test_twoinvalid_check_existing_file_status():
	validation_statusdf = pd.DataFrame({'id':['syn1234','syn2345'],'status':['INVALID','INVALID'],'md5':['3333','44444'],'name':['first.txt','second.txt']})
	error_trackerdf = pd.DataFrame({'id':['syn1234','syn2345'],'errors':['Invalid file format','Invalid formatting issues']})
	first_entity = synapseclient.Entity(id='syn1234',md5='3333')
	second_entity = synapseclient.Entity(id='syn2345',md5='44444')
	entities = [first_entity, second_entity]
	input_filenames = ['first.txt','second.txt']
	file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	assert not file_status['to_validate']
	assert file_status['status_list'] == ['INVALID','INVALID']
	assert file_status['error_list'] == ['Invalid file format','Invalid formatting issues']	

def test_error_check_existing_file_status():
	with pytest.raises(ValueError, match='There should never be more than 2 files being validated.'):
		validation_statusdf = pd.DataFrame(columns=['id'],dtype=str)
		error_trackerdf = pd.DataFrame(columns=['id'],dtype=str)
		entities = ['foo','doo','boo']
		input_filenames = ['first.txt','doo','roo']
		file_status = input_to_database.check_existing_file_status(validation_statusdf, error_trackerdf, entities, input_filenames)
	
def test_create_and_archive_maf_database():
	table_ent = synapseclient.Entity(parentId="syn123",name="foo",primaryKey=['annot'],id='syn12345')
	new_maf_ent = synapseclient.Entity(id="syn2222")
	database_synid_mappingdf = pd.DataFrame({'Database':['vcf2maf','main'], 'Id':['syn12345','syn23455']})
	with mock.patch.object(syn, "store", return_value=new_maf_ent) as patch_syn_store,\
		 mock.patch.object(syn, "setPermissions", return_value=None) as patch_syn_set_permissions,\
		 mock.patch.object(syn, "get", return_value=table_ent) as patch_syn_get,\
		 mock.patch.object(syn, "getTableColumns", return_value=['foo','ddooo']) as patch_syn_get_table_columns:

		database_mappingdf = input_to_database.create_and_archive_maf_database(syn, database_synid_mappingdf)

		assert database_mappingdf['Id'][database_mappingdf['Database'] == 'vcf2maf'].values[0] == new_maf_ent.id
		assert database_mappingdf['Id'][database_mappingdf['Database'] == 'main'].values[0] == 'syn23455'


def test_row_validatefile():
	validation_statusdf = pd.DataFrame({'id':['syn1234','syn2345'],'status':['VALID','INVALID'],'md5':['3333','44444'],'name':['first.txt','second.txt']})
	error_trackerdf = pd.DataFrame(columns=['id'],dtype=str)
	entity = synapseclient.Entity(id='syn2345',md5='44444')
	entity['modifiedOn']='2019-03-24T12:00:00.Z'
	entity.modifiedBy='333'
	entity.createdBy='444'
	#entities = [entity]
	center = 'SAGE'
	threads=0
	testing=False

	fileinfo = {'filePaths':['/path/to/file.csv'],
				'synId':['syn1234']}
	with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get,\
		 mock.patch.object(validate, "validate", return_value=('valid',True)) as patch_validate:
		input_to_database.validatefile(fileinfo, syn, validation_statusdf, error_trackerdf, center, threads, testing, oncotreeurl)


