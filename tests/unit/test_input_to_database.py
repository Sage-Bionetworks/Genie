import pytest
import os
import synapseclient
import synapseutils as synu
import mock
import pandas as pd
from genie import input_to_database
# from genie import validate

syn = mock.create_autospec(synapseclient.Synapse)
sample_clinical_synid = 'syn2222'
patient_clinical_synid = 'syn11111'
vcf1synid = 'syn6666'
vcf2synid = 'syn8888'
first = (
    [('inputs', "syn12345")],
    [('vcfs', 'syn33333')],
    [('data_clinical_supp_sample_SAGE.txt', sample_clinical_synid),
     ('data_clinical_supp_patient_SAGE.txt', patient_clinical_synid)])
second = (
    [('vcfs', "syn33333")],
    [],
    [('GENIE-SAGE-000-1111.vcf', vcf1synid),
     ('GENIE-SAGE-111-2222.vcf', vcf2synid)])
center = "SAGE"
oncotreeurl = "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"


def test_samename_rename_file():
    '''
    Test file isn't renamed
    '''
    filename = synapseclient.utils.make_bogus_data_file()
    entity = synapseclient.Entity(
        path=filename, name=os.path.basename(filename))
    expectedpath = filename
    with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get:
        path = input_to_database.rename_file(syn, "syn12345")
        assert path == expectedpath
        assert os.stat(path) == os.stat(expectedpath)
        patch_syn_get.assert_called_once_with("syn12345")

    os.remove(filename)


def test_diffname_rename_file():
    '''
    Test renaming of the file
    '''
    filename = synapseclient.utils.make_bogus_data_file()
    entity = synapseclient.Entity(path=filename, name="testname")
    expectedpath = os.path.join(os.path.dirname(filename), "testname")
    with mock.patch.object(syn, "get", return_value=entity) as patch_syn_get:
        path = input_to_database.rename_file(syn, "syn12345")
        assert path == expectedpath
        assert os.stat(path) == os.stat(expectedpath)
        patch_syn_get.assert_called_once_with("syn12345")
    os.remove(filename)
    os.remove(path)


def walk_return():
    '''
    Generator returned by synu.walk
    '''
    yield first
    yield second


def walk_return_empty():
    '''
    Generator returned by synu.walk
    '''
    yield ([], [], [])


def test_main_get_center_input_files():
    '''
    Test to make sure center input files are gotten
    excluding the vcf files since process main is specified
    '''
    filename = synapseclient.utils.make_bogus_data_file()
    expected_center_file_list = [(
        [sample_clinical_synid, patient_clinical_synid],
        [filename, filename])]
    calls = [
        mock.call(syn, sample_clinical_synid),
        mock.call(syn, patient_clinical_synid)]
    with mock.patch(
            "genie.input_to_database.rename_file",
            return_value=filename) as patch_rename,\
        mock.patch.object(
            synu, "walk", return_value=walk_return()) as patch_synu_walk:
        center_file_list = input_to_database.get_center_input_files(
            syn, "syn12345", center, process="main")
        assert center_file_list == expected_center_file_list
        patch_synu_walk.assert_called_once_with(syn, 'syn12345')
        patch_rename.assert_has_calls(calls)
    os.remove(filename)


def test_vcf_get_center_input_files():
    '''
    Test to make sure center input files are gotten
    including the vcf files since process vcf is specified
    '''
    filename = synapseclient.utils.make_bogus_data_file()
    expected_center_file_list = [
        ([sample_clinical_synid, patient_clinical_synid],
         [filename, filename]),
        ([vcf1synid], [filename]), ([vcf2synid], [filename])]
    calls = [
        mock.call(syn, sample_clinical_synid),
        mock.call(syn, patient_clinical_synid),
        mock.call(syn, vcf1synid),
        mock.call(syn, vcf2synid)]

    with mock.patch(
            "genie.input_to_database.rename_file",
            return_value=filename) as patch_rename,\
        mock.patch.object(
            synu, "walk", return_value=walk_return()) as patch_synu_walk:
        center_file_list = input_to_database.get_center_input_files(
            syn, "syn12345", center, process="vcf")
        assert center_file_list == expected_center_file_list
        patch_synu_walk.assert_called_once_with(syn, 'syn12345')
        patch_rename.assert_has_calls(calls)
    os.remove(filename)


def test_empty_get_center_input_files():
    '''
    Test that center input files is empty if directory
    pass in is empty
    '''
    filename = synapseclient.utils.make_bogus_data_file()
    with mock.patch(
            "genie.input_to_database.rename_file",
            return_value=filename) as patch_rename,\
        mock.patch.object(
            synu, "walk", return_value=walk_return_empty()) as patch_synu_walk:
        center_file_list = input_to_database.get_center_input_files(
            syn, "syn12345", center, process="vcf")
        assert center_file_list == []
        patch_synu_walk.assert_called_once_with(syn, 'syn12345')
        assert not patch_rename.called
    os.remove(filename)


# @pytest.fixture(params=[
#     # tuple with (input, expectedOutput)
#     (["data_CNA_SAGE.txt"], "cna"),
#     (["data_clinical_supp_SAGE.txt"], "clinical"),
#     (["data_clinical_supp_sample_SAGE.txt",
#       "data_clinical_supp_patient_SAGE.txt"], "clinical")])
# def filename_fileformat_map(request):
#     return request.param


# def test_perfect_get_filetype(filename_fileformat_map):
#     (filepath_list, fileformat) = filename_fileformat_map
#     assert input_to_database.get_filetype(
#         syn, filepath_list, center) == fileformat


# def test_wrongfilename_get_filetype():
#     assert input_to_database.get_filetype(syn, ['wrong.txt'], center) is None


def test_unvalidatedinput_check_existing_file_status():
    validation_statusdf = pd.DataFrame(columns=['id'], dtype=str)
    error_trackerdf = pd.DataFrame(columns=['id'], dtype=str)
    entity = synapseclient.Entity(id='syn1234')
    entities = [entity]
    input_filenames = ['first.txt']
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, input_filenames)
    assert file_status['to_validate']
    assert file_status['status_list'] == []
    assert file_status['error_list'] == []


def test_valid_check_existing_file_status():
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'status': ['VALID', 'INVALID'],
        'md5': ['3333', '44444'],
        'name': ['first.txt', 'second.txt']})
    error_trackerdf = pd.DataFrame({
        'id': ['syn2345'],
        'errors': ['Invalid file format']})
    entity = synapseclient.Entity(id='syn1234', md5='3333')
    entities = [entity]
    input_filenames = ['first.txt']
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, input_filenames)
    assert not file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_invalid_check_existing_file_status():
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'status': ['VALID', 'INVALID'],
        'md5': ['3333', '44444'],
        'name': ['first.txt', 'second.txt']})
    error_trackerdf = pd.DataFrame({
        'id': ['syn2345'],
        'errors': ['Invalid file format']})
    entity = synapseclient.Entity(id='syn2345', md5='44444')
    entities = [entity]
    input_filenames = ['second.txt']
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, input_filenames)
    assert not file_status['to_validate']
    assert file_status['status_list'] == ['INVALID']
    assert file_status['error_list'] == ['Invalid file format']


def test_nostorederrors_check_existing_file_status():
    '''
    If there is no error uploaded, must re-validate file
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'status': ['VALID', 'INVALID'],
        'md5': ['3333', '44444'],
        'name': ['first.txt', 'second.txt']})
    error_trackerdf = pd.DataFrame(columns=['id'], dtype=str)
    entity = synapseclient.Entity(id='syn2345', md5='44444')
    entities = [entity]
    input_filenames = ['second.txt']
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, input_filenames)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['INVALID']
    assert file_status['error_list'] == []


def test_diffmd5validate_check_existing_file_status():
    '''
    If md5 is different from stored md5, must re-validate file
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'status': ['VALID', 'INVALID'],
        'md5': ['3333', '44444'],
        'name': ['first.txt', 'second.txt']})
    error_trackerdf = pd.DataFrame(columns=['id'], dtype=str)
    entity = synapseclient.Entity(id='syn1234', md5='44444')
    entities = [entity]
    input_filenames = ['first.txt']
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, input_filenames)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_diffnametovalidate_check_existing_file_status():
    '''
    If name is different from stored name, must re-validate file
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'status': ['VALID', 'INVALID'],
        'md5': ['3333', '44444'],
        'name': ['first.txt', 'second.txt']})
    error_trackerdf = pd.DataFrame(columns=['id'], dtype=str)
    entity = synapseclient.Entity(id='syn1234', md5='3333')
    entities = [entity]
    input_filenames = ['second.txt']
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, input_filenames)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_twoinvalid_check_existing_file_status():
    '''
    Test to make sure both invalid statues get passed back
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'status': ['INVALID', 'INVALID'],
        'md5': ['3333', '44444'],
        'name': ['first.txt', 'second.txt']})
    error_trackerdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'errors': ['Invalid file format', 'Invalid formatting issues']})
    first_entity = synapseclient.Entity(id='syn1234', md5='3333')
    second_entity = synapseclient.Entity(id='syn2345', md5='44444')
    entities = [first_entity, second_entity]
    input_filenames = ['first.txt', 'second.txt']
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, input_filenames)
    assert not file_status['to_validate']
    assert file_status['status_list'] == [
        'INVALID', 'INVALID']
    assert file_status['error_list'] == [
        'Invalid file format', 'Invalid formatting issues']


def test_error_check_existing_file_status():
    '''
    With exception of clinical files, there should never be
    more than 2 files being passed in for validation.
    '''
    with pytest.raises(
            ValueError,
            match='There should never be more than 2 files being validated.'):
        validation_statusdf = pd.DataFrame(columns=['id'], dtype=str)
        error_trackerdf = pd.DataFrame(columns=['id'], dtype=str)
        entities = ['foo', 'doo', 'boo']
        input_filenames = ['first.txt', 'doo', 'roo']
        input_to_database.check_existing_file_status(
            validation_statusdf, error_trackerdf, entities, input_filenames)


def test_create_and_archive_maf_database():
    '''
    Test the creation and archive of the maf database
    '''
    table_ent = synapseclient.Entity(
        parentId="syn123", name="foo", primaryKey=['annot'], id='syn12345')
    new_maf_ent = synapseclient.Entity(id="syn2222")
    database_synid_mappingdf = pd.DataFrame({
        'Database': ['vcf2maf', 'main'],
        'Id': ['syn12345', 'syn23455']})

    with mock.patch.object(
            syn,
            "store",
            return_value=new_maf_ent) as patch_syn_store,\
        mock.patch.object(
            syn,
            "setPermissions",
            return_value=None) as patch_syn_set_permissions,\
        mock.patch.object(
            syn,
            "get",
            return_value=table_ent) as patch_syn_get,\
        mock.patch.object(
            syn,
            "getTableColumns",
            return_value=['foo', 'ddooo']) as patch_syn_get_table_columns:

        database_mappingdf = input_to_database.create_and_archive_maf_database(
            syn, database_synid_mappingdf)

        assert database_mappingdf['Id'][
            database_mappingdf['Database'] == 'vcf2maf'].values[0] \
            == new_maf_ent.id
        assert database_mappingdf['Id'][
            database_mappingdf['Database'] == 'main'].values[0] == 'syn23455'
        patch_syn_get_table_columns.assert_called_once_with('syn12345')
        patch_syn_get.assert_called_once_with('syn12345')
        assert patch_syn_store.call_count == 3
        patch_syn_set_permissions.assert_called_once_with(
            new_maf_ent.id, 3326313, [])


def test_valid_validatefile():
    '''
    Tests the behavior of a file that gets validated that becomes
    valid
    '''
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame()
    entity = synapseclient.Entity(id='syn1234', md5='44444')
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = '333'
    entity.createdBy = '444'
    # entities = [entity]
    center = 'SAGE'
    threads = 0
    testing = False
    filetype = "clinical"
    fileinfo = {'filePaths': ['/path/to/data_clinical_supp_SAGE.txt'],
                'synId': ['syn1234']}
    # status_list = check_file_status['status_list']
    # error_list = check_file_status['error_list']
    # filetype = get_filetype(syn, filepaths, center)
    # if check_file_status['to_validate']:
    with mock.patch.object(
            syn, "get", return_value=entity) as patch_syn_get,\
        mock.patch(
            "genie.validate.determine_filetype",
            return_value=filetype) as patch_determine_filetype,\
        mock.patch(
            "genie.input_to_database.check_existing_file_status",
            return_value={
                'status_list': [],
                'error_list': [],
                'to_validate': True}) as patch_check, \
        mock.patch(
            "genie.validate.validate_single_file_workflow",
            return_value=(True, 'valid', "clinical")) as patch_validate:
        validate_results = input_to_database.validatefile(
            fileinfo, syn, validation_statusdf,
            error_trackerdf, center, threads, testing, oncotreeurl)
        expected_validate_results = ([[
            fileinfo['synId'][0],
            fileinfo['filePaths'][0],
            '44444',
            'VALIDATED',
            'data_clinical_supp_SAGE.txt',
            1553428800000,
            'clinical']], None)
        assert expected_validate_results == validate_results
        patch_validate.assert_called_once()
        patch_syn_get.assert_called_once()
        patch_check.assert_called_once()
        patch_determine_filetype.assert_called_once_with(
            syn, fileinfo['filePaths'], center, raise_error=False)


def test_invalid_validatefile():
    '''
    Tests the behavior of a file that gets validated that becomes
    invalid
    '''
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame(columns=['id'], dtype=str)
    entity = synapseclient.Entity(id='syn2345', md5='44444')
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    entity.modifiedBy = '333'
    entity.createdBy = '333'
    # entities = [entity]
    center = 'SAGE'
    threads = 0
    testing = False
    check_file_status_dict = {
        'status_list': [],
        'error_list': [],
        'to_validate': True}
    filetype = "clinical"
    fileinfo = {'filePaths': ['/path/to/data_clinical_supp_SAGE.txt'],
                'synId': ['syn1234']}
    with mock.patch.object(
            syn, "get", return_value=entity) as patch_syn_get,\
        mock.patch(
            "genie.validate.determine_filetype",
            return_value=filetype) as patch_determine_filetype,\
        mock.patch(
            "genie.input_to_database.check_existing_file_status",
            return_value=check_file_status_dict) as patch_check,\
        mock.patch.object(
            syn, "getUserProfile",
            return_value={'userName': 'trial'}) as patch_syn_getuserprofile,\
        mock.patch.object(
            syn, "sendMessage") as patch_syn_sendmessage,\
        mock.patch(
            "genie.validate.validate_single_file_workflow",
            return_value=(False, 'invalid', "clinical")) as patch_validate:
        foo = input_to_database.validatefile(
            fileinfo, syn, validation_statusdf,
            error_trackerdf, center, threads, testing, oncotreeurl)
        patch_validate.assert_called_once()
        patch_check.assert_called_once()
        patch_determine_filetype.assert_called_once_with(
            syn, fileinfo['filePaths'], center, raise_error=False)
        error_message = (
            "Dear trial,\n\n"
            "Your files (data_clinical_supp_SAGE.txt) are invalid! "
            "Here are the reasons why:\n\ninvalid")
        patch_syn_sendmessage.assert_called_once_with(
            ['333'], "GENIE Validation Error", error_message)
        patch_syn_get.call_count == 3
        patch_syn_getuserprofile.call_count == 2
        expected_validate_results = (
            [['syn2345',
              '/path/to/data_clinical_supp_SAGE.txt',
              '44444',
              'INVALID',
              'data_clinical_supp_SAGE.txt',
              1553428800000,
              'clinical']],
            [['syn1234', 'invalid', 'data_clinical_supp_SAGE.txt']])
        assert foo == expected_validate_results


def test_already_validated_validatefile():
    '''
    Test already validated files
    '''
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame()
    entity = synapseclient.Entity(id='syn1234', md5='44444')
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = '333'
    entity.createdBy = '444'
    # entities = [entity]
    center = 'SAGE'
    threads = 0
    testing = False

    fileinfo = {'filePaths': ['/path/to/data_clinical_supp_SAGE.txt'],
                'synId': ['syn1234']}
    filetype = "markdown"
    # status_list = check_file_status['status_list']
    # error_list = check_file_status['error_list']
    # filetype = get_filetype(syn, filepaths, center)
    # if check_file_status['to_validate']:
    check_file_status_dict = {
        'status_list': ["INVALID"],
        'error_list': ["invalid file"],
        'to_validate': False}
    with mock.patch.object(
            syn, "get", return_value=entity) as patch_syn_get,\
        mock.patch(
            "genie.validate.determine_filetype",
            return_value=filetype) as patch_determine_filetype,\
        mock.patch(
            "genie.input_to_database.check_existing_file_status",
            return_value=check_file_status_dict) as patch_check:
        validate_results = input_to_database.validatefile(
            fileinfo, syn, validation_statusdf,
            error_trackerdf, center, threads, testing, oncotreeurl)
        expected_validate_results = (
            [[fileinfo['synId'][0],
              fileinfo['filePaths'][0],
              entity.md5,
              check_file_status_dict['status_list'][0],
              'data_clinical_supp_SAGE.txt',
              1553428800000,
              filetype]],
            [[fileinfo['synId'][0],
              check_file_status_dict['error_list'][0],
              'data_clinical_supp_SAGE.txt']])

        assert expected_validate_results == validate_results
        patch_syn_get.assert_called_once()
        patch_check.assert_called_once()
        patch_determine_filetype.assert_called_once_with(
            syn, fileinfo['filePaths'], center, raise_error=False)


# def test_filetypenone__check_valid():
#     input_to_database._check_valid(
#         syn, filepaths, center, filetype, filenames,
#          oncotree_link, threads, testing)
