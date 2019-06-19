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
    '''Test that the file path is not renamed.
    '''
    filename = synapseclient.utils.make_bogus_data_file()
    entity = synapseclient.File(path=filename,
                                id='syn012345',
                                parentId='syn45678',
                                name=os.path.basename(filename))
    expectedpath = filename
    new_entity = input_to_database.rename_file(syn, entity)
    assert new_entity.annotations.expectedPath == expectedpath
    os.remove(filename)


def test_diffname_rename_file():
    '''Test that the file path is renamed.
    '''
    filename = synapseclient.utils.make_bogus_data_file()
    entity = synapseclient.File(path=filename,
                                id='syn012345',
                                parentId='syn45678',
                                name='testname')

    expectedpath = os.path.join(os.path.dirname(filename), "testname")
    new_entity = input_to_database.rename_file(syn, entity)
    assert new_entity.annotations.expectedPath == expectedpath
    os.remove(filename)

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
            "genie.validate.validate_single_file",
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
            syn, fileinfo['filePaths'], center)


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
            "genie.validate.validate_single_file",
            return_value=(False, 'invalid', "clinical")) as patch_validate:
        foo = input_to_database.validatefile(
            fileinfo, syn, validation_statusdf,
            error_trackerdf, center, threads, testing, oncotreeurl)
        patch_validate.assert_called_once()
        patch_check.assert_called_once()
        patch_determine_filetype.assert_called_once_with(
            syn, fileinfo['filePaths'], center)
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
            syn, fileinfo['filePaths'], center)


# def test_filetypenone__check_valid():
#     input_to_database._check_valid(
#         syn, filepaths, center, filetype, filenames,
#          oncotree_link, threads, testing)


def test_dups_get_duplicated_files():
    '''
    Test get all duplicates
    cbs/seg
    clinical
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345', 'syn5555', 'syn1224', 'syn34444'],
        'name': ['first.cbs', 'second.seg', 'data_clinical_supp_1',
                 'data_clinical_supp_2', 'data_clinical_supp_3']})
    expected_dup = validation_statusdf.copy()
    expected_dup['errors'] = ''
    dupsdf = input_to_database.get_duplicated_files(
        syn, validation_statusdf, "")
    assert dupsdf.equals(expected_dup)


def test_nodups_get_duplicated_files():
    '''
    Test no duplicated
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345', 'syn5555', 'syn1224', 'syn34444'],
        'name': ['cbs.txt', 'second.seg', 'no_clinical.txt',
                 'data_clinical_supp_2', 'data_clinical_supp_3']})
    dupsdf = input_to_database.get_duplicated_files(
        syn, validation_statusdf, "")
    assert dupsdf.empty


def test_dups_email_duplication_error():
    '''
    Test duplicated email sent
    '''
    duplicated_filesdf = pd.DataFrame({
        'id': ['syn1234'],
        'name': ['first.cbs']})
    entity = synapseclient.Entity(id='syn1234')
    entity.modifiedBy = '333'
    entity.createdBy = '333'
    error_email = (
        "Dear %s,\n\n"
        "Your files (%s) are duplicated!  FILES SHOULD BE UPLOADED AS "
        "NEW VERSIONS AND THE ENTIRE DATASET SHOULD BE "
        "UPLOADED EVERYTIME".format("trial", "first.cbs"))
    with mock.patch.object(
            syn, "get", return_value=entity) as patch_syn_get,\
        mock.patch.object(
            syn, "getUserProfile",
            return_value={'userName': 'trial'}) as patch_syn_profile,\
        mock.patch.object(
            syn, "sendMessage") as patch_send:
        input_to_database.email_duplication_error(syn, duplicated_filesdf)
        patch_syn_get.assert_called_once_with('syn1234')
        patch_syn_profile.assert_called_once_with('333')
        patch_send.assert_called_once_with(
            ['333'], "GENIE Validation Error", error_email)


def test_nodups_email_duplication_error():
    '''
    Test no email sent
    '''
    duplicated_filesdf = pd.DataFrame()
    with mock.patch.object(syn, "get") as patch_syn_get,\
            mock.patch.object(syn, "getUserProfile") as patch_syn_profile,\
            mock.patch.object(syn, "sendMessage") as patch_send:
        input_to_database.email_duplication_error(syn, duplicated_filesdf)
        patch_syn_get.assert_not_called()
        patch_syn_profile.assert_not_called()
        patch_send.assert_not_called()


def test_valid__get_status_and_error_list():
    '''
    Tests the correct status and error lists received
    when file is valid.
    '''
    entity = synapseclient.Entity(id='syn1234', md5='44444')
    entities = [entity]
    filetype = "clinical"
    fileinfo = {'filePaths': ['/path/to/data_clinical_supp_SAGE.txt'],
                'synId': ['syn1234']}
    modified_ons = [1553428800000]
    filenames = ['data_clinical_supp_SAGE.txt']

    valid = True
    message = 'valid'
    filetype = 'clinical'
    input_status_list, invalid_errors_list = \
        input_to_database._get_status_and_error_list(
           syn, fileinfo, valid, message, filetype,
           entities, fileinfo['filePaths'], filenames, modified_ons)
    assert input_status_list == [
        [entity.id, fileinfo['filePaths'][0], entity.md5,
         'VALIDATED', filenames[0], modified_ons[0],
         filetype]]
    assert invalid_errors_list is None


def test_invalid__get_status_and_error_list():
    '''
    Tests the correct status and error lists received
    when file is invalid.
    '''
    entity = synapseclient.Entity(id='syn1234', md5='44444')
    entities = [entity]
    filetype = "clinical"
    fileinfo = {'filePaths': ['/path/to/data_clinical_supp_SAGE.txt'],
                'synId': ['syn1234']}
    modified_ons = [1553428800000]
    filenames = ['data_clinical_supp_SAGE.txt']
    # This valid variable control the validation status
    valid = False
    message = 'invalid file content'
    filetype = 'clinical'

    input_status_list, invalid_errors_list = \
        input_to_database._get_status_and_error_list(
            syn, fileinfo, valid, message, filetype,
            entities, fileinfo['filePaths'], filenames, modified_ons)
    assert input_status_list == [
        [entity.id, fileinfo['filePaths'][0], entity.md5,
            'INVALID', filenames[0], modified_ons[0],
            filetype]]
    assert invalid_errors_list == [
        ['syn1234', message, 'data_clinical_supp_SAGE.txt']]


def test__send_validation_email():
    message = 'invalid error message here'
    filenames = ['data_clinical_supp_SAGE.txt']
    file_users = ['333', '444']
    with mock.patch.object(
        syn, "getUserProfile",
        return_value={'userName': 'trial'}) as patch_syn_getuserprofile,\
        mock.patch.object(
            syn, "sendMessage") as patch_syn_sendmessage:
        input_to_database._send_validation_error_email(
            syn, filenames, message, file_users)
        error_message = (
            "Dear trial, trial,\n\n"
            "Your files (data_clinical_supp_SAGE.txt) are invalid! "
            "Here are the reasons why:\n\n%s" % message)
        patch_syn_sendmessage.assert_called_once_with(
            ['333', '444'], "GENIE Validation Error", error_message)
        patch_syn_getuserprofile.call_count == 2
