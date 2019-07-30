import mock
import os
import pytest

import pandas as pd
import synapseclient
import synapseutils

from genie import input_to_database
from genie.clinical import clinical
from genie.mafSP import mafSP
from genie.maf import maf
from genie.vcf import vcf



syn = mock.create_autospec(synapseclient.Synapse)
sample_clinical_synid = 'syn2222'

sample_clinical_entity = synapseclient.File(path='data_clinical_supp_sample_SAGE.txt',
                                            id=sample_clinical_synid,
                                            parentId='syn45678',
                                            name='data_clinical_supp_sample_SAGE.txt')

patient_clinical_synid = 'syn11111'
patient_clinical_entity = synapseclient.File(path='data_clinical_supp_patient_SAGE.txt',
                                             id=patient_clinical_synid,
                                             parentId='syn45678',
                                             name='data_clinical_supp_patient_SAGE.txt')

vcf1synid = 'syn6666'
vcf1_entity = synapseclient.File(path='GENIE-SAGE-1-1.vcf',
                                 id=vcf1synid,
                                 parentId='syn45678',
                                 name='GENIE-SAGE-1-1.vcf')
vcf2synid = 'syn8888'
vcf2_entity = synapseclient.File(path='GENIE-SAGE-2-1.vcf',
                                 id=vcf2synid,
                                 parentId='syn45678',
                                 name='GENIE-SAGE-2-1.vcf')
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
oncotree_link = "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21"
center_input_synid = "syn9999"
center_staging_synid = "syn9999"
center_mapping = {'inputSynId': [center_input_synid],
                  'stagingSynId': [center_input_synid],
                  'center': [center]}
center_mapping_df = pd.DataFrame(center_mapping)
validation_statusdf = pd.DataFrame({
    'id': ['syn1234', 'syn2345'],
    'status': ['VALID', 'INVALID'],
    'md5': ['3333', '44444'],
    'name': ['first.txt', 'second.txt']})
error_trackerdf = pd.DataFrame({
    'id': ['syn2345'],
    'errors': ['Invalid file format']})
emptydf = pd.DataFrame(columns=['id'], dtype=str)
# def test_samename_rename_file():
#     '''Test that the file path is not renamed.
#     '''
#     filename = synapseclient.utils.make_bogus_data_file()
#     entity = synapseclient.File(path=filename,
#                                 id='syn012345',
#                                 parentId='syn45678',
#                                 name=os.path.basename(filename))
#     expectedpath = filename
#     new_entity = input_to_database.rename_file(entity)
#     assert new_entity.annotations.expectedPath == expectedpath
#     os.remove(filename)


# def test_diffname_rename_file():
#     '''Test that the file path is renamed.
#     '''
#     filename = synapseclient.utils.make_bogus_data_file()
#     entity = synapseclient.File(path=filename,
#                                 id='syn012345',
#                                 parentId='syn45678',
#                                 name='testname')

#     expectedpath = os.path.join(os.path.dirname(filename), "testname")
#     new_entity = input_to_database.rename_file(entity)
#     assert new_entity.annotations.expectedPath == expectedpath
#     os.remove(filename)


def walk_return():
    '''
    Generator returned by synapseutils.walk
    '''
    yield first
    yield second


def walk_return_empty():
    '''
    Generator returned by synapseutils.walk
    '''
    yield ([], [], [])


def test_main_get_center_input_files():
    '''
    Test to make sure center input files are gotten
    excluding the vcf files since process main is specified
    '''
    syn_get_effects = [sample_clinical_entity, patient_clinical_entity]
    expected_center_file_list = [syn_get_effects]

    calls = [mock.call(sample_clinical_synid),
             mock.call(patient_clinical_synid)]

    with mock.patch.object(synapseutils, "walk",
                           return_value=walk_return()) as patch_synapseutils_walk,\
        mock.patch.object(syn, "get",
                          side_effect=syn_get_effects) as patch_syn_get:
        center_file_list = input_to_database.get_center_input_files(syn,
                                                                    "syn12345",
                                                                    center)

        assert len(center_file_list) == len(expected_center_file_list)
        assert len(center_file_list[0]) == 2
        assert center_file_list == expected_center_file_list
        patch_synapseutils_walk.assert_called_once_with(syn, 'syn12345')
        patch_syn_get.assert_has_calls(calls)


def test_vcf_get_center_input_files():
    '''
    Test to make sure center input files are gotten
    including the vcf files since process vcf is specified
    '''
    syn_get_effects = [sample_clinical_entity, patient_clinical_entity,
                       vcf1_entity, vcf2_entity]
    expected_center_file_list = [
        [vcf1_entity], [vcf2_entity],
        [sample_clinical_entity, patient_clinical_entity]]
    calls = [
        mock.call(sample_clinical_synid),
        mock.call(patient_clinical_synid),
        mock.call(vcf1synid),
        mock.call(vcf2synid)]

    with mock.patch.object(synapseutils, "walk",
                           return_value=walk_return()) as patch_synapseutils_walk,\
        mock.patch.object(syn, "get",
                          side_effect=syn_get_effects) as patch_syn_get:
        center_file_list = input_to_database.get_center_input_files(syn,
                                                                    "syn12345",
                                                                    center,
                                                                    process="vcf")
        assert len(center_file_list) == len(expected_center_file_list)
        assert len(center_file_list[2]) == 2
        assert center_file_list == expected_center_file_list
        patch_synapseutils_walk.assert_called_once_with(syn, 'syn12345')
        patch_syn_get.assert_has_calls(calls)


def test_empty_get_center_input_files():
    '''
    Test that center input files is empty if directory
    pass in is empty
    '''
    with mock.patch.object(synapseutils, "walk",
                           return_value=walk_return_empty()) as patch_synapseutils_walk:
        center_file_list = input_to_database.get_center_input_files(
            syn, "syn12345", center, process="vcf")
        assert center_file_list == []
        patch_synapseutils_walk.assert_called_once_with(syn, 'syn12345')


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
    '''
    Test the values returned by input that hasn't be validated
    '''
    entity = synapseclient.Entity(id='syn1234')
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        emptydf, emptydf, entities)
    assert file_status['to_validate']
    assert file_status['status_list'] == []
    assert file_status['error_list'] == []


def test_valid_check_existing_file_status():
    '''
    Test the values returned by input that is already valid
    '''
    entity = synapseclient.Entity(name='first.txt', id='syn1234', md5='3333')
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities)
    assert not file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_invalid_check_existing_file_status():
    '''
    Test the values returned by input that is invalid
    '''
    entity = synapseclient.Entity(name='second.txt', id='syn2345', md5='44444')
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities)
    assert not file_status['to_validate']
    assert file_status['status_list'] == ['INVALID']
    assert file_status['error_list'] == ['Invalid file format']


def test_nostorederrors_check_existing_file_status():
    '''
    If there is no error uploaded, must re-validate file
    '''
    entity = synapseclient.Entity(name='second.txt', id='syn2345', md5='44444')
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, emptydf, entities)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['INVALID']
    assert file_status['error_list'] == []


def test_diffmd5validate_check_existing_file_status():
    '''
    If md5 is different from stored md5, must re-validate file
    '''
    entity = synapseclient.Entity(name='first.txt', id='syn1234', md5='44444')
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, emptydf, entities)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_diffnametovalidate_check_existing_file_status():
    '''
    If name is different from stored name, must re-validate file
    '''
    entity = synapseclient.Entity(name='second.txt', id='syn1234', md5='3333')
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, emptydf, entities)
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
    first_entity = synapseclient.Entity(name='first.txt', id='syn1234', md5='3333')
    second_entity = synapseclient.Entity(name='second.txt', id='syn2345', md5='44444')
    entities = [first_entity, second_entity]
    file_status = input_to_database.check_existing_file_status(
        validation_statusdf, error_trackerdf, entities)
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
        entities = ['foo', 'doo', 'boo']
        input_to_database.check_existing_file_status(
            emptydf, emptydf, entities)


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
    entity = synapseclient.Entity(name="data_clinical_supp_SAGE.txt",
                                  id='syn1234', md5='44444',
                                  path='/path/to/data_clinical_supp_SAGE.txt')
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = '333'
    entity.createdBy = '444'
    entities = [entity]
    threads = 0
    testing = False
    valid = True
    message = "Is valid"
    filetype = "clinical"

    expected_validate_results = ([[
        entity.id,
        entity.path,
        '44444',
        'VALIDATED',
        'data_clinical_supp_SAGE.txt',
        1553428800000,
        'clinical']], None)
    with mock.patch(
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
            return_value=(valid, message, filetype)) as patch_validate,\
        mock.patch(
            "genie.input_to_database._get_status_and_error_list",
            return_value=expected_validate_results) as patch_get_staterror_list,\
        mock.patch("genie.input_to_database._send_validation_error_email") as patch_send_email:

        validate_results = input_to_database.validatefile(
            syn, entities, validation_statusdf,
            error_trackerdf, center, threads, testing, oncotree_link)

        assert expected_validate_results == validate_results
        patch_validate.assert_called_once_with(
            syn,
            [entity.path],
            center,
            filetype=filetype,
            oncotreelink=oncotree_link,
            testing=testing
        )
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities)
        patch_determine_filetype.assert_called_once_with(
            syn, [entity.name], center)
        patch_get_staterror_list.assert_called_once_with(
            syn, valid, message, filetype,
            entities)
        patch_send_email.assert_not_called()


def test_invalid_validatefile():
    '''
    Tests the behavior of a file that gets validated that becomes
    invalid
    '''
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame(columns=['id'], dtype=str)
    entity = synapseclient.Entity(name="data_clinical_supp_SAGE.txt",
                                  id='syn1234', md5='44444',
                                  path='/path/to/data_clinical_supp_SAGE.txt')
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = '333'
    entity.createdBy = '444'
    entities = [entity]
    threads = 0
    testing = False
    valid = False
    message = "Is invalid"
    filetype = "clinical"
    # fileinfo = {'filePaths': ['/path/to/data_clinical_supp_SAGE.txt'],
    #             'synId': ['syn1234']}
    expected_validate_results = ([[
        entity.id, entity.path,
        '44444', 'INVALID',
        entity.name, 1553428800000, 'clinical']],
        [entity.id, message, entity.name])
    with mock.patch(
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
            return_value=(valid, message, filetype)) as patch_validate,\
        mock.patch(
            "genie.input_to_database._get_status_and_error_list",
            return_value=expected_validate_results) as patch_get_staterror_list,\
        mock.patch("genie.input_to_database._send_validation_error_email") as patch_send_email:

        validate_results = input_to_database.validatefile(
            syn, entities, validation_statusdf,
            error_trackerdf, center, threads, testing, oncotree_link)

        assert expected_validate_results == validate_results
        patch_validate.assert_called_once_with(
            syn,
            [entity.path],
            center,
            filetype=filetype,
            oncotreelink=oncotree_link,
            testing=testing
        )
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities)
        patch_determine_filetype.assert_called_once_with(
            syn, [entity.name], center)
        patch_get_staterror_list.assert_called_once_with(
            syn, valid, message, filetype,
            entities)
        patch_send_email.assert_called_once_with(
            syn, [entity.name], message, [entity.modifiedBy, entity.createdBy])


def test_already_validated_validatefile():
    '''
    Test already validated files
    '''
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame()
    entity = synapseclient.Entity(name="data_clinical_supp_SAGE.txt",
                                  id='syn1234', md5='44444',
                                  path='/path/to/data_clinical_supp_SAGE.txt')
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = '333'
    entity.createdBy = '444'
    entities = [entity]
    threads = 0
    testing = False
    valid = True
    message = "Is valid"
    filetype = "markdown"

    check_file_status_dict = {
        'status_list': ["INVALID"],
        'error_list': ["invalid file"],
        'to_validate': False}
    expected_validate_results = (
        [[entity.id,
            entity.path,
            entity.md5,
            check_file_status_dict['status_list'][0],
            entity.name,
            1553428800000,
            filetype]],
        [[entity.id,
            check_file_status_dict['error_list'][0],
            entity.name]])
    with mock.patch(
            "genie.validate.determine_filetype",
            return_value=filetype) as patch_determine_filetype,\
        mock.patch(
            "genie.input_to_database.check_existing_file_status",
            return_value=check_file_status_dict) as patch_check, \
        mock.patch(
            "genie.validate.validate_single_file",
            return_value=(valid, message, filetype)) as patch_validate,\
        mock.patch(
            "genie.input_to_database._get_status_and_error_list",
            return_value=expected_validate_results) as patch_get_staterror_list,\
        mock.patch("genie.input_to_database._send_validation_error_email") as patch_send_email:

        validate_results = input_to_database.validatefile(
            syn, entities, validation_statusdf,
            error_trackerdf, center, threads, testing, oncotree_link)

        assert expected_validate_results == validate_results
        patch_validate.assert_not_called()
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities)
        patch_determine_filetype.assert_called_once_with(
            syn, [entity.name], center)
        patch_get_staterror_list.assert_not_called()
        patch_send_email.assert_not_called()


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
    dupsdf = input_to_database.get_duplicated_files(validation_statusdf, "")
    assert dupsdf.equals(expected_dup)


def test_nodups_get_duplicated_files():
    '''
    Test no duplicated
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345', 'syn5555', 'syn1224', 'syn34444'],
        'name': ['cbs.txt', 'second.seg', 'no_clinical.txt',
                 'data_clinical_supp_2', 'data_clinical_supp_3']})
    dupsdf = input_to_database.get_duplicated_files(validation_statusdf, "")
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
        "Dear {},\n\n"
        "Your files ({}) are duplicated!  FILES SHOULD BE UPLOADED AS "
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
    modified_on = 1561143558000
    modified_on_string = "2019-06-21T18:59:18.456Z"

    entity = synapseclient.Entity(id='syn1234', md5='44444',
                                  path='/path/to/foobar.txt',
                                  name='data_clinical_supp_SAGE.txt')
    entity.properties.modifiedOn = modified_on_string

    entities = [entity]

    valid = True
    message = 'valid'
    filetype = 'clinical'

    input_status_list, invalid_errors_list = \
        input_to_database._get_status_and_error_list(
           syn, valid, message, filetype,
           entities)
    assert input_status_list == [
        [entity.id, entity.path, entity.md5,
         'VALIDATED', entity.name, modified_on,
         filetype]]
    assert invalid_errors_list is None


def test_invalid__get_status_and_error_list():
    '''
    Tests the correct status and error lists received
    when file is invalid.
    '''
    modified_on = 1561143558000
    modified_on_string = "2019-06-21T18:59:18.456Z"
    entity = synapseclient.Entity(id='syn1234', md5='44444',
                                  path='/path/to/foobar.txt',
                                  name='data_clinical_supp_SAGE.txt')
    entity.properties.modifiedOn = modified_on_string

    entities = [entity]
    filetype = "clinical"
    # This valid variable control the validation status
    valid = False
    message = 'invalid file content'

    input_status_list, invalid_errors_list = \
        input_to_database._get_status_and_error_list(
            syn, valid, message, filetype,
            entities)
    assert input_status_list == [
        [entity.id, entity.path, entity.md5,
            'INVALID', entity.name, modified_on,
            filetype]]
    assert invalid_errors_list == [
        ['syn1234', message, 'data_clinical_supp_SAGE.txt']]


def test__send_validation_error_email():
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


class emptytable_mock:
    '''
    Validation status tablequery dataframe mocking and
    error tracking tablequery dataframe mocking
    This is used because assert_called_once_with has a hard
    time with comparing pandas dataframes
    '''
    tableId = "syn555"

    def asDataFrame(self):
        return([])


def test_update_status_and_error_tables():
    '''
    Test updating validation status and error table
    '''
    validation_status_table = emptytable_mock()
    error_tracker_table = emptytable_mock()
    input_valid_statuses = [[
        sample_clinical_entity.id,
        sample_clinical_entity.path,
        '44444',
        'VALIDATED',
        'data_clinical_supp_SAGE.txt',
        1553428800000,
        'clinical']]
    invalid_errors = []
    input_valid_statusdf = pd.DataFrame(input_valid_statuses,
                                        columns=["id", 'path', 'md5', 'status',
                                                 'name', 'modifiedOn',
                                                 'fileType'])
    input_valid_statusdf['center'] = center
    with mock.patch(
            "genie.input_to_database.get_duplicated_files",
            return_value=pd.DataFrame(columns=['id', 'errors', 'name'], dtype=str)) as mock_get_duplicated,\
        mock.patch(
            "genie.input_to_database.email_duplication_error") as mock_email,\
        mock.patch(
            "genie.process_functions.updateDatabase") as mock_update:
        input_validdf = input_to_database.update_status_and_error_tables(
            syn,
            center,
            input_valid_statuses,
            invalid_errors,
            validation_status_table,
            error_tracker_table)
        mock_get_duplicated.assert_called_once()
        mock_email.assert_not_called()
        assert mock_update.call_count == 2
        assert input_validdf.equals(input_valid_statusdf)


def test_validation():
    '''
    Test validation steps
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234'],
        'status': ['VALIDATED'],
        'path': ["/path/to/file"],
        'fileType': ['clinical']})

    thread = 2
    testing = False
    modified_on = 1561143558000
    process = "main"
    databaseToSynIdMapping = {'Database': ["clinical", 'validationStatus', 'errorTracker'],
                              'Id': ['syn222', 'syn333', 'syn444']}
    databaseToSynIdMappingDf = pd.DataFrame(databaseToSynIdMapping)
    entity = synapseclient.Entity(id='syn1234', md5='44444',
                                  path='/path/to/foobar.txt',
                                  name='data_clinical_supp_SAGE.txt')
    entities = [entity]
    filetype = "clinical"
    input_status_list = [
        [entity.id, entity.path, entity.md5,
         'VALIDATED', entity.name, modified_on,
         filetype]]
    invalid_errors_list = None
    validationstatus_mock = emptytable_mock()
    errortracking_mock = emptytable_mock()
    with mock.patch(
            "genie.input_to_database.get_center_input_files",
            return_value=entities) as patch_get_center,\
        mock.patch.object(
            syn, "tableQuery",
            side_effect=[validationstatus_mock,
                         errortracking_mock]) as patch_tablequery,\
        mock.patch(
            "genie.input_to_database.validatefile",
            return_value=(input_status_list, invalid_errors_list)) as patch_validatefile,\
        mock.patch(
            "genie.input_to_database.update_status_and_error_tables",
            return_value=validation_statusdf) as patch_update_status:
        valid_filedf = input_to_database.validation(
            syn, center, process,
            center_mapping_df, databaseToSynIdMappingDf,
            thread, testing, oncotree_link)
        patch_get_center.assert_called_once_with(
            syn, center_input_synid, center, process)
        assert patch_tablequery.call_count == 2
        patch_validatefile.assert_called_once_with(
            syn, entity,
            validationstatus_mock.asDataFrame(),
            errortracking_mock.asDataFrame(),
            center='SAGE', threads=1,
            testing=False,
            oncotree_link=oncotree_link)
        patch_update_status.assert_called_once_with(
            syn,
            center,
            input_status_list,
            [],
            validationstatus_mock,
            errortracking_mock)

        assert valid_filedf.equals(validation_statusdf[['id', 'path', 'fileType']])


@pytest.mark.parametrize(
    'process, genieclass, filetype', [
        ('main', clinical, 'clinical'),
        ('maf', maf, 'maf'),
        ('mafSP', mafSP, 'mafSP')
    ]
)
def test_main_processfile(process, genieclass, filetype):
    validfiles = {'id': ['syn1'],
                  'path': ['/path/to/data_clinical_supp_SAGE.txt'],
                  'fileType': [filetype]}
    validfilesdf = pd.DataFrame(validfiles)
    center = "SAGE"
    path_to_genie = "./"
    threads = 2
    oncotreeLink = "www.google.com"
    center_mapping = {'stagingSynId': ["syn123"],
                      'center': [center]}
    center_mapping_df = pd.DataFrame(center_mapping)
    databaseToSynIdMapping = {'Database': [filetype],
                              'Id': ['syn222']}
    databaseToSynIdMappingDf = pd.DataFrame(databaseToSynIdMapping)

    with mock.patch.object(genieclass, "process") as patch_class:
        input_to_database.processfiles(
            syn, validfilesdf, center, path_to_genie, threads,
            center_mapping_df, oncotreeLink, databaseToSynIdMappingDf,
            validVCF=None, vcf2mafPath=None,
            veppath=None, vepdata=None,
            processing=process, test=False, reference=None)
        patch_class.assert_called_once()


def test_mainnone_processfile():
    '''
    If file type is None, the processing function is not called
    '''
    validfiles = {'id': ['syn1'],
                  'path': ['/path/to/data_clinical_supp_SAGE.txt'],
                  'fileType': [None]}
    validfilesdf = pd.DataFrame(validfiles)
    center = "SAGE"
    path_to_genie = "./"
    threads = 2
    oncotreeLink = "www.google.com"
    center_mapping = {'stagingSynId': ["syn123"],
                      'center': [center]}
    center_mapping_df = pd.DataFrame(center_mapping)
    databaseToSynIdMapping = {'Database': ["clinical"],
                              'Id': ['syn222']}
    databaseToSynIdMappingDf = pd.DataFrame(databaseToSynIdMapping)

    with mock.patch.object(clinical, "process") as patch_clin:
        input_to_database.processfiles(
            syn, validfilesdf, center, path_to_genie, threads,
            center_mapping_df, oncotreeLink, databaseToSynIdMappingDf,
            validVCF=None, vcf2mafPath=None,
            veppath=None, vepdata=None,
            processing="main", test=False, reference=None)
        patch_clin.assert_not_called()


def test_notvcf_processfile():
    '''
    Make sure vcf, maf, mafSP is called correctly
    '''
    validfiles = {'id': ['syn1'],
                  'path': ['/path/to/data_clinical_supp_SAGE.txt'],
                  'fileType': [None]}
    validfilesdf = pd.DataFrame(validfiles)
    center = "SAGE"
    path_to_genie = "./"
    threads = 2
    oncotreeLink = "www.google.com"
    center_mapping = {'stagingSynId': ["syn123"],
                      'center': [center]}
    center_mapping_df = pd.DataFrame(center_mapping)
    databaseToSynIdMapping = {'Database': ['vcf'],
                              'Id': ['syn222']}
    databaseToSynIdMappingDf = pd.DataFrame(databaseToSynIdMapping)

    with mock.patch.object(vcf, "process") as patch_process:
        input_to_database.processfiles(
            syn, validfilesdf, center, path_to_genie, threads,
            center_mapping_df, oncotreeLink, databaseToSynIdMappingDf,
            validVCF=None, vcf2mafPath=None,
            veppath=None, vepdata=None,
            processing='vcf', test=False, reference=None)
        patch_process.assert_called_once()
