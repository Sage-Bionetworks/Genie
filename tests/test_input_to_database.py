import mock
from mock import patch
import os
import pytest

import pandas as pd
import synapseclient
import synapseutils

import genie.config
from genie import input_to_database, process_functions
from genie.clinical import clinical
import genie.config
from genie.mafSP import mafSP
from genie.maf import maf
from genie.vcf import vcf
from genie.validate import GenieValidationHelper


syn = mock.create_autospec(synapseclient.Synapse)
sample_clinical_synid = 'syn2222'

sample_clinical_entity = synapseclient.File(path='data_clinical_supp_sample_SAGE.txt',
                                            id=sample_clinical_synid,
                                            parentId='syn45678',
                                            name='data_clinical_supp_sample_SAGE.txt',
                                            modifiedOn = '2019-03-24T12:00:00.Z',
                                            md5='44444')

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
    'name': ['first.txt', 'second.txt'],
    'fileType': ['filetype1', 'filetype2']})
error_trackerdf = pd.DataFrame({
    'id': ['syn2345'],
    'errors': ['Invalid file format'],
    'fileType': ['filetype1']})
emptydf = pd.DataFrame(columns=['id'], dtype=str)

class mock_csv_query_result(object):
    def __init__(self, df):
        self.df = df
    def asDataFrame(self):
        return self.df

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
    '''Test to make sure center input files are retrieved.
    '''
    syn_get_effects = [sample_clinical_entity, patient_clinical_entity,
                       vcf1_entity, vcf2_entity]
    expected_center_file_list = [
        [vcf1_entity], [vcf2_entity],
        [sample_clinical_entity, patient_clinical_entity]]
    calls = [
        mock.call(sample_clinical_synid, downloadFile=True),
        mock.call(patient_clinical_synid, downloadFile=True),
        mock.call(vcf1synid, downloadFile=True),
        mock.call(vcf2synid, downloadFile=True)]

    with patch.object(synapseutils, "walk",
                      return_value=walk_return()) as patch_synapseutils_walk,\
         patch.object(syn, "get",
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
    with patch.object(synapseutils, "walk",
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
    entity.properties.versionNumber = '1'

    entities = [entity]

    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(emptydf), mock_csv_query_result(emptydf), entities)
    assert file_status['to_validate']
    assert file_status['status_list'] == []
    assert file_status['error_list'] == []


def test_valid_check_existing_file_status():
    '''
    Test the values returned by input that is already valid
    '''
    entity = synapseclient.Entity(name='first.txt', id='syn1234', md5='3333')
    entity.properties.versionNumber = '1'

    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf), mock_csv_query_result(error_trackerdf), entities)
    assert not file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_invalid_check_existing_file_status():
    '''
    Test the values returned by input that is invalid
    '''
    entity = synapseclient.Entity(name='second.txt', id='syn2345', md5='44444')
    entity.properties.versionNumber = '1'
    entities = [entity]

    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf), mock_csv_query_result(error_trackerdf), entities)
    assert not file_status['to_validate']
    assert file_status['status_list'] == ['INVALID']
    assert file_status['error_list'] == ['Invalid file format']


def test_nostorederrors_check_existing_file_status():
    '''
    If there is no error uploaded, must re-validate file
    '''
    entity = synapseclient.Entity(name='second.txt', id='syn2345', md5='44444')
    entity.properties.versionNumber = '1'
    entities = [entity]

    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf), mock_csv_query_result(emptydf), entities)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['INVALID']
    assert file_status['error_list'] == []


def test_diffmd5validate_check_existing_file_status():
    '''
    If md5 is different from stored md5, must re-validate file
    '''
    entity = synapseclient.Entity(name='first.txt', id='syn1234', md5='44444')
    entity.properties.versionNumber = '1'
    entities = [entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf), mock_csv_query_result(emptydf), entities)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_diffnametovalidate_check_existing_file_status():
    '''
    If name is different from stored name, must re-validate file
    '''
    entity = synapseclient.Entity(name='second.txt', id='syn1234', md5='3333')
    entity.properties.versionNumber = '1'
    entities = [entity]

    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf), mock_csv_query_result(emptydf), entities)
    assert file_status['to_validate']
    assert file_status['status_list'] == ['VALID']
    assert file_status['error_list'] == []


def test_twoinvalid_check_existing_file_status():
    '''
    Test to make sure both invalid statues get passed back
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'versionNumber': ['1', '1'],
        'status': ['INVALID', 'INVALID'],
        'md5': ['3333', '44444'],
        'name': ['first.txt', 'second.txt']})
    error_trackerdf = pd.DataFrame({
        'id': ['syn1234', 'syn2345'],
        'versionNumber': ['1', '1'],
        'errors': ['Invalid file format', 'Invalid formatting issues']})

    first_entity = synapseclient.Entity(name='first.txt', id='syn1234', md5='3333')
    first_entity.properties.versionNumber = '1'

    second_entity = synapseclient.Entity(name='second.txt', id='syn2345', md5='44444')
    second_entity.properties.versionNumber = '1'

    entities = [first_entity, second_entity]
    file_status = input_to_database.check_existing_file_status(
        mock_csv_query_result(validation_statusdf), mock_csv_query_result(error_trackerdf), entities)
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

    with patch.object(syn, "store",
                      return_value=new_maf_ent) as patch_syn_store,\
         patch.object(syn, "setPermissions",
                      return_value=None) as patch_syn_set_permissions,\
         patch.object(syn, "get",
                      return_value=table_ent) as patch_syn_get,\
         patch.object(syn, "getTableColumns",
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
    entity.properties.versionNumber = '1'
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = '333'
    entity.createdBy = '444'
    entities = [entity]
    valid = True
    message = "Is valid"
    filetype = "clinical"
    # Only a list is returned as opposed a list of lists when there are
    # invalid errors
    status_error_list_results = ([{'entity': entity, 'status': 'VALIDATED'}],
                                 [])
    expected_results = ([{'entity': entity, 'status': 'VALIDATED',
                          'fileType': filetype, 'center': center}],
                        [], [])
    with patch.object(GenieValidationHelper, "determine_filetype",
                      return_value=filetype) as patch_determine_filetype,\
         patch.object(input_to_database, "check_existing_file_status",
                      return_value={'status_list': [],
                                    'error_list': [],
                                    'to_validate': True}) as patch_check, \
         patch.object(GenieValidationHelper,"validate_single_file",
                      return_value=(valid, message,
                                    filetype)) as patch_validate,\
         patch.object(input_to_database, "_get_status_and_error_list",
                      return_value=status_error_list_results) as patch_get_staterror_list,\
         patch.object(input_to_database,
                      "_send_validation_error_email") as patch_send_email:

        validate_results = input_to_database.validatefile(
            syn, None, entities, validation_statusdf,
            error_trackerdf, center, oncotree_link,
            format_registry=genie.config.PROCESS_FILES)

        assert expected_results == validate_results
        patch_validate.assert_called_once_with(
            oncotree_link=oncotree_link, nosymbol_check=False)
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities)
        patch_determine_filetype.assert_called_once()
        patch_get_staterror_list.assert_called_once_with(
            valid, message, entities)
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
    valid = False
    message = "Is invalid"
    filetype = "clinical"
    status_error_list_results = ([{'entity': entity, 'status': 'INVALID'}],
                                 [{'entity': entity, 'errors': message}])
    expected_results = ([{'entity': entity, 'status': 'INVALID',
                          'fileType': filetype, 'center': center}],
                        [{'entity': entity, 'errors': message,
                          'fileType': filetype, 'center': center}],
                          [(['data_clinical_supp_SAGE.txt'], 'Is invalid', ['333', '444'])])

    with patch.object(GenieValidationHelper, "determine_filetype",
                      return_value=filetype) as patch_determine_filetype,\
         patch.object(input_to_database, "check_existing_file_status",
                      return_value={'status_list': [],
                                    'error_list': [],
                                    'to_validate': True}) as patch_check, \
         patch.object(GenieValidationHelper, "validate_single_file",
                      return_value=(valid, message,
                                    filetype)) as patch_validate,\
         patch.object(input_to_database, "_get_status_and_error_list",
                      return_value=status_error_list_results) as patch_get_staterror_list:

        validate_results = input_to_database.validatefile(
            syn, None, entities, validation_statusdf,
            error_trackerdf, center, oncotree_link,
            format_registry=genie.config.PROCESS_FILES)

        assert expected_results == validate_results
        patch_validate.assert_called_once_with(
            oncotree_link=oncotree_link, nosymbol_check=False)
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities)
        patch_determine_filetype.assert_called_once()
        patch_get_staterror_list.assert_called_once_with(
            valid, message, entities)


def test_already_validated_validatefile():
    '''
    Test already validated files
    '''
    validation_statusdf = pd.DataFrame()
    error_trackerdf = pd.DataFrame()
    entity = synapseclient.Entity(name="data_clinical_supp_SAGE.txt",
                                  id='syn1234', md5='44444',
                                  path='/path/to/data_clinical_supp_SAGE.txt')
    entity.properties.versionNumber = '1'
    entity['modifiedOn'] = '2019-03-24T12:00:00.Z'
    # This modifiedOn translates to: 1553428800000
    entity.modifiedBy = '333'
    entity.createdBy = '444'
    entities = [entity]
    valid = False
    errors = "Invalid file"
    filetype = "markdown"
    status = "INVALID"
    check_file_status_dict = {
        'status_list': [status],
        'error_list': [errors],
        'to_validate': False}

    status_error_list_results = ([{'entity': entity, 'status': status}],
                                 [{'entity': entity, 'errors': errors}])
    expected_results = ([{'entity': entity, 'status': status,
                          'fileType': filetype, 'center': center}],
                        [{'entity': entity, 'errors': errors,
                          'fileType': filetype, 'center': center}],
                        [])
    with patch.object(GenieValidationHelper, "determine_filetype",
                      return_value=filetype) as patch_determine_filetype,\
         patch.object(input_to_database, "check_existing_file_status",
                      return_value=check_file_status_dict) as patch_check, \
         patch.object(GenieValidationHelper, "validate_single_file",
                      return_value=(valid, errors, filetype)) as patch_validate,\
         patch.object(input_to_database, "_get_status_and_error_list",
                      return_value=status_error_list_results) as patch_get_staterror_list,\
         patch.object(input_to_database,
                      "_send_validation_error_email") as patch_send_email:

        validate_results = input_to_database.validatefile(
            syn, None, entities, validation_statusdf,
            error_trackerdf, center, oncotree_link,
            format_registry=genie.config.PROCESS_FILES)

        assert expected_results == validate_results

        patch_validate.assert_not_called()
        patch_check.assert_called_once_with(
            validation_statusdf, error_trackerdf, entities)
        patch_determine_filetype.assert_called_once()
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
    with patch.object(syn, "get", return_value=entity) as patch_syn_get,\
         patch.object(syn, "getUserProfile",
                      return_value={'userName':
                                    'trial'}) as patch_syn_profile,\
         patch.object(syn, "sendMessage") as patch_send:
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
    with patch.object(syn, "get") as patch_syn_get,\
         patch.object(syn, "getUserProfile") as patch_syn_profile,\
         patch.object(syn, "sendMessage") as patch_send:
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
    entity.properties.versionNumber = '1'
    entity.properties.modifiedOn = modified_on_string

    entities = [entity]

    valid = True
    message = 'valid'
    filetype = 'clinical'

    input_status_list, invalid_errors_list = \
        input_to_database._get_status_and_error_list(
           valid, message, entities)
    assert input_status_list == [{'entity': entity, 'status': 'VALIDATED'}]
    assert not invalid_errors_list


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
    entity.properties.versionNumber = '1'
    entity.properties.modifiedOn = modified_on_string

    entities = [entity]
    filetype = "clinical"
    # This valid variable control the validation status
    valid = False
    errors = 'invalid file content'

    input_status_list, invalid_errors_list = \
        input_to_database._get_status_and_error_list(
            valid, errors, entities)
    assert input_status_list == [{'entity': entity, 'status': 'INVALID'}]
    assert invalid_errors_list == [{'entity': entity, 'errors': errors}]

from datetime import datetime
def test__send_validation_error_email():
    message = "invalid error message here"
    filenames = ['data_clinical_supp_SAGE.txt']
    file_user = '333'
    message_objs = [dict(filenames=filenames, messages=message)]
    print(message_objs)
    with patch.object(syn, "getUserProfile",
                      return_value={'userName':
                                    'trial'}) as patch_syn_getuserprofile,\
         patch.object(syn, "sendMessage") as patch_syn_sendmessage,\
         patch('genie.input_to_database.datetime') as mock_datetime:
        mock_datetime.datetime.today.return_value = datetime(2019, 11, 1, 00, 00, 00, 0)
        mock_datetime.side_effect = lambda *args, **kw: date(*args, **kw)        

        input_to_database._send_validation_error_email(
            syn, file_user, message_objs)
        error_message = (
            "Dear trial,\n\n"
            "You have invalid files! "
            "Here are the reasons why:\n\n"
            "Filenames: data_clinical_supp_SAGE.txt, Errors: %s\n" % message)
        patch_syn_sendmessage.assert_called_once_with(
            userIds=[file_user], messageBody=error_message,
            messageSubject='GENIE Validation Error - 2019-11-01 00:00:00')
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

    input_valid_statuses = [{'entity': sample_clinical_entity,
                             'status': 'VALIDATED',
                             'fileType': 'clinical',
                             'center': center}]

    expected = [[sample_clinical_entity.id, sample_clinical_entity.path,
                 sample_clinical_entity.md5, 'VALIDATED',
                 sample_clinical_entity.name,
                 1553428800000, 'clinical', center]]
    invalid_errors = []
    expected_statusdf = pd.DataFrame(expected,
                                     columns=["id", 'path', 'md5', 'status',
                                              'name', 'modifiedOn',
                                              'fileType', 'center'])
    #input_valid_statusdf['center'] = center
    empty_errorsdf = pd.DataFrame(columns=['id', 'errors', 'name',
                                           'fileType', 'center'], dtype=str)
    with patch.object(input_to_database, "get_duplicated_files",
                      return_value=empty_errorsdf) as mock_get_duplicated,\
         patch.object(input_to_database,
                      "email_duplication_error") as mock_email,\
         patch.object(process_functions, "updateDatabase") as mock_update:
        input_validdf = input_to_database.update_status_and_error_tables(
            syn,
            input_valid_statuses,
            invalid_errors,
            validation_status_table,
            error_tracker_table)
        mock_get_duplicated.assert_called_once()
        mock_email.assert_not_called()
        assert mock_update.call_count == 2
        assert input_validdf.equals(expected_statusdf[input_validdf.columns])
        assert expected_statusdf.equals(input_validdf[expected_statusdf.columns])


def test_validation():
    '''
    Test validation steps
    '''
    validation_statusdf = pd.DataFrame({
        'id': ['syn1234'],
        'status': ['VALIDATED'],
        'path': ["/path/to/file"],
        'fileType': ['clinical']})

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
         filetype, center]]
    invalid_errors_list = []
    messages = []
    validationstatus_mock = emptytable_mock()
    errortracking_mock = emptytable_mock()
    with patch.object(input_to_database, "get_center_input_files",
                      return_value=entities) as patch_get_center,\
         patch.object(syn, "tableQuery",
                      side_effect=[validationstatus_mock,
                                   errortracking_mock]) as patch_tablequery,\
         patch.object(input_to_database, "validatefile",
                      return_value=(input_status_list,
                                    invalid_errors_list, 
                                    messages)) as patch_validatefile,\
         patch.object(input_to_database, "update_status_and_error_tables",
                      return_value=validation_statusdf) as patch_update_status:
        valid_filedf = input_to_database.validation(
            syn, None, center, process,
            center_mapping_df, databaseToSynIdMappingDf,
            oncotree_link, genie.config.PROCESS_FILES)
        patch_get_center.assert_called_once_with(
            syn, center_input_synid, center, process)
        assert patch_tablequery.call_count == 2
        patch_validatefile.assert_called_once_with(
            syn, None, entity,
            validationstatus_mock,
            errortracking_mock,
            center='SAGE',
            oncotree_link=oncotree_link,
            format_registry=genie.config.PROCESS_FILES)
        patch_update_status.assert_called_once_with(
            syn,
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
    oncotree_link = "www.google.com"
    center_mapping = {'stagingSynId': ["syn123"],
                      'center': [center]}
    center_mapping_df = pd.DataFrame(center_mapping)
    databaseToSynIdMapping = {'Database': [filetype],
                              'Id': ['syn222']}
    databaseToSynIdMappingDf = pd.DataFrame(databaseToSynIdMapping)

    with patch.object(genieclass, "process") as patch_class:
        input_to_database.processfiles(
            syn, validfilesdf, center, path_to_genie,
            center_mapping_df, oncotree_link, databaseToSynIdMappingDf,
            validVCF=None, vcf2mafPath=None,
            veppath=None, vepdata=None,
            processing=process, reference=None)
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
    oncotree_link = "www.google.com"
    center_mapping = {'stagingSynId': ["syn123"],
                      'center': [center]}
    center_mapping_df = pd.DataFrame(center_mapping)
    databaseToSynIdMapping = {'Database': ["clinical"],
                              'Id': ['syn222']}
    databaseToSynIdMappingDf = pd.DataFrame(databaseToSynIdMapping)

    with patch.object(clinical, "process") as patch_clin:
        input_to_database.processfiles(
            syn, validfilesdf, center, path_to_genie,
            center_mapping_df, oncotree_link, databaseToSynIdMappingDf,
            validVCF=None, vcf2mafPath=None,
            veppath=None, vepdata=None,
            processing="main", reference=None)
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
    oncotree_link = "www.google.com"
    center_mapping = {'stagingSynId': ["syn123"],
                      'center': [center]}
    center_mapping_df = pd.DataFrame(center_mapping)
    databaseToSynIdMapping = {'Database': ['vcf'],
                              'Id': ['syn222']}
    databaseToSynIdMappingDf = pd.DataFrame(databaseToSynIdMapping)

    with patch.object(vcf, "process") as patch_process:
        input_to_database.processfiles(
            syn, validfilesdf, center, path_to_genie,
            center_mapping_df, oncotree_link, databaseToSynIdMappingDf,
            validVCF=None, vcf2mafPath=None,
            veppath=None, vepdata=None,
            processing='vcf', reference=None)
        patch_process.assert_called_once()
