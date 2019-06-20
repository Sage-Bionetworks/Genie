#!/usr/bin/env python3
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import synapseclient
import synapseutils
import os
import pandas as pd
import datetime
import shutil
import time
# Configuration file
from genie import PROCESS_FILES
from genie import process_functions
from genie import validate
from genie import toRetract
from genie import input_to_database

'''
TODO:
Could potentially get all the inforamation of the file entity right here
To avoid the syn.get rest call later which doesn't actually download the file
'''


def rename_file(syn, ent):
    '''
    Gets file from synapse and renames the file if necessary.

    Adds the expected name as an annotation to a Synapse File object.

    Args:
        syn: Synapse object
        synid : Synapse id or entity

    Returns:
        entity with annotation set for path of corrected file
    '''
    dirpath = os.path.dirname(ent.path)
    expectedpath = os.path.join(dirpath, ent.name)

    ent.annotations.expectedPath = expectedpath
    return ent


def get_center_input_files(syn, synid, center, process="main"):
    '''
    This function walks through each center's input directory
    to get a list of tuples of center files

    Args:
        syn: Synapse object
        synid: Center input folder synid
        center: Center name
        process: Process type includes, main, vcf, maf and mafSP.
                 Defaults to main such that the vcf

    Returns:
        List of entities with the correct format to pass into validation
    '''
    logger.info("GETTING {center} INPUT FILES".format(center=center))
    clinical_pair_name = [
        "data_clinical_supp_sample_{center}.txt".format(center=center),
        "data_clinical_supp_patient_{center}.txt".format(center=center)]

    center_files = synapseutils.walk(syn, synid)
    clinicalpair = []
    prepared_center_file_list = []
    
    for _, _, entities in center_files:
        for name, ent_synid in entities:
            # This is to remove vcfs from being validated during main 
            # processing. Often there are too many vcf files, and it is 
            # not necessary for them to be run everytime.
            if name.endswith(".vcf") and process != "vcf":
                continue

            ent = syn.get(ent_synid)
            logger.debug(ent)

            # Clinical file can come as two files.
            # The two files need to be merged together which is
            # why there is this format

            if name in clinical_pair_name:
                clinicalpair.append(ent)
                continue

            prepared_center_file_list.append([rename_file(syn, ent)])
    
    if clinicalpair:
        clinicalpair_entities = [rename_file(syn, x) for x in clinicalpair] 
        prepared_center_file_list.append(clinicalpair_entities)

    return prepared_center_file_list


def check_existing_file_status(validation_statusdf, error_trackerdf, entities):
    '''
    This function checks input files against the existing validation and error
    tracking dataframe

    Args:
        validation_statusdf: Validation status dataframe
        error_trackerdf: Error tracking dataframe
        entities: list of center input entites

    Returns:
        dict: Input file status
            status_list: file validation status
            error_list: Errors of the files if they exist,
            to_validate: Boolean value for whether of not an input
                         file needs to be validated
    '''
    if len(entities) > 2:
        raise ValueError(
            "There should never be more than 2 files being validated.")

    statuses = []
    errors = []

    for ent in entities:
        to_validate = False
        # Get the current status and errors from the tables.
        current_status = validation_statusdf[validation_statusdf['id'] == ent.id]
        current_error = error_trackerdf[error_trackerdf['id'] == ent.id]

        if current_status.empty:
            to_validate = True
        else:
            # This to_validate is here, because the following is a
            # sequential check of whether files need to be validated
            statuses.append(current_status['status'].values[0])
            if current_error.empty:
                to_validate = \
                    current_status['status'].values[0] == "INVALID"
            else:
                errors.append(current_error['errors'].values[0])
            # Add Name check here (must add name of the entity as a column)
            if current_status['md5'].values[0] != ent.md5 or \
               current_status['name'].values[0] != ent.name:
                to_validate = True
            else:
                status_str = "{filename} ({id}) FILE STATUS IS: {filestatus}"
                logger.info(status_str.format(filename=ent.name, id=ent.id,
                                              filestatus=current_status['status'].values[0]))

    return({
        'status_list': statuses,
        'error_list': errors,
        'to_validate': to_validate})


def _check_valid(syn, filepaths, center, filetype, filenames,
                 oncotree_link, threads, testing):
    '''
    Function to validate a file
    '''
    # If no filetype set, means the file was named incorrectly
    try:
        valid, message, filetype = validate.validate_single_file(
            syn,
            filepaths,
            center,
            oncotreelink=oncotree_link,
            testing=testing)
        logger.info("VALIDATION COMPLETE")
    except ValueError as e:
        # Specify this as None for the single case where filename
        # validation fails
        filetype = None
        message = str(e)
        logger.error(message)
        valid = False
    return(valid, filetype, message)


def _send_validation_error_email(syn, filenames, message, file_users):
    '''
    Sends validation error email

    Args:
        syn: Synapse object
        filenames: invalid filenames
        message: error message
        file_users: List of unique synapse user profiles of
                    users that created and most recently
                    modified the file
    '''
    # Send email the first time the file is invalid
    incorrect_files = ", ".join(filenames)
    usernames = ", ".join([
        syn.getUserProfile(user)['userName']
        for user in file_users])
    email_message = (
        "Dear {username},\n\n"
        "Your files ({filenames}) are invalid! "
        "Here are the reasons why:\n\n{error_message}".format(
            username=usernames,
            filenames=incorrect_files,
            error_message=message))
    syn.sendMessage(
        file_users, "GENIE Validation Error", email_message)


def _get_status_and_error_list(syn, valid, message, filetype, entities, modified_ons):
    '''
    Helper function to return the status and error list of the
    files based on validation result.

    Args:
        syn: Synapse object
        valid: Boolean value of results of validation
        message: Validation message
        filetype: File type
        entities: List of Synapse Entities
        modified_ons: List of modified on dates

    Returns:
        tuple: input_status_list - status of input files list,
               invalid_errors_list - error list
    '''
    if valid:
        input_status_list = [
            [ent.id, ent.expectedPath, ent.md5, "VALIDATED",
             os.path.basename(ent.expectedPath), modifiedon, filetype]
            for ent, modifiedon in
            zip(entities, modified_ons)]
        invalid_errors_list = None
    else:
        input_status_list = [
            [ent.id, ent.expectedPath, ent.md5, "INVALID",
            os.path.basename(ent.expectedPath), modifiedon, filetype]
            for ent, modifiedon in
            zip(entities, modified_ons)]
        invalid_errors_list = [
            [ent.id, message, os.path.basename(ent.expectedPath)]
            for ent in entities]
    return(input_status_list, invalid_errors_list)


def validatefile(syn, entities, validation_statusdf, error_trackerdf,
                 center, threads, testing, oncotree_link):
    '''Validate a list of entities.

    If a file has not changed, then it doesn't need to be validated.

    Args:
        entitylist: A list of entities for a single file 'type' (usually a single file, but clinical can have two)
        syn: Synapse object
        validation_statusdf: Validation status dataframe
        error_trackerdf: Invalid files error tracking dataframe
        center: Center of interest
        testing: Boolean determining whether using testing parameter
        oncotree_link: Oncotree url

    Returns:
        tuple: input_status_list - status of input files,
               invalid_errors_list - error list

    '''

    filepaths = [entity.path for entity in entities]
    filenames = [os.path.basename(path) for path in filepaths]

    logger.info(
        "VALIDATING {filenames}".format(filenames=", ".join(filenames)))

    modified_ons = [
        synapseclient.utils.to_unix_epoch_time(
            datetime.datetime.strptime(
                entity.modifiedOn.split(".")[0], "%Y-%m-%dT%H:%M:%S"))
        for entity in entities]
    file_users = list(set([entities[0].modifiedBy, entities[0].createdBy]))

    check_file_status = check_existing_file_status(
        validation_statusdf, error_trackerdf, entities)

    status_list = check_file_status['status_list']
    error_list = check_file_status['error_list']
    # Need to figure out to how to remove this

    filetype = validate.determine_filetype(syn, filepaths, center)
    if check_file_status['to_validate']:
        try:
            valid, message, filetype = validate.validate_single_file(
                syn,
                filepaths,
                center,
                oncotreelink=oncotree_link,
                testing=testing)
            logger.info("VALIDATION COMPLETE")
        except ValueError as e:
            # Specify this as None for the single case where filename
            # validation fails
            filetype = None
            message = str(e)
            logger.error(message)
            valid = False
        input_status_list, invalid_errors_list = _get_status_and_error_list(
            syn, valid, message, filetype,
            entities, modified_ons)
        # Send email the first time the file is invalid
        if invalid_errors_list is not None:
            _send_validation_error_email(syn, filenames, message, file_users)
    else:
        input_status_list = [
            [ent.id, path, ent.md5, status, filename, modifiedon, filetype]
            for ent, path, status, filename, modifiedon in
            zip(entities, filepaths, status_list, filenames, modified_ons)]
        invalid_errors_list = [
            [entity.id, error, filename]
            for entity, error, filename in
            zip(entities, error_list, filenames)]
    return(input_status_list, invalid_errors_list)


def processFiles(syn, validFiles, center, path_to_GENIE, threads,
                 center_mapping_df, oncotreeLink, databaseToSynIdMappingDf,
                 validVCF=None, vcf2mafPath=None,
                 veppath=None, vepdata=None,
                 processing="main", test=False, reference=None):
    '''
    Processing validated files

    Args:
        syn: Synapse object
        validFiles: pandas dataframe containing validated files
                    has 'id', 'path', and 'fileType' column
        center: GENIE center name
        path_to_GENIE: Path to GENIE workdir
        threads: Threads used
        center_mapping_df: Center mapping dataframe
        oncotreeLink: Link to oncotree
        databaseToSynIdMappingDf: Database to synapse id mapping dataframe
        validVCF: Valid vcf files
        vcf2mafPath: Path to vcf2maf
        veppath: Path to vep
        vepdata: Path to vep index files
        processing: Processing type. Defaults to main
        test: Test flag
        reference: Reference file for vcf2maf
    '''
    logger.info("PROCESSING {} FILES: {}".format(center, len(validFiles)))
    centerStagingFolder = os.path.join(path_to_GENIE, center)
    centerStagingSynId = center_mapping_df['stagingSynId'][
        center_mapping_df['center'] == center][0]

    if not os.path.exists(centerStagingFolder):
        os.makedirs(centerStagingFolder)
    if processing == "main":
        for fileSynId, filePath, fileType in \
                zip(validFiles['id'],
                    validFiles['path'],
                    validFiles['fileType']):

            filename = os.path.basename(filePath)
            newPath = os.path.join(centerStagingFolder, filename)
            # store = True
            synId = databaseToSynIdMappingDf.Id[
                databaseToSynIdMappingDf['Database'] == fileType]
            if len(synId) == 0:
                synId = None
            else:
                synId = synId[0]
            # if fileType not in [None,"cna"]:
            if fileType is not None:
                processor = PROCESS_FILES[fileType](syn, center, threads)
                processor.process(
                    filePath=filePath, newPath=newPath,
                    parentId=centerStagingSynId, databaseSynId=synId,
                    oncotreeLink=oncotreeLink,
                    fileSynId=fileSynId, validVCF=validVCF,
                    path_to_GENIE=path_to_GENIE, vcf2mafPath=vcf2mafPath,
                    veppath=veppath, vepdata=vepdata,
                    processing=processing,
                    databaseToSynIdMappingDf=databaseToSynIdMappingDf,
                    reference=reference, test=test)

    elif processing in ["vcf", "maf", "mafSP"]:
        filePath = None
        newPath = None
        fileType = None
        synId = databaseToSynIdMappingDf.Id[
            databaseToSynIdMappingDf['Database'] == processing][0]
        fileSynId = None
        processor = PROCESS_FILES[processing](syn, center, threads)
        processor.process(
            filePath=filePath, newPath=newPath,
            parentId=centerStagingSynId, databaseSynId=synId,
            oncotreeLink=oncotreeLink,
            fileSynId=fileSynId, validVCF=validVCF,
            path_to_GENIE=path_to_GENIE, vcf2mafPath=vcf2mafPath,
            veppath=veppath, vepdata=vepdata,
            processing=processing,
            databaseToSynIdMappingDf=databaseToSynIdMappingDf,
            reference=reference)

    logger.info("ALL DATA STORED IN DATABASE")

# def _create_maf_db(syn, foo):
#     maf_database_ent = syn.get(maf_database_synid)
#     print(maf_database_ent)
#     maf_columns = list(syn.getTableColumns(maf_database_synid))
#     schema = synapseclient.Schema(
#         name='Narrow MAF {current_time} Database'.format(
#             current_time=time.time()),
#         columns=maf_columns,
#         parent=process_functions.getDatabaseSynId(
#             syn, "main",
#             databaseToSynIdMappingDf=database_synid_mappingdf))
#     schema.primaryKey = maf_database_ent.primaryKey
#     new_maf_database = syn.store(schema)

# TODO: Should split this into 3 funcitons
# so that unit tests are easier to write


def create_and_archive_maf_database(syn, database_synid_mappingdf):
    '''
    Creates new MAF database and archives the old database in the staging site

    Args:
        syn: Synapse object
        databaseToSynIdMappingDf: Database to synapse id mapping dataframe

    Return:
        Editted database to synapse id mapping dataframe
    '''
    maf_database_synid = process_functions.getDatabaseSynId(
        syn, "vcf2maf", databaseToSynIdMappingDf=database_synid_mappingdf)
    maf_database_ent = syn.get(maf_database_synid)
    maf_columns = list(syn.getTableColumns(maf_database_synid))
    schema = synapseclient.Schema(
        name='Narrow MAF {current_time} Database'.format(
            current_time=time.time()),
        columns=maf_columns,
        parent=process_functions.getDatabaseSynId(
            syn, "main", databaseToSynIdMappingDf=database_synid_mappingdf))
    schema.primaryKey = maf_database_ent.primaryKey
    new_maf_database = syn.store(schema)
    # Store in the new database synid
    database_synid_mappingdf['Id'][
        database_synid_mappingdf[
            'Database'] == 'vcf2maf'] = new_maf_database.id

    vcf2maf_mappingdf = database_synid_mappingdf[
        database_synid_mappingdf['Database'] == 'vcf2maf']
    # vcf2maf_mappingdf['Id'][0] = newMafDb.id
    # Update this synid later
    syn.store(synapseclient.Table("syn12094210", vcf2maf_mappingdf))
    # Move and archive old mafdatabase (This is the staging synid)
    maf_database_ent.parentId = "syn7208886"
    maf_database_ent.name = "ARCHIVED " + maf_database_ent.name
    syn.store(maf_database_ent)
    # maf_database_synid = new_maf_database.id
    # Remove can download permissions from project GENIE team
    syn.setPermissions(new_maf_database.id, 3326313, [])
    return(database_synid_mappingdf)


def email_duplication_error(syn, duplicated_filesdf):
    '''
    Sends an email if there is a duplication error

    Args:
        syn: Synapse object
        duplicated_filesdf: dataframe with 'id', 'name' column
    '''
    if not duplicated_filesdf.empty:
        incorrect_files = [
            name for synId, name in zip(duplicated_filesdf['id'],
                                        duplicated_filesdf['name'])]
        incorrect_filenames = ", ".join(incorrect_files)
        incorrect_ent = syn.get(duplicated_filesdf['id'].iloc[0])
        send_to_users = set([incorrect_ent.modifiedBy,
                             incorrect_ent.createdBy])
        usernames = ", ".join(
            [syn.getUserProfile(user)['userName'] for user in send_to_users])
        error_email = (
            "Dear %s,\n\n"
            "Your files (%s) are duplicated!  FILES SHOULD BE UPLOADED AS "
            "NEW VERSIONS AND THE ENTIRE DATASET SHOULD BE "
            "UPLOADED EVERYTIME".format(usernames, incorrect_filenames))
        syn.sendMessage(
            list(send_to_users), "GENIE Validation Error", error_email)


def get_duplicated_files(syn, validation_statusdf, duplicated_error_message):
    '''
    Check for duplicated files.  There should be no duplication,
    files should be uploaded as new versions and the entire dataset
    should be uploaded everytime

    Args:
        syn: Synapse object
        validation_statusdf: dataframe with 'name' and 'id' column
        duplicated_error_message: Error message for duplicated files

    Returns:
        dataframe with 'id', 'name' and 'errors' of duplicated files
    '''
    logger.info("CHECK FOR DUPLICATED FILES")
    duplicated_filesdf = validation_statusdf[
        validation_statusdf['name'].duplicated(keep=False)]
    # cbs/seg files should not be duplicated.
    cbs_seg_files = validation_statusdf.query(
        'name.str.endswith("cbs") or name.str.endswith("seg")')
    if len(cbs_seg_files) > 1:
        duplicated_filesdf = duplicated_filesdf.append(cbs_seg_files)
    # clinical files should not be duplicated.
    clinical_files = validation_statusdf.query(
        'name.str.startswith("data_clinical_supp")')
    if len(clinical_files) > 2:
        duplicated_filesdf = duplicated_filesdf.append(clinical_files)
    duplicated_filesdf.drop_duplicates("id", inplace=True)
    logger.info("THERE ARE {} DUPLICATED FILES".format(
        len(duplicated_filesdf)))
    duplicated_filesdf['errors'] = duplicated_error_message
    return(duplicated_filesdf)


def validation(syn, center, process,
               center_mapping_df, databaseToSynIdMappingDf,
               thread, testing, oncotreeLink):
    '''
    Validation of all center files

    Args:
        syn: Synapse object
        center: Center name
        process: main, vcf, maf
        center_mapping_df: center mapping dataframe
        thread: Unused parameter for now
        testing: True if testing
        oncotreeLink: Link to oncotree

    Returns:
        dataframe: Valid files
    '''
    centerInputSynId = center_mapping_df['inputSynId'][
        center_mapping_df['center'] == center][0]
    logger.info("Center: " + center)
    allFiles = get_center_input_files(syn, centerInputSynId, center, process)

    # If a center has no files, then return empty list
    if not allFiles:
        logger.info("%s has not uploaded any files" % center)
        return([])
    else:
        # Make sure the vcf validation statuses don't get wiped away
        if process != "vcf":
            addToQuery = "and name not like '%.vcf'"
        else:
            addToQuery = ''
        validationStatus = syn.tableQuery(
            "SELECT id,md5,status,name,center,modifiedOn FROM {} where center = '{}' {}".format(
                process_functions.getDatabaseSynId(
                    syn, "validationStatus",
                    databaseToSynIdMappingDf=databaseToSynIdMappingDf),
                center,
                addToQuery))
        validation_statusdf = validationStatus.asDataFrame()
        errorTracker = syn.tableQuery(
            "SELECT id,center,errors,name FROM {} where center = '{}' {}".format(
                process_functions.getDatabaseSynId(
                    syn, "errorTracker",
                    databaseToSynIdMappingDf=databaseToSynIdMappingDf),
                center,
                addToQuery))
        error_trackerdf = errorTracker.asDataFrame()

        inputValidStatus = []
        invalidErrors = []
        
        for ents in allFiles:
            status, errors = input_to_database.validatefile(syn, ents,
                                                            validation_statusdf, 
                                                            error_trackerdf, 
                                                            center='SAGE', threads=1, 
                                                            testing=False, 
                                                            oncotree_link=None)
            inputValidStatus.extend(status)
            if errors is not None:
                invalidErrors.extend(errors)

        inputValidStatus = pd.DataFrame(
            inputValidStatus,
            columns=[
                "id", 'path', 'md5', 'status',
                'name', 'modifiedOn', 'fileType'])

        duplicatedFileError = (
            "DUPLICATED FILENAME! FILES SHOULD BE UPLOADED AS NEW VERSIONS "
            "AND THE ENTIRE DATASET SHOULD BE UPLOADED EVERYTIME")
        duplicatedFiles = get_duplicated_files(
            syn, inputValidStatus, duplicatedFileError)
        # Send an email if there are any duplicated files
        if not duplicatedFiles.empty:
            email_duplication_error(syn, duplicatedFiles)

        inputValidStatus['status'][
            inputValidStatus['id'].isin(duplicatedFiles['id'])] = "INVALID"
        # Create invalid error synapse table
        logger.info("UPDATE INVALID FILE REASON DATABASE")
        invalidErrors = pd.DataFrame(
            invalidErrors, columns=["id", 'errors', 'name'])
        # Remove fixed duplicated files
        # This makes sure that the files removed actually had duplicated file
        # errors and not some other error
        dupIds = invalidErrors['id'][
            invalidErrors['errors'] == duplicatedFileError]
        removeIds = dupIds[~dupIds.isin(duplicatedFiles['id'])]
        invalidErrors = invalidErrors[~invalidErrors['id'].isin(removeIds)]
        # Append duplicated file errors
        invalidErrors = invalidErrors.append(
            duplicatedFiles[['id', 'errors', 'name']])
        invalidErrors['center'] = center
        invalidIds = inputValidStatus['id'][
            inputValidStatus['status'] == "INVALID"]
        invalidErrors = invalidErrors[invalidErrors['id'].isin(invalidIds)]
        process_functions.updateDatabase(
            syn, errorTracker.asDataFrame(), invalidErrors,
            process_functions.getDatabaseSynId(
                syn, "errorTracker",
                databaseToSynIdMappingDf=databaseToSynIdMappingDf),
            ["id"], to_delete=True)

        paths = inputValidStatus['path']
        # filenames = [os.path.basename(name) for name in paths]
        del inputValidStatus['path']
        logger.info("UPDATE VALIDATION STATUS DATABASE")
        inputValidStatus['center'] = center
        # Remove fixed duplicated files
        inputValidStatus = inputValidStatus[
            ~inputValidStatus['id'].isin(removeIds)]

        process_functions.updateDatabase(
            syn,
            validationStatus.asDataFrame(),
            inputValidStatus[
                ["id", 'md5', 'status', 'name', 'center', 'modifiedOn']],
            process_functions.getDatabaseSynId(
                syn, "validationStatus",
                databaseToSynIdMappingDf=databaseToSynIdMappingDf),
            ["id"],
            to_delete=True)

        inputValidStatus['path'] = paths
        validFiles = inputValidStatus[['id', 'path', 'fileType']][
            inputValidStatus['status'] == "VALIDATED"]
        return(validFiles)


def center_input_to_database(
        syn, center, process, testing,
        only_validate, vcf2maf_path, vep_path,
        vep_data, database_to_synid_mappingdf,
        center_mapping_df, reference=None,
        delete_old=False, oncotree_link=None, thread=1):
    if only_validate:
        log_path = os.path.join(
            process_functions.SCRIPT_DIR,
            "{}_validation_log.txt".format(center))
    else:
        log_path = os.path.join(
            process_functions.SCRIPT_DIR,
            "{}_{}_log.txt".format(center, process))

    logFormatter = logging.Formatter(
        "%(asctime)s [%(name)s][%(levelname)s] %(message)s")
    fileHandler = logging.FileHandler(log_path, mode='w')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    if testing:
        logger.info("###########################################")
        logger.info("############NOW IN TESTING MODE############")
        logger.info("###########################################")

    # ----------------------------------------
    # Start input to staging process
    # ----------------------------------------
    '''
    # path_to_genie = os.path.realpath(os.path.join(
    #    process_functions.SCRIPT_DIR, "../"))
    Make the synapsecache dir the genie input folder for now
    The main reason for this is because the .synaspecache dir
    is mounted by batch
    '''
    path_to_genie = os.path.expanduser("~/.synapseCache")
    # Create input and staging folders
    if not os.path.exists(os.path.join(path_to_genie, center, "input")):
        os.makedirs(os.path.join(path_to_genie, center, "input"))
    if not os.path.exists(os.path.join(path_to_genie, center, "staging")):
        os.makedirs(os.path.join(path_to_genie, center, "staging"))

    if delete_old:
        process_functions.rmFiles(os.path.join(path_to_genie, center))

    validFiles = validation(
        syn, center, process, center_mapping_df,
        database_to_synid_mappingdf, thread,
        testing, oncotree_link)

    if len(validFiles) > 0 and not only_validate:
        # Reorganize so BED file are always validated and processed first
        validBED = [
            os.path.basename(i).endswith('.bed') for i in validFiles['path']]
        beds = validFiles[validBED]
        validFiles = beds.append(validFiles)
        validFiles.drop_duplicates(inplace=True)
        # Valid vcf files
        validVCF = [
            i for i in validFiles['path']
            if os.path.basename(i).endswith('.vcf')]

        processTrackerSynId = process_functions.getDatabaseSynId(
            syn, "processTracker",
            databaseToSynIdMappingDf=database_to_synid_mappingdf)
        # Add process tracker for time start
        processTracker = syn.tableQuery(
            "SELECT timeStartProcessing FROM {} "
            "where center = '{}' and "
            "processingType = '{}'".format(
                processTrackerSynId, center, process))
        processTrackerDf = processTracker.asDataFrame()
        if len(processTrackerDf) == 0:
            new_rows = [[
                center,
                str(int(time.time()*1000)),
                str(int(time.time()*1000)),
                process]]

            syn.store(synapseclient.Table(
                processTrackerSynId, new_rows))
        else:
            processTrackerDf['timeStartProcessing'][0] = \
                str(int(time.time()*1000))
            syn.store(synapseclient.Table(
                processTrackerSynId, processTrackerDf))

        processFiles(
            syn, validFiles, center, path_to_genie, thread,
            center_mapping_df, oncotree_link, database_to_synid_mappingdf,
            validVCF=validVCF,
            vcf2mafPath=vcf2maf_path,
            veppath=vep_path, vepdata=vep_data,
            test=testing, processing=process, reference=reference)

        # Should add in this process end tracking
        # before the deletion of samples
        processTracker = syn.tableQuery(
            "SELECT timeEndProcessing FROM {synid} where center = '{center}' "
            "and processingType = '{processtype}'".format(
                synid=processTrackerSynId,
                center=center,
                processtype=process))
        processTrackerDf = processTracker.asDataFrame()
        processTrackerDf['timeEndProcessing'][0] = str(int(time.time()*1000))
        syn.store(synapseclient.Table(processTrackerSynId, processTrackerDf))

        logger.info("SAMPLE/PATIENT RETRACTION")
        toRetract.retract(syn, testing)

    else:
        messageOut = \
            "{} does not have any valid files" if not only_validate \
            else "ONLY VALIDATION OCCURED FOR {}"
        logger.info(messageOut.format(center))

    # Store log file
    log_folder_synid = process_functions.getDatabaseSynId(
        syn, "logs", databaseToSynIdMappingDf=database_to_synid_mappingdf)
    syn.store(synapseclient.File(log_path, parentId=log_folder_synid))
    os.remove(log_path)
    logger.info("ALL PROCESSES COMPLETE")
