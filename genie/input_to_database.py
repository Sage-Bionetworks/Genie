#!/usr/bin/env python3
import logging
logger = logging.getLogger("genie")
import synapseclient
import synapseutils
import argparse
import os
import pandas as pd
import datetime
import shutil
import time
# Configuration file
from genie import PROCESS_FILES, process_functions, validate, \
                  toRetract, write_invalid_reasons

'''
TODO:
Could potentially get all the inforamation of the file entity right here
To avoid the syn.get rest call later which doesn't actually download the file
'''


def rename_file(syn, synid):
    '''
    Gets file from synapse and renames the file if necessary.

    Args:
        syn: Synapse object
        synid : Synapse id or entity

    Returns:
        Path of corrected file
    '''
    ent = syn.get(synid)
    dirpath = os.path.dirname(ent.path)
    expectedpath = os.path.join(dirpath, ent.name)
    if expectedpath != ent.path:
        shutil.copyfile(ent.path, expectedpath)
    return(expectedpath)


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
        List of Tuples with the correct format to pass into validation
    '''
    logger.info("GETTING {center} INPUT FILES".format(center=center))
    clinical_pair_name = [
        "data_clinical_supp_sample_{center}.txt".format(center=center),
        "data_clinical_supp_patient_{center}.txt".format(center=center)]

    center_files = synapseutils.walk(syn, synid)
    clinicalpair = []
    prepared_center_file_list = []
    for dirpath, dirname, filenames in center_files:
        for name, ent_synid in filenames:
            logger.info(name)
            paired = False
            '''
            Clinical file can come as two files.
            The two files need to be merged together which is
            why there is this format
            '''
            if name in clinical_pair_name:
                paired = True
                clinicalpair.append(ent_synid)
            if len(clinicalpair) == 2:
                syns = [i for i in clinicalpair]
                paths = [rename_file(syn, i) for i in clinicalpair]
                clinicalpair = []
                prepared_center_file_list.append((syns, paths))
            elif not paired:
                '''
                This is to remove vcfs from being validated
                during main processing.  Often there are too many vcf files.
                Not necessary for them to be run everytime.
                '''
                if process == "vcf":
                    prepared_center_file_list.append(
                        ([ent_synid], [rename_file(syn, ent_synid)]))
                elif not name.endswith(".vcf"):
                    prepared_center_file_list.append(
                        ([ent_synid], [rename_file(syn, ent_synid)]))
    return(prepared_center_file_list)


def get_filetype(syn, path_list, center):
    '''
    Get the file type of the file by validating its filename

    Args:
        syn: Synapse object
        path_list: list of filepaths to center files
        center: Participating Center

    Returns:
        str: File type of input files
    '''
    filetype = None
    for file_format in PROCESS_FILES:
        try:
            filetype = PROCESS_FILES[file_format](
                syn, center).validateFilename(path_list)
        except AssertionError:
            continue
        # If valid filename, return file type.
        if filetype is not None:
            break
    return(filetype)


def check_existing_file_status(validation_statusdf, error_trackerdf,
                               entities, input_filenames):
    '''
    This function checks input files against the existing validation and error
    tracking dataframe

    Args:
        validation_statusdf: Validation status dataframe
        error_trackerdf: Error tracking dataframe
        entities: list of center input entites
        input_filenames: List of center input filenames

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
    for ent, input_filename in zip(entities, input_filenames):
        input_validation_status = \
            validation_statusdf[validation_statusdf['id'] == ent.id]
        input_error = error_trackerdf[error_trackerdf['id'] == ent.id]
        if input_validation_status.empty:
            to_validate = True
        else:
            # This to_validate is here, because the following is a
            # sequential check of whether files need to be validated
            to_validate = False
            statuses.append(input_validation_status['status'].values[0])
            if input_error.empty:
                to_validate = \
                    input_validation_status['status'].values[0] == "INVALID"
            else:
                errors.append(input_error['errors'].values[0])
            # Add Name check here (must add name of the entity as a column)
            if input_validation_status['md5'].values[0] != ent.md5 or \
               input_validation_status['name'].values[0] != input_filename:
                to_validate = True
            else:
                logger.info(
                    "{filename} FILE STATUS IS: {filestatus}".format(
                      filename=input_filename,
                      filestatus=input_validation_status['status'].values[0]))

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
    if filetype is None:
        message = (
            "{filenames}: Incorrect filenaming convention or can't be "
            "processed".format(filenames=", ".join(filenames)))
        logger.error(message)
        valid = False
    else:
        try:
            message, valid = validate.validate(
                syn,
                filetype,
                filepaths,
                center,
                threads,
                oncotree_url=oncotree_link,
                testing=testing)
            logger.info("VALIDATION COMPLETE")
        except ValueError as e:
            logger.error(e)
            message = e
            valid = False
    return(valid, message)


def _get_status_and_error_list(syn, fileinfo, valid, message, filetype,
                               entities, filepaths, filenames, modified_ons):
    '''
    '''
    if valid:
        input_status_list = [
            [ent.id, filepath, ent.md5, "VALIDATED",
             filename, modifiedon, filetype]
            for ent, filepath, filename, modifiedon in
            zip(entities, filepaths, filenames, modified_ons)]

        invalid_errors_list = None
    else:
        # Send email the first time the file is invalid
        incorrect_files = ", ".join(filenames)
        incorrect_ent = syn.get(fileinfo['synId'][0])
        find_file_users = \
            list(set([incorrect_ent.modifiedBy, incorrect_ent.createdBy]))
        usernames = ", ".join([
            syn.getUserProfile(user)['userName']
            for user in find_file_users])
        email_message = (
            "Dear {username},\n\n"
            "Your files ({filenames}) are invalid! "
            "Here are the reasons why:\n\n{error_message}".format(
                username=usernames,
                filenames=incorrect_files,
                error_message=message))
        syn.sendMessage(
            find_file_users, "GENIE Validation Error", email_message)
        input_status_list = [
            [ent.id, path, ent.md5, "INVALID", name, modifiedon, filetype]
            for ent, path, name, modifiedon in
            zip(entities, filepaths, filenames, modified_ons)]

        invalid_errors_list = [
            [synid, message, filename]
            for synid, filename in zip(fileinfo['synId'], filenames)]
    return(input_status_list, invalid_errors_list)


def validatefile(fileinfo,
                 syn,
                 validation_statusdf,
                 error_trackerdf,
                 center,
                 threads,
                 testing,
                 oncotree_link):
    '''
    Function that is applied to a pandas dataframe to
    validates each row. If a file has not changed, then it
    doesn't need to be validated

    Args:
        fileinfo: A row passed in as a Series through the apply function
                  in pandas
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

    filenames = [os.path.basename(i) for i in fileinfo['filePaths']]
    logger.info(
        "VALIDATING {filenames}".format(filenames=", ".join(filenames)))
    filepaths = fileinfo['filePaths']
    entities = [
        syn.get(synid, downloadFile=False) for synid in fileinfo['synId']]
    modified_ons = [
        synapseclient.utils.to_unix_epoch_time(
            datetime.datetime.strptime(
                entity.modifiedOn.split(".")[0], "%Y-%m-%dT%H:%M:%S"))
        for entity in entities]

    check_file_status = check_existing_file_status(
        validation_statusdf, error_trackerdf, entities, filenames)

    status_list = check_file_status['status_list']
    error_list = check_file_status['error_list']
    filetype = get_filetype(syn, filepaths, center)
    if check_file_status['to_validate']:

        valid, message = _check_valid(
            syn, filepaths, center, filetype, filenames,
            oncotree_link, threads, testing)

        input_status_list, invalid_errors_list = _get_status_and_error_list(
            syn, fileinfo, valid, message, filetype,
            entities, filepaths, filenames, modified_ons)

    else:
        input_status_list = [
            [ent.id, path, ent.md5, status, filename, modifiedon, filetype]
            for ent, path, status, filename, modifiedon in
            zip(entities, filepaths, status_list, filenames, modified_ons)]
        invalid_errors_list = [
            [synid, error, filename]
            for synid, error, filename in
            zip(fileinfo['synId'], error_list, filenames)]
    return(input_status_list, invalid_errors_list)


def processFiles(syn, validFiles, center, path_to_GENIE, threads,
                 center_mapping_df, oncotreeLink, databaseToSynIdMappingDf,
                 validVCF=None, vcf2mafPath=None,
                 veppath=None, vepdata=None,
                 processing="main", test=False, reference=None):
    '''
    Processing a single file
    '''
    logger.info("PROCESSING %s FILES: %d" % (center, len(validFiles)))
    centerStagingFolder = os.path.join(path_to_GENIE, center)
    centerStagingSynId = center_mapping_df['stagingSynId'][
        center_mapping_df['center'] == center][0]
    # PROCESS_FILES is in config_process_scripts.py
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
                PROCESS_FILES[fileType](syn, center, threads).process(
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
        PROCESS_FILES[processing](syn, center, threads).process(
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

    allFiles = pd.DataFrame(allFiles, columns=['synId', 'filePaths'])
    # If a center has no files, then return empty list
    if allFiles.empty:
        logger.info("%s has not uploaded any files" % center)
        return([])
    else:
        # Make sure the vcf validation statuses don't get wiped away
        if process != "vcf":
            addToQuery = "and name not like '%.vcf'"
        else:
            addToQuery = ''
        validationStatus = syn.tableQuery(
            "SELECT * FROM {} where center = '{}' {}".format(
                process_functions.getDatabaseSynId(
                    syn, "validationStatus",
                    databaseToSynIdMappingDf=databaseToSynIdMappingDf),
                center,
                addToQuery))

        errorTracker = syn.tableQuery(
            "SELECT * FROM {} where center = '{}' {}".format(
                process_functions.getDatabaseSynId(
                    syn, "errorTracker",
                    databaseToSynIdMappingDf=databaseToSynIdMappingDf),
                center,
                addToQuery))

        # VALIDATE FILES
        validationStatusDf = validationStatus.asDataFrame()
        errorTrackerDf = errorTracker.asDataFrame()
        validated = allFiles.apply(
            lambda fileinfo: validatefile(
                fileinfo, syn, validationStatusDf,
                errorTrackerDf, center, thread,
                testing, oncotreeLink), axis=1)

        inputValidStatus = []
        invalidErrors = []
        for inputStat, invalErrors in validated:
            inputValidStatus.extend(inputStat)
            if invalErrors is not None:
                invalidErrors.extend(invalErrors)
        inputValidStatus = pd.DataFrame(
            inputValidStatus,
            columns=[
                "id", 'path', 'md5', 'status',
                'name', 'modifiedOn', 'fileType'])

        logger.info("CHECK FOR DUPLICATED FILES")
        '''
        Check for duplicated filenames.
        There should be no duplication, files should be uploaded as
        new versions and the entire dataset should be uploaded everytime
        cbs and seg files should not be duplicated.  There can only be one
        '''
        duplicatedFiles = inputValidStatus[
            inputValidStatus['name'].duplicated(keep=False)]
        cbsSegBool = [
            os.path.basename(i).endswith('.cbs') or
            os.path.basename(i).endswith('.seg')
            for i in inputValidStatus['name']]
        cbsSegFiles = inputValidStatus[cbsSegBool]
        if len(cbsSegFiles) > 1:
            duplicatedFiles = duplicatedFiles.append(cbsSegFiles)

        duplicatedFiles.drop_duplicates("id", inplace=True)
        inputValidStatus['status'][
            inputValidStatus['id'].isin(duplicatedFiles['id'])] = "INVALID"
        duplicatedFileError = (
            "DUPLICATED FILENAME! FILES SHOULD BE UPLOADED AS NEW VERSIONS "
            "AND THE ENTIRE DATASET SHOULD BE UPLOADED EVERYTIME")
        duplicatedFiles['errors'] = duplicatedFileError
        # Send an email if there are any duplicated files
        if not duplicatedFiles.empty:
            incorrectFiles = ", ".join(
                [name for synId, name in
                 zip(duplicatedFiles['id'], duplicatedFiles['name'])])
            incorrectEnt = syn.get(duplicatedFiles['id'].iloc[0])
            sendEmail = set([incorrectEnt.modifiedBy, incorrectEnt.createdBy])
            userNames = ", ".join(
                [syn.getUserProfile(user).userName for user in sendEmail])
            errorEmail = (
                "Dear %s,\n\n"
                "Your files (%s) are duplicated!  FILES SHOULD BE UPLOADED AS "
                "NEW VERSIONS AND THE ENTIRE DATASET SHOULD BE "
                "UPLOADED EVERYTIME".format(userNames, incorrectFiles))
            syn.sendMessage(
                list(sendEmail), "GENIE Validation Error", errorEmail)

        logger.info("THERE ARE %d DUPLICATED FILES" % len(duplicatedFiles))

        # Create invalid error synapse table
        logger.info("UPDATE INVALID FILE REASON DATABASE")
        invalidErrors = pd.DataFrame(
            invalidErrors, columns=["id", 'errors', 'name'])
        # Remove fixed duplicated files
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

        # Resolve with https://github.com/Sage-Bionetworks/Genie/issues/94
        # logger.info("SAMPLE/PATIENT RETRACTION")
        # toRetract.retract(syn, testing)

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
