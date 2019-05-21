#! /usr/bin/env python
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import synapseclient
from synapseclient import File, Table
import synapseutils
import argparse
import os
from multiprocessing import Pool
import pandas as pd
import subprocess
import datetime
import subprocess
import shutil
import time
import toRetract
import write_invalid_reasons

#Configuration file
from genie import PROCESS_FILES, process_functions, validate

def reNameFile(syn, synId):
    temp = syn.get(synId)
    dirname = os.path.dirname(temp.path)
    newPath = os.path.join(dirname, temp.name)
    if newPath != temp.path:
        shutil.copyfile(temp.path, newPath)
    return(newPath)

def getCenterInputFiles(syn, synId, center, process="main"):
    """
    This function walks through each center's input directory and validates every file
    """
    ################################################################
    ##If a file has not changed than it does not need to be processed!!!!
    ################################################################
    logger.info("GETTING %s INPUT FILES" % center)
    CLINICAL_PAIR_NAME = ["data_clinical_supp_sample_%s.txt" % center, "data_clinical_supp_patient_%s.txt" % center]
    walked = synapseutils.walk(syn, synId)
    clinicalpair = []
    allFiles = []
    for dirpath, dirname, filenames in walked:
        for name, synid in filenames:
            logger.info(name)
            paired=False
            if name in CLINICAL_PAIR_NAME:
                paired = True
                clinicalpair.append(synid)
            if len(clinicalpair) == 2:
                syns = [i for i in clinicalpair]
                paths = [reNameFile(syn, i) for i in clinicalpair]
                clinicalpair = []
                allFiles.append((syns,paths))
            elif not paired:
                if process == "vcf":
                    allFiles.append(([synid],[reNameFile(syn, synid)]))
                elif not name.endswith(".vcf"):
                    allFiles.append(([synid],[reNameFile(syn, synid)]))
    return(allFiles)

def getFileType(syn, paths, name, center):
    fileType = None
    for fileNameFormat in PROCESS_FILES:
        try:
            fileType = PROCESS_FILES[fileNameFormat](syn, center).validateFilename(paths)
        except AssertionError as e:
            continue
        if fileType is not None:
            break
    return(fileType)

def validateFile(syn, validationStatusDf, errorTracker, center, threads, x, testing, oncotreeLink):
    names = [os.path.basename(i) for i in x['filePaths']]
    logger.info("VALIDATING %s" % ", ".join(names))

    paths = x['filePaths']
    name = names[0]
    entities = [syn.get(synId,downloadFile=False) for synId in x['synId']]
    md5s = [entity.md5 for entity in entities]
    modifiedOns = [synapseclient.utils.to_unix_epoch_time(datetime.datetime.strptime(entity.modifiedOn.split(".")[0], "%Y-%m-%dT%H:%M:%S")) for entity in entities]

    toValidate = False
    statuses = []
    errors = []
    # Check validation status and md5 of file
    for synId,md5,filename in zip(x['synId'],md5s,names):
        checkValid = validationStatusDf[validationStatusDf['id'] == synId]
        checkError = errorTracker[errorTracker['id'] == synId]
        if checkValid.empty:
            toValidate = True
        else:
            statuses.append(checkValid['status'].values[0])
            if checkError.empty:
                if checkValid['status'].values[0] == "INVALID":
                    toValidate =True
            else:
                errors.append(checkError['errors'].values[0])
            #Add Name check here (must add name of the entity as a column)
            if checkValid['md5'].values[0] != md5 or checkValid['name'].values[0] != filename:
                toValidate = True
            else:
                logger.info("%s FILE STATUS IS: %s" % (filename, checkValid['status'].values[0]))
    fileType = getFileType(syn, paths, name, center)
    if toValidate:
        # If no filetype set, means the file was named incorrectly
        if fileType is None:
            message = "%s: Incorrect filenaming convention or can't be processed" % name
            logger.error(message)
            valid=False
        else:
            try:
                message, valid = validate.validate(syn, fileType, paths, center, threads, oncotree_url=oncotreeLink, testing=testing)
                logger.info("VALIDATION COMPLETE")
            except ValueError as e:
                logger.error(e)
                message = e
                valid = False
        if valid:
            return([[synId,path,md5,"VALIDATED",name, modifiedOn,fileType] for synId, path, md5, name, modifiedOn in zip(x['synId'],paths, md5s, names, modifiedOns)],None)
        else:
            #Send email the first time the file is invalid
            incorrectFiles = ", ".join([name for synId, name in zip(x['synId'],names)])
            incorrectEnt = syn.get(x['synId'][0])
            sendEmail = set([incorrectEnt.modifiedBy, incorrectEnt.createdBy])
            userNames = ", ".join([syn.getUserProfile(user).userName for user in sendEmail])
            syn.sendMessage(list(sendEmail), "GENIE Validation Error", "Dear %s,\n\nYour files (%s) are invalid! Here are the reasons why:\n\n%s" % (userNames, incorrectFiles, message))
            return([[synId,path,md5,"INVALID",name,modifiedOn,fileType] for synId, path, md5, name, modifiedOn in zip(x['synId'],paths, md5s, names, modifiedOns)],[[synId, message, name] for synId, name in zip(x['synId'],names)])
    else:
        return([[synId,path,md5,status,name,modifiedOn,fileType] for synId, path, md5, status, name, modifiedOn in zip(x['synId'],paths, md5s, statuses, names, modifiedOns)], [[synId, errorMes, name] for synId, errorMes, name in zip(x['synId'],errors,names)])

########################################################################
#Processing files
########################################################################
#Processing single file
def processFiles(syn, validFiles, center, path_to_GENIE, threads, 
                 center_mapping_df, oncotreeLink, databaseToSynIdMappingDf, 
                 validVCF=None, vcf2mafPath=None,
                 veppath=None, vepdata=None,
                 processing="main", test=False, reference=None):

    logger.info("PROCESSING %s FILES: %d" % (center, len(validFiles)))
    centerStagingFolder = os.path.join(path_to_GENIE, center)
    centerStagingSynId = center_mapping_df['stagingSynId'][center_mapping_df['center'] == center][0]
    #PROCESS_FILES is in config_process_scripts.py
    if not os.path.exists(centerStagingFolder):
        os.makedirs(centerStagingFolder)
    if processing == "main":
        for fileSynId, filePath, fileType in zip(validFiles['id'],validFiles['path'],validFiles['fileType']):
            filename = os.path.basename(filePath)
            newPath = os.path.join(centerStagingFolder, filename)
            store = True
            synId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == fileType]
            if len(synId) == 0:
                synId = None
            else:
                synId = synId[0]
            if fileType is not None:
            #if fileType not in [None,"cna"]:
                PROCESS_FILES[fileType](syn, center, threads).process(filePath=filePath, newPath=newPath, 
                                    parentId=centerStagingSynId, databaseSynId=synId, oncotreeLink=oncotreeLink, 
                                    fileSynId=fileSynId, validVCF=validVCF, 
                                    path_to_GENIE=path_to_GENIE, vcf2mafPath=vcf2mafPath,
                                    veppath=veppath,vepdata=vepdata,
                                    processing=processing,databaseToSynIdMappingDf=databaseToSynIdMappingDf, reference=reference, test=test)

    elif processing in ["vcf","maf","mafSP"]:
        filePath = None
        newPath = None
        fileType = None
        synId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == processing][0]
        fileSynId = None
        PROCESS_FILES[processing](syn, center, threads).process(filePath=filePath, newPath=newPath, 
                                    parentId=centerStagingSynId, databaseSynId=synId, oncotreeLink=oncotreeLink, 
                                    fileSynId=fileSynId, validVCF=validVCF,
                                    path_to_GENIE=path_to_GENIE, vcf2mafPath=vcf2mafPath,
                                    veppath=veppath,vepdata=vepdata,
                                    processing=processing,databaseToSynIdMappingDf=databaseToSynIdMappingDf, reference=reference)

    logger.info("ALL DATA STORED IN DATABASE")

#Create and archive maf database
def createMafDatabase(syn, databaseToSynIdMappingDf,testing=False,staging=False):
    mafDatabaseSynId = process_functions.getDatabaseSynId(syn, "vcf2maf", databaseToSynIdMappingDf=databaseToSynIdMappingDf)
    mafDatabaseEnt = syn.get(mafDatabaseSynId)
    mafCols = list(syn.getTableColumns(mafDatabaseSynId))
    schema = synapseclient.Schema(name='Narrow MAF %s Database' % time.time(), columns=mafCols, parent=process_functions.getDatabaseSynId(syn, "main", databaseToSynIdMappingDf=databaseToSynIdMappingDf))
    schema.primaryKey = mafDatabaseEnt.primaryKey
    newMafDb = syn.store(schema)
    #Store in the new database synid
    databaseToSynIdMappingDf['Id'][0] = newMafDb.id
    syn.store(synapseclient.Table(process_functions.getDatabaseSynId(syn, "dbMapping", test=testing),databaseToSynIdMappingDf))
    if not staging and not testing:
        #Make sure to store the newly created maf db synid into the staging synapse mapping
        databaseToSynIdMapping = syn.tableQuery("SELECT * FROM syn12094210 where Database = 'vcf2maf'")
        databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()
        databaseToSynIdMappingDf['Id'][0] = newMafDb.id
        syn.store(synapseclient.Table("syn12094210",databaseToSynIdMappingDf))
    #Move and archive old mafdatabase
    mafDatabaseEnt.parentId = "syn7208886"
    mafDatabaseEnt.name = "ARCHIVED " + mafDatabaseEnt.name
    syn.store(mafDatabaseEnt)
    mafDatabaseSynId = newMafDb.id
    #Remove can download permissions from project GENIE team
    syn.setPermissions(mafDatabaseSynId, 3326313, [])

#Validation of all center files
def validation(syn, center, process, center_mapping_df, databaseToSynIdMappingDf, thread, testing, oncotreeLink):
    centerInputSynId = center_mapping_df['inputSynId'][center_mapping_df['center'] == center][0]
    logger.info("Center: " + center)
    allFiles = getCenterInputFiles(syn, centerInputSynId, center, process)

    allFiles = pd.DataFrame(allFiles,columns=['synId','filePaths'])
    #If a center has no files, then return empty list
    if allFiles.empty:
        logger.info("%s has not uploaded any files" % center)
        return([])
    else:
        #Make sure the vcf validation statuses don't get wiped away
        if process != "vcf":
            addToQuery = "and name not like '%.vcf'"
        else:
            addToQuery = ''
        validationStatus = syn.tableQuery("SELECT * FROM %s where center = '%s' %s" % (process_functions.getDatabaseSynId(syn, "validationStatus", databaseToSynIdMappingDf=databaseToSynIdMappingDf), center, addToQuery))
        errorTracker = syn.tableQuery("SELECT * FROM %s where center = '%s' %s"  % (process_functions.getDatabaseSynId(syn, "errorTracker", databaseToSynIdMappingDf=databaseToSynIdMappingDf), center,addToQuery))
        #VALIDATE FILES
        validationStatusDf = validationStatus.asDataFrame()
        errorTrackerDf = errorTracker.asDataFrame()
        validated = allFiles.apply(lambda x: validateFile(syn, validationStatusDf, errorTrackerDf, center, thread, x, testing, oncotreeLink), axis=1)
        inputValidStatus = []
        invalidErrors = []
        for inputStat, invalErrors in validated:
            inputValidStatus.extend(inputStat)
            if invalErrors is not None:
                invalidErrors.extend(invalErrors)
        inputValidStatus = pd.DataFrame(inputValidStatus, columns = ["id",'path','md5','status','name','modifiedOn','fileType'])
        logger.info("CHECK FOR DUPLICATED FILES")
        ##### DUPLICATED FILES ######
        #check for duplicated filenames.  There should be no duplication, files should be uploaded as new versions and the entire dataset should be uploaded everytime
        #cbs and seg files should not be duplicated.  There can only be one
        duplicatedFiles = inputValidStatus[inputValidStatus['name'].duplicated(keep=False)]
        cbsSegBool = [os.path.basename(i).endswith('.cbs') or os.path.basename(i).endswith('.seg') for i in inputValidStatus['name']]
        cbsSegFiles = inputValidStatus[cbsSegBool]
        if len(cbsSegFiles) >1:
            duplicatedFiles = duplicatedFiles.append(cbsSegFiles)
        # nodups = ["data_mutations_extended"]
        # allDuplicatedFiles = []
        # for nodup in nodups:
        #   checkDups = [name for name in inputValidStatus['name'] if name.startswith(nodup)]
        #   if len(checkDups) > 1:
        #       allDuplicatedFiles.extend(checkDups)
        # duplicatedFiles = duplicatedFiles.append(inputValidStatus[inputValidStatus['name'].isin(allDuplicatedFiles)])

        duplicatedFiles.drop_duplicates("id",inplace=True)
        inputValidStatus['status'][inputValidStatus['id'].isin(duplicatedFiles['id'])] = "INVALID"
        duplicatedFiles['errors'] = "DUPLICATED FILENAME! FILES SHOULD BE UPLOADED AS NEW VERSIONS AND THE ENTIRE DATASET SHOULD BE UPLOADED EVERYTIME"
        #Send an email if there are any duplicated files
        if not duplicatedFiles.empty:
            incorrectFiles = ", ".join([name for synId, name in zip(duplicatedFiles['id'],duplicatedFiles['name'])])
            incorrectEnt = syn.get(duplicatedFiles['id'].iloc[0])
            sendEmail = set([incorrectEnt.modifiedBy, incorrectEnt.createdBy])
            userNames = ", ".join([syn.getUserProfile(user).userName for user in sendEmail])
            syn.sendMessage(list(sendEmail), "GENIE Validation Error", "Dear %s,\n\nYour files (%s) are duplicated!  FILES SHOULD BE UPLOADED AS NEW VERSIONS AND THE ENTIRE DATASET SHOULD BE UPLOADED EVERYTIME" % (userNames, incorrectFiles))
        logger.info("THERE ARE %d DUPLICATED FILES" % len(duplicatedFiles))
        ##### DUPLICATED FILES ######

        #Create invalid error synapse table
        logger.info("UPDATE INVALID FILE REASON DATABASE")
        invalidErrors = pd.DataFrame(invalidErrors, columns = ["id",'errors','name'])
        # Remove fixed duplicated files
        dupIds = invalidErrors['id'][invalidErrors['errors'] == "DUPLICATED FILENAME! FILES SHOULD BE UPLOADED AS NEW VERSIONS AND THE ENTIRE DATASET SHOULD BE UPLOADED EVERYTIME"]
        removeIds = dupIds[~dupIds.isin(duplicatedFiles['id'])]
        invalidErrors = invalidErrors[~invalidErrors['id'].isin(removeIds)]
        # Append duplicated file errors
        invalidErrors = invalidErrors.append(duplicatedFiles[['id','errors','name']])
        invalidErrors['center'] = center
        invalidIds = inputValidStatus['id'][inputValidStatus['status'] == "INVALID"]
        invalidErrors = invalidErrors[invalidErrors['id'].isin(invalidIds)]
        process_functions.updateDatabase(syn, errorTracker.asDataFrame(), invalidErrors, process_functions.getDatabaseSynId(syn, "errorTracker", databaseToSynIdMappingDf=databaseToSynIdMappingDf), ["id"], to_delete=True)

        paths = inputValidStatus['path']
        filenames = [os.path.basename(name) for name in paths]
        del inputValidStatus['path']
        logger.info("UPDATE VALIDATION STATUS DATABASE")
        inputValidStatus['center'] = center
        #Remove fixed duplicated files
        inputValidStatus = inputValidStatus[~inputValidStatus['id'].isin(removeIds)]

        process_functions.updateDatabase(syn, validationStatus.asDataFrame(), inputValidStatus[["id",'md5','status','name','center','modifiedOn']], process_functions.getDatabaseSynId(syn, "validationStatus", databaseToSynIdMappingDf=databaseToSynIdMappingDf), ["id"], to_delete=True)
        inputValidStatus['path'] = paths
        validFiles = inputValidStatus[['id','path','fileType']][inputValidStatus['status'] == "VALIDATED"]
        return(validFiles)

def input_to_database(syn, center, process, testing, only_validate, vcf2maf_path, vep_path, vep_data, database_to_synid_mappingdf, center_mapping_df, reference=None, delete_old=False, oncotree_link=None, thread=1):
    if only_validate:
        log_path = os.path.join(process_functions.SCRIPT_DIR, "%s_validation_log.txt" % center)
    else:
        log_path = os.path.join(process_functions.SCRIPT_DIR, "%s_%s_log.txt" % (center, process))

    logFormatter = logging.Formatter("%(asctime)s [%(name)s][%(levelname)s] %(message)s")
    fileHandler = logging.FileHandler(log_path,mode='w')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    if testing:
        logger.info("###########################################")
        logger.info("############NOW IN TESTING MODE############")
        logger.info("###########################################")

    # ----------------------------------------
    # Start input to staging process
    # ----------------------------------------

    #path_to_genie = os.path.realpath(os.path.join(process_functions.SCRIPT_DIR,"../"))
    #Make the synapsecache dir the genie input folder for now
    #The main reason for this is because the .synaspecache dir is mounted by batch
    path_to_genie = os.path.expanduser("~/.synapseCache")
    #Create input and staging folders
    if not os.path.exists(os.path.join(path_to_genie,center,"input")):
        os.makedirs(os.path.join(path_to_genie,center,"input"))
    if not os.path.exists(os.path.join(path_to_genie,center,"staging")):
        os.makedirs(os.path.join(path_to_genie,center,"staging"))
     
    if delete_old:
        process_functions.rmFiles(os.path.join(path_to_genie,center))

    validFiles = validation(syn, center, process, center_mapping_df, database_to_synid_mappingdf, thread, testing, oncotree_link)

    if len(validFiles) > 0 and not only_validate:
        #Reorganize so BED file are always validated and processed first
        validBED = [os.path.basename(i).endswith('.bed') for i in validFiles['path']]
        beds = validFiles[validBED]
        validFiles = beds.append(validFiles)
        validFiles.drop_duplicates(inplace=True)
        #Valid vcf files
        validVCF = [i for i in validFiles['path'] if os.path.basename(i).endswith('.vcf')]
        #validCBS = [i for i in validFiles['path'] if os.path.basename(i).endswith('.cbs')]

        processTrackerSynId = process_functions.getDatabaseSynId(syn, "processTracker", databaseToSynIdMappingDf = database_to_synid_mappingdf)
        #Add process tracker for time start
        processTracker = syn.tableQuery("SELECT timeStartProcessing FROM %s where center = '%s' and processingType = '%s'" % (processTrackerSynId, center, process))
        processTrackerDf = processTracker.asDataFrame()
        if len(processTrackerDf) == 0:
            new_rows = [[center,str(int(time.time()*1000)), str(int(time.time()*1000)), process]] 
            table = syn.store(synapseclient.Table(processTrackerSynId, new_rows))
        else:
            processTrackerDf['timeStartProcessing'][0] = str(int(time.time()*1000))
            syn.store(synapseclient.Table(processTrackerSynId,processTrackerDf))

        processFiles(syn, validFiles, center, path_to_genie, thread, 
                     center_mapping_df, oncotree_link, database_to_synid_mappingdf, 
                     validVCF=validVCF,
                     vcf2mafPath=vcf2maf_path,
                     veppath=vep_path,vepdata=vep_data,
                     test=testing, processing=process,reference=reference)

        #Should add in this process end tracking before the deletion of samples
        processTracker = syn.tableQuery("SELECT timeEndProcessing FROM %s where center = '%s' and processingType = '%s'" % (processTrackerSynId, center, process))
        processTrackerDf = processTracker.asDataFrame()
        processTrackerDf['timeEndProcessing'][0] = str(int(time.time()*1000))
        syn.store(synapseclient.Table(processTrackerSynId,processTrackerDf))

        logger.info("SAMPLE/PATIENT RETRACTION")
        toRetract.retract(syn, testing)

    else:
        messageOut = "%s does not have any valid files" if not only_validate else "ONLY VALIDATION OCCURED FOR %s"
        logger.info(messageOut % center)

    # Store log file
    log_folder_synid = process_functions.getDatabaseSynId(
        syn, "logs", databaseToSynIdMappingDf=database_to_synid_mappingdf)
    syn.store(synapseclient.File(log_path, parentId=log_folder_synid))
    os.remove(log_path)
    logger.info("ALL PROCESSES COMPLETE")


def main():
    """Set up argument parser and returns"""
    parser = argparse.ArgumentParser(
        description='GENIE center inputs to database')
    parser.add_argument("process",choices=['vcf','maf','main','mafSP'],
                        help='Process vcf, maf or the rest of the files')
    parser.add_argument('--center', help='The centers')
    parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")
    parser.add_argument("--deleteOld", action='store_true', help = "Delete all old processed and temp files")
    parser.add_argument("--onlyValidate", action='store_true', help = "Only validate the files, don't process")
    parser.add_argument("--oncotreeLink", type=str, help="Link to oncotree code")
    parser.add_argument("--createNewMafDatabase", action='store_true', help = "Creates a new maf database")
    parser.add_argument("--testing", action='store_true', help = "Testing the infrastructure!")
    parser.add_argument("--debug", action='store_true', help = "Add debug mode to synapse")
    parser.add_argument("--reference", type=str, help="Path to VCF reference file")

    #DEFAULT PARAMS
    parser.add_argument("--vcf2mafPath", type=str, help="Path to vcf2maf", default="~/vcf2maf-1.6.14")
    parser.add_argument("--vepPath", type=str, help="Path to VEP", default="~/vep")
    parser.add_argument("--vepData", type=str, help="Path to VEP data", default="~/.vep")
    parser.add_argument('--thread', type=int, help="Number of threads to use for validation", default=1)

    args = parser.parse_args()
    syn = process_functions.synLogin(args.pemFile, debug=args.debug)
    #Must specify path to vcf2maf, VEP and VEP data is these types are specified
    if args.process in ['vcf','maf','mafSP'] and not args.onlyValidate:
        assert os.path.exists(args.vcf2mafPath), "Path to vcf2maf (--vcf2mafPath) must be specified if `--process {vcf,maf,mafSP}` is used"
        assert os.path.exists(args.vepPath), "Path to VEP (--vepPath) must be specified if `--process {vcf,maf,mafSP}` is used"
        assert os.path.exists(args.vepData), "Path to VEP data (--vepData) must be specified if `--process {vcf,maf,mafSP}` is used"
    
    if args.testing:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn11600968')
    else:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn10967259')

    databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()

    center_mapping_id = process_functions.getDatabaseSynId(syn, "centerMapping", databaseToSynIdMappingDf=databaseToSynIdMappingDf)
    center_mapping = syn.tableQuery('SELECT * FROM %s' % center_mapping_id)
    center_mapping_df = center_mapping.asDataFrame()

    if args.center is not None:
        assert args.center in center_mapping_df.center.tolist(), "Must specify one of these centers: %s" % ", ".join(center_mapping_df.center)
        centers = [args.center]
    else:
        center_mapping_df = center_mapping_df[~center_mapping_df['inputSynId'].isnull()]
        center_mapping_df = center_mapping_df[center_mapping_df['release'] == True]
        centers = center_mapping_df.center
    
    if args.oncotreeLink is None:
        onco_link = databaseToSynIdMappingDf['Id'][databaseToSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
        onco_link_ent = syn.get(onco_link)
        args.oncotreeLink = onco_link_ent.externalURL
    #Check if you can connect to oncotree link, if not then don't run validation / processing
    process_functions.checkUrl(args.oncotreeLink)

    center_mapping_ent = syn.get(center_mapping_id)
    if center_mapping_ent.get('isProcessing',['True'])[0] == 'True':
        raise Exception("Processing/validation is currently happening.  Please change/add the 'isProcessing' annotation on %s to False to enable processing" % center_mapping_id)
    else:
        center_mapping_ent.isProcessing="True"
        center_mapping_ent = syn.store(center_mapping_ent)
    #remove this query timeout and see what happens
    #syn.table_query_timeout = 50000

    #Create new maf database, should only happen once if its specified
    if args.createNewMafDatabase:
        createMafDatabase(syn, databaseToSynIdMappingDf, testing=args.testing)

    for center in centers:
        input_to_database(syn, center, args.process, args.testing, args.onlyValidate, args.vcf2mafPath, args.vepPath, args.vepData, databaseToSynIdMappingDf, center_mapping_df, reference=args.reference, delete_old=args.deleteOld, oncotree_link=args.oncotreeLink, thread=args.thread)

    # To ensure that this is the new entity
    center_mapping_ent = syn.get(center_mapping_id)
    center_mapping_ent.isProcessing="False"
    center_mapping_ent = syn.store(center_mapping_ent)

    error_tracker_synid = process_functions.getDatabaseSynId(syn, "errorTracker", databaseToSynIdMappingDf=databaseToSynIdMappingDf)
    #Only write out invalid reasons if the center isnt specified and if only validate
    if args.center is None and args.onlyValidate:
        logging.info("WRITING INVALID REASONS TO CENTER STAGING DIRS")
        write_invalid_reasons.write_invalid_reasons(syn, center_mapping_df, error_tracker_synid)

if __name__ == "__main__":
    main()
    
