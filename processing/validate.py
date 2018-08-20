#!/usr/bin/env python
import pandas as pd
assert pd.__version__ >= "0.20.0", "Please make sure your pandas version is at least 0.20.0.  If not, please do pip install pandas --upgrade"
import synapseclient
import os
#import subprocess
import argparse
import getpass
import string
import httplib2 as http
import json
from multiprocessing import Pool
import datetime
from functools import partial
import re
import warnings
from config_process_scripts import *
warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)
try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Validates annotations on Synapse
# def validateAnnotations(fileList):
#     logger.info("VALIDATING ANNOTATIONS")
#     notcorrect = []
#     for i,ID in enumerate(fileList['entity.id']):
#         foo = syn.get(ID, downloadFile=False)
#         required_annot = ["center","dataType","fileType","disease","consortium",
#         "platform","tissueSource","organism","dataSubType"]
#         check = [annot for annot in required_annot if foo.annotations.has_key(annot)]
#         if len(check) != len(required_annot):
#             notcorrect.append(fileList.iloc[i]['entity.id'])
#     if len(notcorrect) >0:
#         return(False)
#     else:
#         return(True)

def synapse_login():
    """
    This function logs into synapse for you if credentials are saved.  
    If not saved, then user is prompted username and password.

    :returns:     Synapseclient object
    """
    try:
        syn = synapseclient.login(silent=True)
    except Exception as e:
        logger.info("Please provide your synapse username/email and password (You will only be prompted once), or write your username and password to your bash profile under GENIE_USER and GENIE_PASS")
        Username = os.getenv("GENIE_USER")
        Password = os.getenv("GENIE_PASS")
        if Username is None or Password is None:
            Username = raw_input("Username: ")
            Password = getpass.getpass()
        syn = synapseclient.login(email=Username, password=Password,rememberMe=True,silent=True)
    return(syn)

def main(syn, fileType, filePath, center, threads, oncotreeLink="http://oncotree.mskcc.org/oncotree/api/tumor_types.txt?version=oncotree_2017_06_21", offline=False, uploadToSynapse=None, testing=False, noSymbolCheck=False):
    """
    This performs the validation of files

    :returns:   Text with the errors of the chosen file
    """
    #CHECK: Fail if filename is incorrect
    if not offline:
        try:
            PROCESS_FILES[fileType](syn, center, threads).validateFilename(filePath)
           # PROCESS_FILES[fileType].validateFilename(filePath, center)
        except AssertionError as e:
            raise ValueError("Your filename is incorrect!\n%s\nPlease change your filename before you run the validator again."  % e)
    total_error, warning = PROCESS_FILES[fileType](syn, center, threads).validate(filePathList=filePath, oncotreeLink=oncotreeLink, testing=testing, noSymbolCheck=noSymbolCheck)

    #Complete error message
    message = "----------------ERRORS----------------\n"
    if total_error == "":
        message = "YOUR FILE IS VALIDATED!\n"
        logger.info(message)
        valid = True
    else:
        for errors in total_error.split("\n"):
            if errors!='':
                logger.error(errors)
        message += total_error
        valid = False
    if warning != "":
        for warn in warning.split("\n"):  
            if warn!='':      
                logger.warning(warn)
        message += "-------------WARNINGS-------------\n" + warning
    if valid and uploadToSynapse is not None:
        logger.info("Uploading file to %s" % uploadToSynapse)
        [syn.store(synapseclient.File(path, parent=uploadToSynapse)) for path in filePath]
    return(message, valid)


def perform_main(syn, args):
    #pool = Pool(args.thread)
    message = main(syn, args.fileType, args.file, args.center, args.thread, args.oncotreeLink, args.offline, args.uploadToSynapse, args.testing, args.noSymbolCheck)
    #pool.close()
    #pool.join()

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Validate GENIE files')

    parser.add_argument("fileType", type=str, choices = PROCESS_FILES.keys(),
                        help='Filetypes that you are validating. Note, the filetypes with SP at the end are for special sponsored projects')
    parser.add_argument("file", type=str, nargs="+",
                        help='File(s) that you are validating.  If you validation your clinical files and you have both sample and patient files, you must provide both')
    parser.add_argument("center", type=str,
                        help='Contributing Centers')
    parser.add_argument("--thread", type=int, required=False, default=1,
                        help='Number of threads used in validation symbols')
    parser.add_argument("--offline", action='store_true',
                        help='No validation of filenames')
    parser.add_argument("--uploadToSynapse", type=str, default=None,
                        help='Will upload the file to the synapse directory of users choice')
    parser.add_argument("--oncotreeLink", type=str,
                        help="Link to oncotree code")
    parser.add_argument("--noSymbolCheck", action='store_true',
                        help='Do not check hugo symbols of fusion and cna file')
    parser.add_argument("--testing", action='store_true',
                        help='Put in testing mode')
    args = parser.parse_args()

    syn = synapse_login()

    if args.testing:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn11600968')
    else:
        databaseToSynIdMapping = syn.tableQuery('SELECT * FROM syn10967259')

    databaseToSynIdMappingDf = databaseToSynIdMapping.asDataFrame()
    synId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "centerMapping"]
    center_mapping = syn.tableQuery('SELECT * FROM %s' % synId[0])
    center_mapping_df = center_mapping.asDataFrame()
    assert args.center in center_mapping_df.center.tolist(), "Must specify one of these centers: %s" % ", ".join(center_mapping_df.center)
    
    if args.oncotreeLink is None:
        oncoLink = databaseToSynIdMappingDf['Id'][databaseToSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
        oncoLinkEnt = syn.get(oncoLink)
        args.oncotreeLink = oncoLinkEnt.externalURL

    if args.uploadToSynapse is not None:
        if args.offline:
            raise ValueError("If you specify the uploadToSynapse option, your filename must be named correctly")
        else:
            try:
                syn.get(args.uploadToSynapse)
            except synapseclient.exceptions.SynapseHTTPError as e:
                raise ValueError("Provided Synapse id must be your input folder Synapse id or a Synapse Id of a folder inside your input directory")
    perform_main(syn, args)
