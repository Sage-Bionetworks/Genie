import argparse
import datetime
import logging
import os
import re
import shutil
import subprocess
import time

import synapseclient
import synapseutils
import pandas as pd

from . import process_functions
from . import database_to_staging
from . import create_case_lists
from . import dashboard_table_updater

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def storeFile(syn, filePath, parentId, anonymizeCenterDf,
              genie_version, name=None):
    #process.center_anon(filePath, anonymizeCenterDf)
    if name is None:
        name = os.path.basename(filePath)
    return(syn.store(synapseclient.File(filePath, name=name, parent = parentId, versionComment=genie_version)))
    #process.center_convert_back(filePath, anonymizeCenterDf)

#This is the only filter that returns mutation columns to keep
def commonVariantFilter(mafDf):
    mafDf['FILTER'] = mafDf['FILTER'].fillna("")
    toKeep = ["common_variant" not in i for i in mafDf['FILTER']]
    mafDf = mafDf[toKeep]
    return(mafDf)

def consortiumToPublic(syn, processingDate, genie_version, releaseId, databaseSynIdMappingDf, publicReleaseCutOff=365, staging=False):

    ANONYMIZE_CENTER = syn.tableQuery('SELECT * FROM syn10170510')
    ANONYMIZE_CENTER_DF = ANONYMIZE_CENTER.asDataFrame()
    CNA_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,"data_CNA_%s.txt" % genie_version)
    CLINICAL_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,'data_clinical_%s.txt' % genie_version)
    CLINICAL_SAMPLE_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,'data_clinical_sample_%s.txt' % genie_version)
    CLINICAL_PATIENT_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,'data_clinical_patient_%s.txt' % genie_version)
    DATA_GENE_PANEL_PATH = os.path.join(dbTodatabase_to_stagingstaging.GENIE_RELEASE_DIR,'data_gene_matrix_%s.txt' % genie_version)
    MUTATIONS_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,'data_mutations_extended_%s.txt' % genie_version)
    FUSIONS_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,'data_fusions_%s.txt' % genie_version)
    SEG_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,'genie_public_data_cna_hg19_%s.seg' % genie_version)
    COMBINED_BED_PATH = os.path.join(database_to_staging.GENIE_RELEASE_DIR,'genomic_information_%s.txt' % genie_version)

    if not os.path.exists(database_to_staging.GENIE_RELEASE_DIR):
        os.mkdir(database_to_staging.GENIE_RELEASE_DIR)
    if not os.path.exists(database_to_staging.CASE_LIST_PATH):
        os.mkdir(database_to_staging.CASE_LIST_PATH)

    # if staging:
    #   #public release staging
    #   PUBLIC_RELEASE_PREVIEW = "syn7871696"
    #   PUBLIC_RELEASE_PREVIEW_CASELIST = "syn9689659"
    # else:
    #public release preview
    PUBLIC_RELEASE_PREVIEW =  databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'public'].values[0]
    PUBLIC_RELEASE_PREVIEW_CASELIST = database_to_staging.find_caselistid(syn, PUBLIC_RELEASE_PREVIEW)

    ##############################################################################################################################
    ## Sponsored projects filter
    ##############################################################################################################################
    ## if before release date -> go into staging consortium
    ## if after date -> go into public 
    # sponsoredReleaseDate = syn.tableQuery('SELECT * FROM syn8545108')
    # sponsoredReleaseDateDf = sponsoredReleaseDate.asDataFrame()
    # sponsoredProjectSamples = syn.tableQuery('SELECT * FROM syn8545106')
    # sponsoredProjectSamplesDf = sponsoredProjectSamples.asDataFrame()
    # sponsoredProjectsDf = sponsoredProjectSamplesDf.merge(sponsoredReleaseDateDf, left_on="sponsoredProject", right_on="sponsoredProjects")
    # dates = sponsoredProjectsDf['releaseDate'].apply(lambda date: datetime.datetime.strptime(date, '%b-%Y'))
    # publicReleaseSamples = sponsoredProjectsDf['genieSampleId'][dates < processingDate]
    ##############################################################################################################################
    
    # SEQ_DATE filter
    # Jun-2015, given processing date (today) -> public release (processing date - Jun-2015 > 12 months)
    consortiumReleaseWalk = synapseutils.walk(syn, releaseId)

    consortiumRelease = next(consortiumReleaseWalk)
    clinical = [syn.get(synid, followLink=True) for filename, synid in consortiumRelease[2] if filename == "data_clinical.txt"][0]
    gene_matrix = [syn.get(synid, followLink=True) for filename, synid in consortiumRelease[2] if filename == "data_gene_matrix.txt"][0]

    clinicalDf = pd.read_csv(clinical.path, sep="\t", comment="#")
    gene_matrixdf = pd.read_csv(gene_matrix.path, sep="\t")

    removeForPublicSamples = process_functions.seqDateFilter(clinicalDf,processingDate,publicReleaseCutOff)
    #comment back in when public release filter back on
    #publicReleaseSamples = publicReleaseSamples.append(keepForPublicSamples)
    #Make sure all null oncotree codes are removed
    clinicalDf = clinicalDf[~clinicalDf['ONCOTREE_CODE'].isnull()]
    publicReleaseSamples = clinicalDf.SAMPLE_ID[~clinicalDf.SAMPLE_ID.isin(removeForPublicSamples)]

    logger.info("SEQ_DATES for public release: " + ", ".join(set(clinicalDf.SEQ_DATE[clinicalDf.SAMPLE_ID.isin(publicReleaseSamples)].astype(str))))

    #Clinical release scope filter
    #If consortium -> Don't release to public
    clinicalReleaseScope = syn.tableQuery("SELECT * FROM syn8545211 where releaseScope = 'public'")
    publicRelease = clinicalReleaseScope.asDataFrame()

    allClin = clinicalDf[clinicalDf['SAMPLE_ID'].isin(publicReleaseSamples)]
    allClin.to_csv(CLINICAL_PATH, sep="\t", index=False)

    gene_matrixdf = gene_matrixdf[gene_matrixdf['SAMPLE_ID'].isin(publicReleaseSamples)]
    gene_matrixdf.to_csv(DATA_GENE_PANEL_PATH,sep="\t",index=False)
    storeFile(syn, DATA_GENE_PANEL_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_gene_matrix.txt")
    storeFile(syn, CLINICAL_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_clinical.txt")
    
    create_case_lists.main(CLINICAL_PATH, DATA_GENE_PANEL_PATH, database_to_staging.CASE_LIST_PATH, "genie_public")

    caseListFiles = os.listdir(database_to_staging.CASE_LIST_PATH)
    caseListEntities = []
    for casePath in caseListFiles:
        casePath = os.path.join(database_to_staging.CASE_LIST_PATH, casePath)
        caseListEntities.append(storeFile(syn, casePath, PUBLIC_RELEASE_PREVIEW_CASELIST, ANONYMIZE_CENTER_DF, genie_version))

    #Grab mapping table to fill in clinical headers
    mapping_table = syn.tableQuery('SELECT * FROM syn9621600')
    mapping = mapping_table.asDataFrame()
    genePanelEntities = []
    for entName, entId in consortiumRelease[2]:
        if "data_linear" in entName or "meta_" in entName:
            continue
        elif entName == "data_clinical.txt":
            patientCols = publicRelease['fieldName'][publicRelease['level'] == "patient"].tolist()
            sampleCols = ["PATIENT_ID"]
            sampleCols.extend(publicRelease['fieldName'][publicRelease['level'] == "sample"].tolist())
            #clinicalDf is defined on line 36
            # clinicalDf['AGE_AT_SEQ_REPORT'] = [int(math.floor(int(float(i))/365.25)) if process.checkInt(i) else i for i in clinicalDf['AGE_AT_SEQ_REPORT']]
            # clinicalDf['AGE_AT_SEQ_REPORT'][clinicalDf['AGE_AT_SEQ_REPORT'] == ">32485"] = ">89"
            # clinicalDf['AGE_AT_SEQ_REPORT'][clinicalDf['AGE_AT_SEQ_REPORT'] == "<6570"] = "<18"

            clinicalDf = clinicalDf[clinicalDf['SAMPLE_ID'].isin(publicReleaseSamples)]

            #Delete columns that are private scope
            # for private in privateRelease:
            #   del clinicalDf[private]
            process_functions.addClinicalHeaders(clinicalDf, mapping, patientCols, sampleCols, CLINICAL_SAMPLE_PATH, CLINICAL_PATIENT_PATH)

            storeFile(syn, CLINICAL_SAMPLE_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_clinical_sample.txt")
            storeFile(syn, CLINICAL_PATIENT_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_clinical_patient.txt")

        elif "mutation" in entName:
            mutation = syn.get(entId, followLink=True)
            mutationDf = pd.read_csv(mutation.path, sep="\t", comment="#")
            mutationDf = commonVariantFilter(mutationDf)
            mutationDf['FILTER'] = "PASS"
            mutationDf = mutationDf[mutationDf['Tumor_Sample_Barcode'].isin(publicReleaseSamples)]
            text = process_functions.removeFloat(mutationDf)
            with open(MUTATIONS_PATH, 'w') as f:
                f.write(text)
            storeFile(syn, MUTATIONS_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_mutations_extended.txt")

        elif "fusion" in entName:
            fusion = syn.get(entId, followLink=True)
            fusionDf = pd.read_csv(fusion.path, sep="\t")
            #remove = ["Entrez_Gene_Id","Method"]
            #fusionDf = fusionDf[fusionDf.columns[~fusionDf.columns.isin(remove)]]
            fusionDf = fusionDf[fusionDf['Tumor_Sample_Barcode'].isin(publicReleaseSamples)]
            fusionDf.to_csv(FUSIONS_PATH,sep="\t",index=False)
            storeFile(syn, FUSIONS_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_fusions.txt")
        elif "CNA" in entName:
            cna = syn.get(entId, followLink=True)
            cnaDf = pd.read_csv(cna.path, sep="\t")
            cnaDf = cnaDf[cnaDf.columns[cnaDf.columns.isin(publicReleaseSamples.append(pd.Series("Hugo_Symbol")))]]
            text = process_functions.removeFloat(cnaDf)
            text = text.replace("\t\t","\tNA\t").replace("\t\t","\tNA\t").replace('\t\n',"\tNA\n")
            with open(CNA_PATH, "w") as cnaFile:
                cnaFile.write(text)
            storeFile(syn, CNA_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_CNA.txt")
        elif entName.endswith(".seg"):
            seg = syn.get(entId, followLink=True)
            segDf = pd.read_csv(seg.path, sep="\t")
            segDf = segDf[segDf['ID'].isin(publicReleaseSamples)]
            text = process_functions.removeFloat(segDf)
            with open(SEG_PATH, "w") as segFile:
                segFile.write(text)
            storeFile(syn, SEG_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="genie_public_data_cna_hg19.seg")
        elif entName == "data_gene_matrix.txt":
            pass
            # This file was processed above because it had to be used for generating caselists
            # panel = syn.get(entId, followLink=True)
            # panelDf = pd.read_csv(panel.path, sep="\t")
            # panelDf = panelDf[panelDf['SAMPLE_ID'].isin(publicReleaseSamples)]
            # panelDf.to_csv(DATA_GENE_PANEL_PATH,sep="\t",index=False)
            # storeFile(syn, DATA_GENE_PANEL_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="data_gene_matrix.txt")
        elif entName == "genomic_information.txt":
            bed = syn.get(entId, followLink=True)
            bedDf = pd.read_csv(bed.path, sep="\t")
            bedDf = bedDf[bedDf.SEQ_ASSAY_ID.isin(allClin.SEQ_ASSAY_ID)]
            bedDf.to_csv(COMBINED_BED_PATH,sep="\t",index=False)
            storeFile(syn, COMBINED_BED_PATH, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name="genomic_information.txt")
        elif entName in ["data_clinical_sample.txt", "data_clinical_patient.txt"] or entName.endswith(".html"):
            continue
        elif entName.startswith("data_gene_panel"):
            genePanel = syn.get(entId, followLink=True)
            #Create new gene panel naming and store
            fileName = os.path.basename(genePanel.path)
            newFileList = fileName.split("_")
            newFileList[-1] = genie_version + ".txt"
            newFileName = "_".join(newFileList)
            genePanelPath = os.path.join(database_to_staging.GENIE_RELEASE_DIR, newFileName)
            shutil.copy(genePanel.path, genePanelPath)
            del newFileList[-1]
            entName = "_".join(newFileList)
            entName = entName + ".txt"
            genePanelEntities.append(storeFile(syn, genePanelPath, PUBLIC_RELEASE_PREVIEW, ANONYMIZE_CENTER_DF, genie_version, name=entName))
        else:
            ent = syn.get(entId, followLink=True, downloadFile=False)
            copiedId = synapseutils.copy(syn, ent, PUBLIC_RELEASE_PREVIEW, version=ent.versionNumber, updateExisting=True, setProvenance = None, skipCopyAnnotations=True)
            copiedEnt = syn.get(copiedId[ent.id],downloadFile=False)
            #Set version comment
            copiedEnt.versionComment=genie_version
            syn.store(copiedEnt, forceVersion=False)
    return((caseListEntities,genePanelEntities))

def perform_consortiumToPublic(syn, args, databaseSynIdMappingDf):
    try:
        processingDate = datetime.datetime.strptime(args.processingDate, '%b-%Y')
    except ValueError as e:
        raise ValueError("Process date must be in the format abbreviated_month-YEAR ie. Oct-2017")
    return(consortiumToPublic(syn, processingDate, args.genieVersion, args.releaseId, databaseSynIdMappingDf, publicReleaseCutOff = args.publicReleaseCutOff, staging = args.staging))

def command_reviseMetadataFiles(syn, args, databaseSynIdMappingDf):
    reviseMetadataFiles(syn, args.staging, databaseSynIdMappingDf, args.genieVersion)

def reviseMetadataFiles(syn, staging, databaseSynIdMappingDf, genieVersion=None):
    ANONYMIZE_CENTER = syn.tableQuery('SELECT * FROM syn10170510')
    ANONYMIZE_CENTER_DF = ANONYMIZE_CENTER.asDataFrame()
    # if staging:
    #   parent = "syn7871696"
    # else:
    parent =  databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'public'].values[0]
    allFiles = syn.getChildren(parent)
    metadataEnts = [syn.get(i['id'], downloadLocation=database_to_staging.GENIE_RELEASE_DIR, ifcollision="overwrite.local") for i in allFiles if 'meta' in i['name']]
    for metaEnt in metadataEnts:
        with open(metaEnt.path, "r+") as meta:
            metaText = meta.read()
            if "meta_study" not in metaEnt.path:
                version = ''
            else:
                version = re.search(".+GENIE.+v(.+)", metaText).group(1)
            #Fix this line
            genieVersion = version if genieVersion is None else genieVersion
            dataFileVersion = re.search(".+data_(.+)[.]txt",metaText)
            if dataFileVersion is None:
                dataFileVersion = re.search(".+data_(.+)[.]seg",metaText)
            if dataFileVersion is not None:
                dataFileVersion = dataFileVersion.group(1).split("_")[-1]

            if version != genieVersion:
                metaText = metaText.replace("GENIE Cohort v%s" % version,"GENIE Cohort v%s" % genieVersion)
                metaText = metaText.replace("GENIE v%s" % version,"GENIE v%s" % genieVersion)
                if dataFileVersion is not None:
                    metaText = metaText.replace(dataFileVersion, genieVersion)
                    metaText = metaText.replace(dataFileVersion, genieVersion)
                meta.seek(0)
                meta.write(metaText)
                meta.truncate()

                storeFile(syn, metaEnt.path, parent, ANONYMIZE_CENTER_DF, genieVersion)


def createLinkVersion(syn, genie_version, caseListEntities, genePanelEntities, databaseSynIdMappingDf):
    versioning = genie_version.split(".")
    logger.info(genie_version)
    main = versioning[0]
    releaseSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'release'].values[0]
    publicSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'public'].values[0]
    #second = ".".join(versioning[1:])
    releases = synapseutils.walk(syn, releaseSynId)
    mainReleaseFolders = next(releases)[1]
    releaseFolderSynId = [synId for folderName, synId in mainReleaseFolders if folderName == "Release %s" % main] 
    if len(releaseFolderSynId) > 0:
        secondRelease = synapseutils.walk(syn, releaseFolderSynId[0])
        secondReleaseFolders = next(secondRelease)[1]
        secondReleaseFolderSynIdList = [synId for folderName, synId in secondReleaseFolders if folderName == genie_version] 
        if len(secondReleaseFolderSynIdList) > 0:
            secondReleaseFolderSynId = secondReleaseFolderSynIdList[0]
        else:
            secondReleaseFolderSynId = syn.store(synapseclient.Folder(genie_version, parent = releaseFolderSynId[0])).id
    else:
        mainReleaseFolderId = syn.store(synapseclient.Folder("Release %s" % main, parent = releaseSynId)).id
        secondReleaseFolderSynId = syn.store(synapseclient.Folder(genie_version, parent = mainReleaseFolderId)).id

    caselistId = database_to_staging.find_caselistid(syn, secondReleaseFolderSynId)

    publicRelease = syn.getChildren(publicSynId)
    [syn.store(synapseclient.Link(ents['id'], parent=secondReleaseFolderSynId, targetVersion=ents['versionNumber'])) for ents in publicRelease if ents['type'] != "org.sagebionetworks.repo.model.Folder" and ents['name'] != "data_clinical.txt"  and not ents['name'].startswith("data_gene_panel")]
    [syn.store(synapseclient.Link(ents.id, parent=caselistId, targetVersion=ents.versionNumber)) for ents in caseListEntities]
    #Store gene panels
    [syn.store(synapseclient.Link(ents.id, parent=secondReleaseFolderSynId, targetVersion=ents.versionNumber)) for ents in genePanelEntities]


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("processingDate", type=str, metavar="Jan-2017",
                        help="The process date of GENIE in Month-Year format (ie. Apr-2017)")
    parser.add_argument("cbioportalPath", type=str, metavar="/path/to/cbioportal",
                        help="Make sure you clone the cbioportal github: git clone https://github.com/cBioPortal/cbioportal.git")
    parser.add_argument("genieVersion", type=str,
                        help="GENIE public release version")    
    parser.add_argument("--publicReleaseCutOff", type=int, metavar=366, default=366,
                        help="Public release cut off time in days (Must account for leap year, 366)")
    parser.add_argument("--staging", action='store_true',
                        help="Store into staging folder")
    parser.add_argument("--test", action='store_true',
                        help="Store into staging folder")
    parser.add_argument("--pemFile", type=str, 
                        help="Path to PEM file (genie.pem)")
    parser.add_argument("--debug", action='store_true', 
                        help="Synapse debug feature")
    args = parser.parse_args()
    cbioValidatorPath = os.path.join(args.cbioportalPath,"core/src/main/scripts/importer/validateData.py")
    assert os.path.exists(cbioValidatorPath), "Please specify correct cbioportalPath"
    assert not (args.test and args.staging), "You can only specify --test or --staging, not both"
        
    syn = process_functions.synLogin(args.pemFile, debug=args.debug)
    #Get all the possible public releases
    if args.test:
        databaseSynIdMappingId = 'syn11600968'
        args.genieVersion = "TESTpublic"
    elif args.staging:
        databaseSynIdMappingId = 'syn12094210'  
    else:
        databaseSynIdMappingId = 'syn10967259'  
    databaseSynIdMapping = syn.tableQuery('select * from %s' % databaseSynIdMappingId)
    databaseSynIdMappingDf = databaseSynIdMapping.asDataFrame()
    releaseSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'release'].values[0]
    
    temp = synapseutils.walk(syn, releaseSynId)
    officialPublic = dict()
    for dirpath, dirnames, filenames in temp:
        release = os.path.basename(dirpath[0])
        #checkRelease = release.split(".")
        final = [i.split("-") for i in release.split(".")]
        checkRelease = []
        for i in final:
            checkRelease.extend(i)
        if args.test:
            officialPublic['TESTpublic'] = "syn12299959"
        else:
            if len(checkRelease) == 3 and checkRelease[0] != "0":
                if int(checkRelease[1]) > 0:
                    if checkRelease[0] in ['1','2']:
                        officialPublic[str(int(checkRelease[0])+1)+".0.0"] = dirpath[1]
                    else:
                        officialPublic[str(int(checkRelease[0]))+".0-public"] = dirpath[1]
    assert args.genieVersion in officialPublic.keys(), "genieVersion must be one of these: %s." % ", ".join(officialPublic.keys())
    args.releaseId = officialPublic[args.genieVersion]
    if not args.test and not args.staging:
        processTrackerSynId = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'processTracker'].values[0]
        processTracker = syn.tableQuery("SELECT timeStartProcessing FROM %s where center = 'SAGE' and processingType = 'public'" % processTrackerSynId)
        processTrackerDf = processTracker.asDataFrame()
        processTrackerDf['timeStartProcessing'][0] = str(int(time.time()*1000))
        syn.store(synapseclient.Table(processTrackerSynId,processTrackerDf))

    caseListEntities, genePanelEntities = perform_consortiumToPublic(syn, args, databaseSynIdMappingDf)
    command_reviseMetadataFiles(syn, args, databaseSynIdMappingDf)
    logger.info("CBIO VALIDATION")
    #Must be exit 0 because the validator sometimes fails, but we still want to capture the output
    command = ['python',cbioValidatorPath, '-s', database_to_staging.GENIE_RELEASE_DIR, '-n','; exit 0']
    cbioOutput = subprocess.check_output(" ".join(command), shell=True)
    logger.info(cbioOutput.decode("utf-8"))
    if not args.test and not args.staging:
        log_folder_synid = databaseSynIdMappingDf['Id'][databaseSynIdMappingDf['Database'] == 'logs'].values[0]
        with open("cbioValidatorLogsPublic_%s.txt" % args.genieVersion, "w") as cbioLog:
            cbioLog.write(cbioOutput.decode("utf-8"))
        syn.store(synapseclient.File("cbioValidatorLogsPublic_%s.txt" % args.genieVersion, parentId=log_folder_synid))
        os.remove("cbioValidatorLogsPublic_%s.txt" % args.genieVersion)
    logger.info("REMOVING OLD FILES")
    process_functions.rmFiles(database_to_staging.CASE_LIST_PATH)
    if os.path.exists('%s/genie_public_meta_cna_hg19_seg.txt' % database_to_staging.GENIE_RELEASE_DIR):
        os.unlink('%s/genie_public_meta_cna_hg19_seg.txt' % database_to_staging.GENIE_RELEASE_DIR)

    logger.info("CREATING LINK VERSION")
    createLinkVersion(syn, args.genieVersion, caseListEntities, genePanelEntities, databaseSynIdMappingDf)
    #Don't update process tracker is testing or staging
    if not args.test and not args.staging:
        processTracker = syn.tableQuery("SELECT timeEndProcessing FROM %s where center = 'SAGE' and processingType = 'public'" % processTrackerSynId)
        processTrackerDf = processTracker.asDataFrame()
        processTrackerDf['timeEndProcessing'][0] = str(int(time.time()*1000))
        syn.store(synapseclient.Table(processTrackerSynId,processTrackerDf))

    if not args.test:
        logger.info("DASHBOARD UPDATE")
        dashboard_table_updater.run_dashboard(syn, databaseSynIdMappingDf, args.genieVersion, staging=args.staging, public=True)
        dashboard_markdown_html_commands = ['Rscript', os.path.join(os.path.dirname(os.path.abspath(__file__)),'dashboard_markdown_generator.R'), args.genieVersion]
        if args.staging:
            dashboard_markdown_html_commands.append('--staging')
        subprocess.check_call(dashboard_markdown_html_commands)
        logger.info("DASHBOARD UPDATE COMPLETE")

    logger.info("COMPLETED CONSORTIUM TO PUBLIC")
