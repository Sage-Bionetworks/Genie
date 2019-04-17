#!/usr/local/bin/python
import pandas as pd
import argparse
import synapseclient
import os
import process_functions as process
import logging
import math
import datetime
import time
import re
import subprocess
import synapseutils as synu
import create_case_lists
import dashboard_table_updater

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# GENIE CONSTANTS
# GENIE_RELEASE_DIR = os.path.join(
#    os.path.dirname(os.path.abspath(__file__)),"GENIE_release")
GENIE_RELEASE_DIR = os.path.join(
    os.path.expanduser("~/.synapseCache"), "GENIE_release")
CASE_LIST_PATH = os.path.join(
    GENIE_RELEASE_DIR, 'case_lists')
CNA_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, "data_CNA_%s.txt")
SAMPLE_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, 'data_clinical_supp_sample_%s.txt')
PATIENT_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, 'data_clinical_supp_patient_%s.txt')
MUTATIONS_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, 'data_mutations_extended_%s.txt')
FUSIONS_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, 'data_fusions_%s.txt')
SEG_CENTER_PATH = os.path.join(
    GENIE_RELEASE_DIR, 'genie_data_cna_hg19_%s.seg')
BED_DIFFS_SEQASSAY_PATH = os.path.join(
    GENIE_RELEASE_DIR, 'diff_%s.csv')


def findCaseListId(syn, parentId):
    '''
    Search for case_lists folder based on parentId given

    Args:
        syn: Synapse object
        parentId: Synapse Id of Folder or Project
    '''
    releaseEnts = synu.walk(syn, parentId)
    releaseFolders = next(releaseEnts)
    if len(releaseFolders[1]) == 0:
        caselistId = syn.store(
            synapseclient.Folder(name="case_lists", parent=parentId)).id
    else:
        caselistId = releaseFolders[1][0][1]
    return(caselistId)


def storeFile(
        syn,
        filePath,
        genieVersion="database",
        name=None,
        parent=None,
        fileFormat=None,
        cBioFileFormat=None,
        staging=False,
        caseLists=False,
        centerStaging=False):
    '''
    Convenience function to store Files
    Take care of clinical case (Clinical files go elsewhere)
    '''
    # ANONYMIZE_CENTER = syn.tableQuery('SELECT * FROM syn10170510')
    # ANONYMIZE_CENTER_DF = ANONYMIZE_CENTER.asDataFrame()

    if staging:
        parent = "syn9689654" if caseLists else "syn7551278"
    logger.info("STORING FILE: %s" % os.path.basename(filePath))
    # if not centerStaging:
    #   process.center_anon(filePath, ANONYMIZE_CENTER_DF)
    if name is None:
        name = os.path.basename(filePath)
    temp = synapseclient.File(
        filePath, name=name, parent=parent, versionComment=genieVersion)
    if fileFormat is not None:
        temp.fileFormat = fileFormat
    if cBioFileFormat is not None:
        temp.cBioFileFormat = cBioFileFormat
    temp = syn.store(temp)
    # if not centerStaging:
    #   process.center_convert_back(filePath, ANONYMIZE_CENTER_DF)
    return(temp)


def redaction_phi(values):
    '''
    Boolean vector of pediatric and PHI rows

    Args:
        values: Dataframe column or list of values to check

    Returns:
        tuple: redact and pediatric redaction bool vectors
    '''
    phi_cutoff = 365*89
    pediatric_cutoff = 365*18
    # Some sites submit redacted values already
    values = [
        pediatric_cutoff - 1
        if "<" in str(value) else value
        for value in values]
    values = [
        phi_cutoff + 1
        if ">" in str(value) else value
        for value in values]
    to_redact = [
        int(float(value)) > phi_cutoff
        if value not in ['', 'Unknown', 'Not Applicable', 'Not Collected']
        else False for value in values]
    to_redact_pediatric = [
        int(float(value)) < pediatric_cutoff
        if value not in ['', 'Unknown', 'Not Applicable', 'Not Collected']
        else False for value in values]
    return(to_redact, to_redact_pediatric)


def reAnnotatePHI(mergedClinical):
    # Remove PHI data
    mergedClinical['AGE_AT_SEQ_REPORT'] = \
        mergedClinical['AGE_AT_SEQ_REPORT'].fillna('')
    mergedClinical['BIRTH_YEAR'] = mergedClinical['BIRTH_YEAR'].fillna('')
    toRedact, toRedactPeds = redaction_phi(mergedClinical['AGE_AT_SEQ_REPORT'])
    logger.info("Redacting >89")
    logger.info(mergedClinical[toRedact])
    logger.info("Redacting <18")
    logger.info(mergedClinical[toRedactPeds])
    '''
    Moved to cannotReleaseHIPAA and withheld because the HIPAA
    years would change every single year.
    mergedClinical['BIRTH_YEAR'][toRedact] = "<1926"
    mergedClinical['BIRTH_YEAR'][toRedactPeds] = ">1998"
    '''
    mergedClinical['BIRTH_YEAR'][toRedact] = "cannotReleaseHIPAA"
    mergedClinical['AGE_AT_SEQ_REPORT'][toRedact] = ">32485"
    mergedClinical['BIRTH_YEAR'][toRedactPeds] = "withheld"
    mergedClinical['AGE_AT_SEQ_REPORT'][toRedactPeds] = "<6570"

    mergedClinical['INT_CONTACT'] = mergedClinical['INT_CONTACT'].fillna('')
    toRedact, toRedactPeds = redaction_phi(mergedClinical['INT_CONTACT'])

    mergedClinical['INT_CONTACT'][toRedact] = ">32485"
    mergedClinical['INT_CONTACT'][toRedactPeds] = "<6570"
    mergedClinical['BIRTH_YEAR'][toRedact] = "cannotReleaseHIPAA"
    mergedClinical['BIRTH_YEAR'][toRedactPeds] = "withheld"

    mergedClinical['INT_DOD'] = mergedClinical['INT_DOD'].fillna('')
    toRedact, toRedactPeds = redaction_phi(mergedClinical['INT_DOD'])

    mergedClinical['INT_DOD'][toRedact] = ">32485"
    mergedClinical['INT_DOD'][toRedactPeds] = "<6570"
    mergedClinical['BIRTH_YEAR'][toRedact] = "cannotReleaseHIPAA"
    mergedClinical['BIRTH_YEAR'][toRedactPeds] = "withheld"

    # Redact DFCI inputed fields
    mergedClinical['BIRTH_YEAR'] = [
        "cannotReleaseHIPAA" if ">" in str(year) else year
        for year in mergedClinical['BIRTH_YEAR']]

    mergedClinical['BIRTH_YEAR'] = [
        "withheld" if "<" in str(year) else year
        for year in mergedClinical['BIRTH_YEAR']]

    return(mergedClinical)


# Configure each maf row
def configureMafRow(rowArray, headers, keepSamples, remove_variants):
    chrom = str(rowArray[headers.index('Chromosome')])
    start = str(rowArray[headers.index('Start_Position')])
    end = str(rowArray[headers.index('End_Position')])
    ref = str(rowArray[headers.index('Reference_Allele')])
    seq = str(rowArray[headers.index('Tumor_Seq_Allele2')])
    sampleId = str(rowArray[headers.index('Tumor_Sample_Barcode')])
    variant = chrom+' '+start+' '+end+' '+ref+' '+seq+' '+sampleId
    # if pd.Series(sampleId).isin(keepSamples).any() and \
    # not pd.Series(variant).isin(remove_variants).any():
    if sampleId in keepSamples.tolist() \
            and variant not in remove_variants.tolist():
        fillnas = ['t_depth', 't_ref_count', 't_alt_count',
                   'n_depth', 'n_ref_count', 'n_alt_count']
        for i in fillnas:
            # mutationsDf[i] = mutationsDf[i].fillna("NA")
            # mutationsDf[i] = ["NA" if str(each) == "." else each
            #                    for each in  mutationsDf[i]]
            value = rowArray[headers.index(i)]
            rowArray[headers.index(i)] = "" if str(value) == "." else value

        nDepth = rowArray[headers.index("n_depth")]
        rowArray[headers.index("Match_Norm_Seq_Allele2")] = \
            '' if str(nDepth) in ["NA", "0.0"] else nDepth
        rowArray[headers.index("Match_Norm_Seq_Allele1")] = \
            '' if str(nDepth) in ["NA", "0.0"] else nDepth
        # rowArray.pop(headers.index('inBED'))
        newRow = "\t".join(rowArray)
        newRow += "\n"
        newRow = process.removeStringFloat(newRow)
        return(newRow)
    else:
        return(None)


def runMAFinBED(
        syn,
        center_mappingdf,
        test=False,
        genieVersion="test",
        genie_user=None,
        genie_pass=None):
    '''
    Run MAF in BED script, filter data and update MAFinBED database

    Args:
        syn: Synapse object
        center_mappingdf: center mapping dataframe
        test: Testing parameter. Default is False.
        genieVersion: GENIE version. Default is test.
        genie_user: Synapse username. Default is None.
        genie_pass: Synapse password.  Default is None.

    Returns:
        pd.Series: Variants to remove
    '''
    MAFinBED_script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '../analyses/genomicData/MAFinBED.R')
    notinbed_variant_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '../analyses/genomicData/notinbed.csv')

    command = ['Rscript', MAFinBED_script, notinbed_variant_file]
    if genie_user is not None and genie_pass is not None:
        command.extend(['--syn_user', genie_user, '--syn_pass', genie_pass])
    if test:
        command.append('--testing')
    subprocess.check_call(command)

    # mutationSynId = databaseSynIdMappingDf['Id'][
    #     databaseSynIdMappingDf['Database'] == "vcf2maf"][0]
    # removedVariants = syn.tableQuery(
    #     "select Chromosome, Start_Position, End_Position, Reference_Allele, "
    #     "Tumor_Seq_Allele2, Tumor_Sample_Barcode, Center from {} where inBED"
    #     " is False and Center in ('{}')".format(
    #         mutationSynId, "','".join(center_mappingdf.center)))
    # removedVariantsDf = removedVariants.asDataFrame()
    removedVariantsDf = pd.read_csv(notinbed_variant_file)
    removedVariantsDf['removeVariants'] = \
        removedVariantsDf['Chromosome'].astype(str) + ' ' + \
        removedVariantsDf['Start_Position'].astype(str) + ' ' + \
        removedVariantsDf['End_Position'].astype(str) + ' ' + \
        removedVariantsDf['Reference_Allele'].astype(str) + ' ' + \
        removedVariantsDf['Tumor_Seq_Allele2'].astype(str) + ' ' + \
        removedVariantsDf['Tumor_Sample_Barcode'].astype(str)
    # Store filtered vairants
    for center in removedVariantsDf['Center'].unique():
        center_mutation = removedVariantsDf[
            removedVariantsDf['Center'] == center]

        # mafText = process.removePandasDfFloat(center_mutation)
        center_mutation.to_csv("mafinbed_filtered_variants.csv", index=False)

        storeFile(
            syn,
            "mafinbed_filtered_variants.csv",
            parent=center_mappingdf['stagingSynId'][
                center_mappingdf['center'] == center][0],
            centerStaging=True,
            genieVersion=genieVersion)
        os.unlink("mafinbed_filtered_variants.csv")
    return(removedVariantsDf['removeVariants'])


def seq_date_filter(clinicalDf, processingDate, consortiumReleaseCutOff):
    '''
    Filter samples by seq date

    Args:
        clinicalDf: Clinical dataframe
        processingDate: Processing date in form of Apr-XXXX
        consortiumReleaseCutOff: Release cut off days

    Returns:
        list: Samples to remove
    '''
    removeSeqDateSamples = process.seqDateFilter(
        clinicalDf, processingDate, consortiumReleaseCutOff)
    return(removeSeqDateSamples)


def mutation_in_cis_filter(
        syn,
        skipMutationsInCis,
        variant_filtering_synId,
        center_mappingdf,
        genieVersion,
        test=False,
        genie_user=None,
        genie_pass=None):
    '''
    Run mutation in cis filter, look up samples to remove

    Args:
        syn: Synapse object
        skipMutationsInCis: Skip this filter
        variant_filtering_synId: mergeCheck database dataframe
        center_mappingdf: center mapping dataframe
        genieVersion: GENIE version. Default is test.
        test: Testing parameter. Default is False.
        genie_user: Synapse username. Default is None.
        genie_pass: Synapse password.  Default is None.

    Returns:
        pd.Series: Samples to remove

    '''
    if not skipMutationsInCis:
        mergeCheck_script = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            '../analyses/mergeFlag/mergeCheck.R')
        command = ['Rscript', mergeCheck_script]
        if genie_user is not None and genie_pass is not None:
            command.extend(['--syn_user', genie_user,
                            '--syn_pass', genie_pass])
        if test:
            command.append('--testing')
        subprocess.check_call(command)
        # Store each centers mutations in cis to their staging folder
        mergeCheck = syn.tableQuery(
            "select * from {} where Center in ('{}')".format(
                variant_filtering_synId,
                "','".join(center_mappingdf.center)))
        mergeCheckDf = mergeCheck.asDataFrame()
        for center in mergeCheckDf.Center.unique():
            if not pd.isnull(center):
                stagingSynId = center_mappingdf.stagingSynId[
                    center_mappingdf['center'] == center]
                mergeCheckDf[mergeCheckDf['Center'] == center].to_csv(
                    "mutationsInCis_filtered_samples.csv", index=False)
                storeFile(
                    syn,
                    "mutationsInCis_filtered_samples.csv",
                    parent=stagingSynId[0],
                    centerStaging=True,
                    genieVersion=genieVersion)
                os.unlink("mutationsInCis_filtered_samples.csv")
    variant_filtering = syn.tableQuery(
        "SELECT Tumor_Sample_Barcode FROM {} where Flag = 'TOSS' and "
        "Tumor_Sample_Barcode is not null".format(variant_filtering_synId))

    filtered_samples = variant_filtering.asDataFrame()
    # #Alex script #1 removed patients
    remove_samples = filtered_samples['Tumor_Sample_Barcode'].drop_duplicates()
    return(remove_samples)


def seq_assay_id_filter(clinicaldf):
    '''
    (Deprecated)
    Remove samples that are part of SEQ_ASSAY_IDs with less
    than 50 samples

    Args:
        clinicalDf: Sample clinical dataframe
    '''
    remove_seqassayId = clinicaldf['SEQ_ASSAY_ID'].value_counts()[
        clinicaldf['SEQ_ASSAY_ID'].value_counts() < 50]
    clinicaldf = clinicaldf[clinicaldf['SEQ_ASSAY_ID'].isin(
        remove_seqassayId.keys().values)]
    return(clinicaldf.SAMPLE_ID)


def no_genepanel_filter(clinicaldf, beddf):
    '''
    Remove samples that don't have bed files associated with
    them

    Args:
        clinicaldf:  Clinical data
        beddf: bed data

    Returns:
        series: samples to remove
    '''

    logger.info("NO GENE PANEL FILTER")
    has_seq_assay = clinicaldf['SEQ_ASSAY_ID'].isin(beddf['SEQ_ASSAY_ID'])
    remove_samples = clinicaldf['SAMPLE_ID'][~has_seq_assay]
    logger.info("Removing samples with no bed file: {}".format(
        ",".join(remove_samples)))
    return(remove_samples)


def store_gene_panel_files(
        syn,
        fileviewSynId,
        genieVersion,
        data_gene_panel,
        consortiumReleaseSynId,
        current_release_staging):
    # Only need to upload these files once
    logger.info("STORING GENE PANELS FILES")
    genePanels = syn.tableQuery(
        "select id from %s where cBioFileFormat = 'genePanel' "
        "and fileStage = 'staging'" % fileviewSynId)
    genePanelDf = genePanels.asDataFrame()
    genePanelEntities = []
    for synId in genePanelDf['id']:
        genePanel = syn.get(synId)
        panelNames = set(data_gene_panel['mutations'])
        genePanelName = os.path.basename(genePanel.path)
        newGenePanelPath = os.path.join(
            GENIE_RELEASE_DIR,
            genePanelName.replace(".txt", "_%s.txt" % genieVersion))
        if genePanelName.replace(".txt", "").replace(
                "data_gene_panel_", "") in panelNames:
            os.rename(genePanel.path, newGenePanelPath)
            genePanelEntities.append(storeFile(
                syn, newGenePanelPath,
                parent=consortiumReleaseSynId,
                genieVersion=genieVersion,
                name=genePanelName,
                cBioFileFormat="genePanel",
                staging=current_release_staging))
    return(genePanelEntities)


def store_fusion_files(
        syn,
        consortiumReleaseSynId,
        genieVersion,
        fusionSynId,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        current_release_staging,
        CENTER_MAPPING_DF):
    # FUSIONS
    logger.info("MERING, FILTERING, STORING FUSION FILES")
    Fusions = syn.tableQuery(
        'SELECT HUGO_SYMBOL,ENTREZ_GENE_ID,CENTER,TUMOR_SAMPLE_BARCODE,FUSION,'
        'DNA_SUPPORT,RNA_SUPPORT,METHOD,FRAME FROM %s' % fusionSynId)
    FusionsDf = Fusions.asDataFrame()
    FusionsDf['ENTREZ_GENE_ID'][FusionsDf['ENTREZ_GENE_ID'] == 0] = pd.np.nan

    if not current_release_staging:
        FusionsStagingDf = FusionsDf[FusionsDf['TUMOR_SAMPLE_BARCODE'].isin(
            keepForCenterConsortiumSamples)]
        for center in CENTER_MAPPING_DF.center:
            center_fusion = FusionsStagingDf[
                FusionsStagingDf['CENTER'] == center]
            if not center_fusion.empty:
                center_fusion.to_csv(
                    FUSIONS_CENTER_PATH % center,
                    sep="\t", index=False)
                storeFile(
                    syn, FUSIONS_CENTER_PATH % center,
                    genieVersion=genieVersion,
                    parent=CENTER_MAPPING_DF['stagingSynId'][
                        CENTER_MAPPING_DF['center'] == center][0],
                    centerStaging=True)

    FusionsDf = FusionsDf[FusionsDf['TUMOR_SAMPLE_BARCODE'].isin(
        keepForMergedConsortiumSamples)]
    FusionsDf = FusionsDf.rename(columns={
        'HUGO_SYMBOL': 'Hugo_Symbol',
        'ENTREZ_GENE_ID': 'Entrez_Gene_Id',
        'CENTER': 'Center',
        'TUMOR_SAMPLE_BARCODE': 'Tumor_Sample_Barcode',
        'FUSION': 'Fusion',
        'DNA_SUPPORT': 'DNA_support',
        'RNA_SUPPORT': 'RNA_support',
        'METHOD': 'Method',
        'FRAME': 'Frame'})
    # Remove duplicated Fusions
    FusionsDf = FusionsDf[~FusionsDf[
        ['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Fusion']].duplicated()]
    # FusionsDf.to_csv(FUSIONS_PATH, sep="\t", index=False)
    fusionText = process.removePandasDfFloat(FusionsDf)
    fusions_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_fusions_%s.txt' % genieVersion)
    with open(fusions_path, "w") as fusionFile:
        fusionFile.write(fusionText)
    storeFile(
        syn, fusions_path,
        parent=consortiumReleaseSynId,
        genieVersion=genieVersion,
        name="data_fusions.txt",
        staging=current_release_staging)


def store_maf_files(
        syn,
        genieVersion,
        centerMafFileViewSynId,
        consortiumReleaseSynId,
        clinicalDf,
        CENTER_MAPPING_DF,
        keepForMergedConsortiumSamples,
        keepForCenterConsortiumSamples,
        remove_mafInBed_variants,
        current_release_staging):
    logger.info("FILTERING, STORING MUTATION FILES")
    centerMafSynIds = syn.tableQuery(
        "select id from %s " % centerMafFileViewSynId +
        "where name like '%mutation%'")
    centerMafSynIdsDf = centerMafSynIds.asDataFrame()
    mutations_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_mutations_extended_%s.txt' % genieVersion)
    sequenced_samples = "#sequenced_samples: {}".format(
        " ".join(clinicalDf['SAMPLE_ID']))
    with open(mutations_path, 'w') as f:
        f.write(sequenced_samples + "\n")
    for index, mafSynId in enumerate(centerMafSynIdsDf.id):
        mafEnt = syn.get(mafSynId)
        logger.info(mafEnt.path)
        with open(mafEnt.path, "r") as mafFile:
            header = mafFile.readline()
            headers = header.replace("\n", "").split("\t")
            if index == 0:
                with open(mutations_path, 'a') as f:
                    f.write(header)
                # Create maf file per center for their staging directory
                for center in clinicalDf['CENTER'].unique():
                    with open(MUTATIONS_CENTER_PATH % center, 'w') as f:
                        f.write(header)
        # with open(mafEnt.path,"r") as newMafFile:
        #   newMafFile.readline()
            center = mafEnt.path.split("_")[3]
            # Make sure to only write the centers that release = True
            if center in CENTER_MAPPING_DF.center.tolist():
                for row in mafFile:
                    rowArray = row.replace("\n", "").split("\t")
                    center = rowArray[headers.index('Center')]
                    newMergedRow = configureMafRow(
                        rowArray, headers,
                        keepForMergedConsortiumSamples,
                        remove_mafInBed_variants)
                    if newMergedRow is not None:
                        with open(mutations_path, 'a') as f:
                            f.write(newMergedRow)
                    newCenterRow = configureMafRow(
                        rowArray, headers,
                        keepForCenterConsortiumSamples,
                        remove_mafInBed_variants)
                    if newCenterRow is not None:
                        with open(MUTATIONS_CENTER_PATH % center, 'a') as f:
                            f.write(newCenterRow)
    storeFile(
        syn, mutations_path,
        parent=consortiumReleaseSynId,
        genieVersion=genieVersion,
        name="data_mutations_extended.txt",
        staging=current_release_staging)

    if not current_release_staging:
        for center in clinicalDf['CENTER'].unique():
            storeFile(
                syn, MUTATIONS_CENTER_PATH % center,
                genieVersion=genieVersion,
                parent=CENTER_MAPPING_DF['stagingSynId'][
                    CENTER_MAPPING_DF['center'] == center][0],
                centerStaging=True)


def run_genie_filters(
        syn,
        genie_user,
        genie_pass,
        genieVersion,
        variant_filtering_synId,
        clinicalDf,
        bedDf,
        sampleCols,
        CENTER_MAPPING_DF,
        processingDate,
        skipMutationsInCis,
        consortiumReleaseCutOff,
        test):
    # ADD CHECKS TO CODE BEFORE UPLOAD.
    # Throw error if things don't go through
    logger.info("RUN GENIE FILTERS")
    # STORING CLINICAL FILES INTO CBIOPORTAL

    ''' FILTERING '''
    logger.info("MAF IN BED FILTER")
    remove_mafInBed_variants = runMAFinBED(
        syn, CENTER_MAPPING_DF, test=test,
        genieVersion=genieVersion,
        genie_user=genie_user, genie_pass=genie_pass)

    logger.info("MUTATION IN CIS FILTER")
    remove_mutationInCis_samples = mutation_in_cis_filter(
        syn, skipMutationsInCis, variant_filtering_synId, CENTER_MAPPING_DF,
        genieVersion=genieVersion, test=test,
        genie_user=genie_user, genie_pass=genie_pass)

    remove_no_genepanel_samples = no_genepanel_filter(clinicalDf, bedDf)

    logger.info("SEQ DATE FILTER")
    remove_seqDate_samples = seq_date_filter(
        clinicalDf, processingDate, consortiumReleaseCutOff)
    # Only certain samples are removed for the files that go into
    # staging center folder
    removeForCenterConsortiumSamples = set(remove_mutationInCis_samples).union(
        set(remove_no_genepanel_samples))
    # Most filteres are applied for the files that go into the merged
    # consortium release
    removeForMergedConsortiumSamples = set(remove_seqDate_samples)
    # set(remove_seqAssayId_samples)#.union(set(remove_seqDate_samples))
    removeForMergedConsortiumSamples = \
        removeForMergedConsortiumSamples.union(
            removeForCenterConsortiumSamples)

    return(remove_mafInBed_variants,
           removeForMergedConsortiumSamples,
           removeForCenterConsortiumSamples)


def store_clinical_files(
        syn,
        genieVersion,
        clinicalDf,
        oncotree_url,
        sampleCols,
        patientCols,
        removeForCenterConsortiumSamples,
        removeForMergedConsortiumSamples,
        consortiumReleaseSynId,
        current_release_staging,
        CENTER_MAPPING_DF):
    logger.info("CONFIGURING CLINICAL FILES")
    logger.info("REMOVING PHI")
    clinicalDf = reAnnotatePHI(clinicalDf)
    logger.info("ADD CANCER TYPES")
    # This removes support for both oncotree urls (only support json)
    oncotreeDict = process.get_oncotree_code_mappings(oncotree_url)
    # Add in unknown key which maps to UNKNOWN everything
    oncotreeDict['UNKNOWN'] = {
        'CANCER_TYPE': 'UNKNOWN',
        'CANCER_TYPE_DETAILED': 'UNKNOWN',
        'ONCOTREE_PRIMARY_NODE': 'UNKNOWN',
        'ONCOTREE_SECONDARY_NODE': 'UNKNOWN'}

    clinicalDf['CANCER_TYPE'] = [
        oncotreeDict[code.upper()]["CANCER_TYPE"]
        if code.upper() in oncotreeDict.keys() else float('nan')
        for code in clinicalDf['ONCOTREE_CODE']]

    clinicalDf['CANCER_TYPE_DETAILED'] = [
        oncotreeDict[code.upper()]["CANCER_TYPE_DETAILED"]
        if code.upper() in oncotreeDict.keys() else float('nan')
        for code in clinicalDf['ONCOTREE_CODE']]

    clinicalDf['ONCOTREE_PRIMARY_NODE'] = [
        oncotreeDict[code.upper()]["ONCOTREE_PRIMARY_NODE"]
        if code.upper() in oncotreeDict.keys() else float('nan')
        for code in clinicalDf['ONCOTREE_CODE']]

    clinicalDf['ONCOTREE_SECONDARY_NODE'] = [
        oncotreeDict[code.upper()]["ONCOTREE_SECONDARY_NODE"]
        if code.upper() in oncotreeDict.keys() else float('nan')
        for code in clinicalDf['ONCOTREE_CODE']]

    # All cancer types that are null should have null oncotree codes
    clinicalDf['ONCOTREE_CODE'][
        clinicalDf['CANCER_TYPE'].isnull()] = float('nan')
    # Suggest using AGE_AT_SEQ_REPORT_DAYS instead so that the
    # descriptions can match
    clinicalDf['AGE_AT_SEQ_REPORT_DAYS'] = clinicalDf['AGE_AT_SEQ_REPORT']
    clinicalDf['AGE_AT_SEQ_REPORT'] = [
        int(math.floor(int(float(i))/365.25))
        if process.checkInt(i) else i
        for i in clinicalDf['AGE_AT_SEQ_REPORT']]
    clinicalDf['AGE_AT_SEQ_REPORT'][
        clinicalDf['AGE_AT_SEQ_REPORT'] == ">32485"] = ">89"
    clinicalDf['AGE_AT_SEQ_REPORT'][
        clinicalDf['AGE_AT_SEQ_REPORT'] == "<6570"] = "<18"

    ############################################################
    # CENTER SPECIFIC CODE FOR RIGHT NOW (REMOVE UHN-555-V1)
    ############################################################
    clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'] != "UHN-555-V1"]
    # clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'] != "PHS-TRISEQ-V1"]
    # clinicalDf = clinicalDf[clinicalDf['CENTER'] != "WAKE"]
    # clinicalDf = clinicalDf[clinicalDf['CENTER'] != "CRUK"]
    ############################################################
    ############################################################

    clinicalDf.drop_duplicates("SAMPLE_ID", inplace=True)

    logger.info("STORING CLINICAL FILES")
    # samples must be removed after reading in the clinical file again
    clinicalDfStaging = clinicalDf[
        ~clinicalDf['SAMPLE_ID'].isin(removeForCenterConsortiumSamples)]
    if not current_release_staging:
        for center in CENTER_MAPPING_DF.center:
            center_clinical = clinicalDfStaging[
                clinicalDfStaging['CENTER'] == center]
            center_sample = center_clinical[
                sampleCols].drop_duplicates('SAMPLE_ID')
            center_patient = center_clinical[
                patientCols].drop_duplicates('PATIENT_ID')
            center_sample.to_csv(
                SAMPLE_CENTER_PATH % center, sep="\t", index=False)
            center_patient.to_csv(
                PATIENT_CENTER_PATH % center, sep="\t", index=False)
            storeFile(
                syn, SAMPLE_CENTER_PATH % center,
                genieVersion=genieVersion,
                parent=CENTER_MAPPING_DF['stagingSynId'][
                    CENTER_MAPPING_DF['center'] == center][0],
                centerStaging=True)
            storeFile(
                syn, PATIENT_CENTER_PATH % center,
                genieVersion=genieVersion,
                parent=CENTER_MAPPING_DF['stagingSynId'][
                    CENTER_MAPPING_DF['center'] == center][0],
                centerStaging=True)

    clinicalDf = clinicalDf[~clinicalDf['SAMPLE_ID'].isin(
        removeForMergedConsortiumSamples)]
    # This must happen here because the seq assay filter
    # must happen after all the other filters
    # logger.info("SEQ ASSAY FILTER")
    # remove_seqAssayId_samples = seq_assay_id_filter(clinicalDf)
    # removeForMergedConsortiumSamples = \
    #     removeForMergedConsortiumSamples.union(set(remove_seqAssayId_samples))
    # clinicalDf = clinicalDf[~clinicalDf['SAMPLE_ID'].isin(
    #    remove_seqAssayId_samples)]

    keepForCenterConsortiumSamples = clinicalDfStaging.SAMPLE_ID
    keepForMergedConsortiumSamples = clinicalDf.SAMPLE_ID
    # This mapping table is the GENIE clinical code to description
    # mapping to generate the headers of the clinical file
    mapping_table = syn.tableQuery('SELECT * FROM syn9621600')
    mapping = mapping_table.asDataFrame()
    clinical_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_clinical_%s.txt' % genieVersion)
    clinical_sample_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_clinical_sample_%s.txt' % genieVersion)
    clinical_patient_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_clinical_patient_%s.txt' % genieVersion)
    process.addClinicalHeaders(
        clinicalDf, mapping, patientCols, sampleCols,
        clinical_sample_path, clinical_patient_path)
    storeFile(
        syn, clinical_sample_path,
        parent=consortiumReleaseSynId,
        genieVersion=genieVersion,
        name="data_clinical_sample.txt",
        staging=current_release_staging)

    storeFile(
        syn, clinical_patient_path,
        parent=consortiumReleaseSynId,
        genieVersion=genieVersion,
        name="data_clinical_patient.txt",
        staging=current_release_staging)

    clinicalDf.to_csv(clinical_path, sep="\t", index=False)
    storeFile(
        syn, clinical_path,
        parent=consortiumReleaseSynId,
        name="data_clinical.txt",
        staging=current_release_staging)

    return(clinicalDf,
           keepForCenterConsortiumSamples,
           keepForMergedConsortiumSamples)


def store_cna_files(
        syn,
        centerMafFileViewSynId,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        CENTER_MAPPING_DF,
        genieVersion,
        consortiumReleaseSynId,
        current_release_staging):
    # CNA
    logger.info("MERING, FILTERING, STORING CNA FILES")
    cna_path = os.path.join(
        GENIE_RELEASE_DIR, "data_CNA_%s.txt" % genieVersion)
    centerCNASynIds = syn.tableQuery(
        "select id from %s " % centerMafFileViewSynId +
        "where name like 'data_CNA%'")
    centerCNASynIdsDf = centerCNASynIds.asDataFrame()
    # Grab all unique symbols and form cnaTemplate
    allSymbols = set()

    for cnaSynId in centerCNASynIdsDf.id:
        cnaEnt = syn.get(cnaSynId)
        with open(cnaEnt.path, "r") as cnaFile:
            # Read first line first
            cnaFile.readline()
            # Get all hugo symbols
            allSymbols = allSymbols.union(
                set(line.split("\t")[0] for line in cnaFile))
    cnaTemplate = pd.DataFrame({"Hugo_Symbol": list(allSymbols)})
    cnaTemplate.sort_values("Hugo_Symbol", inplace=True)
    cnaTemplate.to_csv(cna_path, sep="\t", index=False)
    # Loop through to create finalized CNA file
    withCenterHugoSymbol = pd.Series("Hugo_Symbol")
    withCenterHugoSymbol = withCenterHugoSymbol.append(
        pd.Series(keepForCenterConsortiumSamples))

    withMergedHugoSymbol = pd.Series("Hugo_Symbol")
    withMergedHugoSymbol = withMergedHugoSymbol.append(
        pd.Series(keepForMergedConsortiumSamples))

    cnaSamples = []

    for cnaSynId in centerCNASynIdsDf.id:
        cnaEnt = syn.get(cnaSynId)
        center = cnaEnt.name.replace("data_CNA_", "").replace(".txt", "")
        logger.info(cnaEnt.path)
        if center in CENTER_MAPPING_DF.center.tolist():
            centerCNA = pd.read_csv(cnaEnt.path, sep="\t")
            merged = cnaTemplate.merge(
                centerCNA, on="Hugo_Symbol", how="outer")
            merged.sort_values("Hugo_Symbol", inplace=True)

            merged = merged[merged.columns[
                merged.columns.isin(withCenterHugoSymbol)]]

            cnaText = process.removePandasDfFloat(merged)
            # Replace blank with NA's
            cnaText = cnaText.replace("\t\t", "\tNA\t").replace(
                "\t\t", "\tNA\t").replace('\t\n', "\tNA\n")

            # Store center CNA file in staging dir
            with open(CNA_CENTER_PATH % center, "w") as cnaFile:
                cnaFile.write(cnaText)
            storeFile(
                syn, CNA_CENTER_PATH % center,
                genieVersion=genieVersion,
                parent=CENTER_MAPPING_DF['stagingSynId'][
                    CENTER_MAPPING_DF['center'] == center][0],
                centerStaging=True)
            # This is to remove more samples for the final cna file
            merged = merged[merged.columns[
                merged.columns.isin(withMergedHugoSymbol)]]

            cnaText = process.removePandasDfFloat(merged)
            cnaText = cnaText.replace("\t\t", "\tNA\t").replace(
                "\t\t", "\tNA\t").replace('\t\n', "\tNA\n")

            with open(CNA_CENTER_PATH % center, "w") as cnaFile:
                cnaFile.write(cnaText)
            # Join CNA file
            cnaSamples.extend(merged.columns[1:].tolist())
            joinCommand = ["join", cna_path, CNA_CENTER_PATH % center]
            output = subprocess.check_output(joinCommand)
            with open(cna_path, "w") as cnaFile:
                cnaFile.write(output.decode("utf-8").replace(" ", "\t"))

    storeFile(
        syn, cna_path,
        parent=consortiumReleaseSynId,
        genieVersion=genieVersion,
        name="data_CNA.txt",
        staging=current_release_staging)

    return(cnaSamples)


# SEG
def store_seg_files(
        syn,
        genie_version,
        seg_synid,
        release_synid,
        keep_for_center_consortium_samples,
        keep_for_merged_consortium_samples,
        center_mappingdf,
        current_release_staging):
    '''
    Create, filter and store seg file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        seg_synid: Seg database synid
        release_synid: Synapse id to store release file
        keep_for_center_consortium_samples: Samples to keep for center files
        keep_for_merged_consortium_samples: Samples to keep for merged file
        center_mappingdf: Center mapping dataframe
        current_release_staging: Staging flag
    '''
    logger.info("MERING, FILTERING, STORING SEG FILES")
    seg_path = os.path.join(
        GENIE_RELEASE_DIR,
        'genie_private_data_cna_hg19_%s.seg' % genie_version)
    seg = syn.tableQuery(
        'SELECT ID,CHROM,LOCSTART,LOCEND,NUMMARK,SEGMEAN'
        ',CENTER FROM %s' % seg_synid)
    segdf = seg.asDataFrame()
    segdf = segdf.rename(columns={
        'CHROM': 'chrom',
        'LOCSTART': 'loc.start',
        'LOCEND': 'loc.end',
        'SEGMEAN': 'seg.mean',
        'NUMMARK': 'num.mark'})
    if not current_release_staging:
        segStagingDf = segdf[segdf['ID'].isin(
            keep_for_center_consortium_samples)]
        for center in center_mappingdf.center:
            center_seg = segStagingDf[segStagingDf['CENTER'] == center]
            if not center_seg.empty:
                del center_seg['CENTER']
                segtext = process.removePandasDfFloat(center_seg)
                with open(SEG_CENTER_PATH % center, "w") as segfile:
                    segfile.write(segtext)
                storeFile(
                    syn, SEG_CENTER_PATH % center,
                    genieVersion=genie_version,
                    parent=center_mappingdf['stagingSynId'][
                        center_mappingdf['center'] == center][0],
                    centerStaging=True)
    del segdf['CENTER']
    segdf = segdf[segdf['ID'].isin(keep_for_merged_consortium_samples)]
    segtext = process.removePandasDfFloat(segdf)
    with open(seg_path, "w") as segfile:
        segfile.write(segtext)
    storeFile(
        syn, seg_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="genie_private_data_cna_hg19.seg",
        staging=current_release_staging)


def store_data_gene_matrix(
        syn,
        genie_version,
        clinicaldf,
        cna_samples,
        release_synid,
        current_release_staging):
    '''
    Create and store data gene matrix file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        clinicaldf: Clinical dataframe with SAMPLE_ID and SEQ_ASSAY_ID
        cna_samples: Samples with CNA
        release_synid: Synapse id to store release file
        current_release_staging: Staging flag

    Returns:
        dataframe: data gene matrix dataframe
    '''
    logger.info("STORING DATA GENE MATRIX FILE")
    data_gene_matrix_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_gene_matrix_%s.txt' % genie_version)
    # Samples have already been removed
    data_gene_matrix = pd.DataFrame(columns=["SAMPLE_ID", "SEQ_ASSAY_ID"])
    data_gene_matrix = data_gene_matrix.append(
        clinicaldf[['SAMPLE_ID', 'SEQ_ASSAY_ID']])
    data_gene_matrix = data_gene_matrix.rename(
        columns={"SEQ_ASSAY_ID": "mutations"})
    data_gene_matrix = data_gene_matrix[data_gene_matrix['SAMPLE_ID'] != ""]
    data_gene_matrix.drop_duplicates("SAMPLE_ID", inplace=True)
    # Gene panel file is written below CNA, because of the "cna" column
    # Add in CNA column into gene panel file
    cna_seqids = data_gene_matrix['mutations'][
        data_gene_matrix['SAMPLE_ID'].isin(cna_samples)].unique()
    data_gene_matrix['cna'] = data_gene_matrix['mutations']
    data_gene_matrix['cna'][~data_gene_matrix['cna'].isin(cna_seqids)] = "NA"
    data_gene_matrix.to_csv(data_gene_matrix_path, sep="\t", index=False)

    storeFile(
        syn, data_gene_matrix_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="data_gene_matrix.txt",
        staging=current_release_staging)
    return(data_gene_matrix)


def store_bed_files(
        syn,
        genie_version,
        beddf,
        seq_assay_ids,
        center_mappingdf,
        current_release_staging,
        release_synid):
    '''
    Store bed files, store the bed regions that had symbols remapped
    Filters bed file by clinical dataframe seq assays

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        beddf: Bed dataframe
        seq_assay_ids: All SEQ_ASSAY_IDs in the clinical file
        center_mappingdf: Center mapping dataframe
        current_release_staging: Staging flag
        release_synid: Synapse id to store release file
    '''
    logger.info("STORING COMBINED BED FILE")
    combined_bed_path = os.path.join(
        GENIE_RELEASE_DIR, 'genie_combined_%s.bed' % genie_version)
    if not current_release_staging:
        for seq_assay in beddf['SEQ_ASSAY_ID'].unique():
            bed_seq_df = beddf[beddf['SEQ_ASSAY_ID'] == seq_assay]
            center = seq_assay.split("-")[0]
            bed_seq_df = bed_seq_df[bed_seq_df['Hugo_Symbol'] != bed_seq_df['ID']]
            if not bed_seq_df.empty:
                bed_seq_df.to_csv(
                    BED_DIFFS_SEQASSAY_PATH % seq_assay,
                    index=False)
                storeFile(
                    syn, BED_DIFFS_SEQASSAY_PATH % seq_assay,
                    genieVersion=genie_version,
                    parent=center_mappingdf['stagingSynId'][
                        center_mappingdf['center'] == center][0],
                    centerStaging=True)
    # This clinicalDf is already filtered through most of the filters
    beddf = beddf[beddf['SEQ_ASSAY_ID'].isin(seq_assay_ids)]
    beddf.to_csv(combined_bed_path, sep="\t", index=False)
    storeFile(
        syn, combined_bed_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="genie_combined.bed",
        staging=current_release_staging)


def stagingToCbio(
        syn, processingDate, genieVersion,
        CENTER_MAPPING_DF, databaseSynIdMappingDf,
        oncotree_url=None, consortiumReleaseCutOff=183,
        current_release_staging=False,
        skipMutationsInCis=False, test=False,
        genie_user=None, genie_pass=None):
    '''
    Main function that takes the GENIE database and creates release files

    Args:
        syn: Synapse object
        processingDate: Processing date in form of Apr-XXXX
        genieVersion: GENIE version. Default is test.
        CENTER_MAPPING_DF: center mapping dataframe
        databaseSynIdMappingDf: Database to Synapse Id mapping
        oncotree_url: Oncotree link
        consortiumReleaseCutOff: Release cut off days
        current_release_staging: Is it staging. Default is False.
        skipMutationsInCis: Skip mutation in cis filter. Default is False.
        test: Testing parameter. Default is False.
        genie_user: Synapse username. Default is None.
        genie_pass: Synapse password.  Default is None.

    '''
    if not os.path.exists(GENIE_RELEASE_DIR):
        os.mkdir(GENIE_RELEASE_DIR)
    consortiumReleaseSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "consortium"][0]
    centerMafFileViewSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "centerMafView"][0]
    patientSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "patient"][0]
    sampleSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "sample"][0]
    bedSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "bed"][0]
    fileviewSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "fileview"][0]
    segSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "seg"][0]
    variant_filtering_synId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "mutationsInCis"][0]
    fusionSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == "fusions"][0]
    # Using center mapping df to gate centers in release fileStage
    patient = syn.tableQuery(
        "SELECT * FROM {} where CENTER in ('{}')".format(
            patientSynId, "','".join(CENTER_MAPPING_DF.center)))
    sample = syn.tableQuery(
        "SELECT * FROM {} where CENTER in ('{}')".format(
            sampleSynId, "','".join(CENTER_MAPPING_DF.center)))
    bed = syn.tableQuery(
        "SELECT Chromosome,Start_Position,End_Position,Hugo_Symbol,ID,"
        "SEQ_ASSAY_ID,Feature_Type,includeInPanel FROM"
        " {} where CENTER in ('{}')".format(
            bedSynId, "','".join(CENTER_MAPPING_DF.center)))
    patientDf = patient.asDataFrame()
    sampleDf = sample.asDataFrame()
    bedDf = bed.asDataFrame()

    # Clinical release scope filter
    # If private -> Don't release to public
    clinicalReleaseScope = syn.tableQuery(
        "SELECT * FROM syn8545211 where releaseScope <> 'private'")
    clinicalReleaseScopeDf = clinicalReleaseScope.asDataFrame()

    patientCols = clinicalReleaseScopeDf['fieldName'][
        clinicalReleaseScopeDf['level'] == "patient"].tolist()
    sampleCols = clinicalReleaseScopeDf['fieldName'][
        clinicalReleaseScopeDf['level'] == "sample"].tolist()

    # Remove this when these columns are removed from both databases
    if sampleDf.get("AGE_AT_SEQ_REPORT_NUMERICAL") is not None:
        del sampleDf['AGE_AT_SEQ_REPORT_NUMERICAL']
    del sampleDf['CENTER']
    # Remove this when these columns are removed from both databases
    if patientDf.get("BIRTH_YEAR_NUMERICAL") is not None:
        del patientDf['BIRTH_YEAR_NUMERICAL']
    # del patientDf['BIRTH_YEAR_NUMERICAL']

    totalSample = ['PATIENT_ID']
    totalSample.extend(sampleCols)
    sampleCols = totalSample
    # Make sure to only grab samples that have patient information
    sampleDf = sampleDf[sampleDf['PATIENT_ID'].isin(patientDf['PATIENT_ID'])]
    clinicalDf = sampleDf.merge(patientDf, on="PATIENT_ID", how="outer")
    # Remove patients without any sample or patient ids
    clinicalDf = clinicalDf[~clinicalDf['SAMPLE_ID'].isnull()]
    clinicalDf = clinicalDf[~clinicalDf['PATIENT_ID'].isnull()]

    remove_mafInBed_variants, \
        removeForMergedConsortiumSamples, \
        removeForCenterConsortiumSamples = \
        run_genie_filters(
            syn,
            genie_user,
            genie_pass,
            genieVersion,
            variant_filtering_synId,
            clinicalDf,
            bedDf,
            sampleCols,
            CENTER_MAPPING_DF,
            processingDate,
            skipMutationsInCis,
            consortiumReleaseCutOff,
            test)

    clinicalDf, keepForCenterConsortiumSamples, keepForMergedConsortiumSamples\
        = store_clinical_files(
            syn,
            genieVersion,
            clinicalDf,
            oncotree_url,
            sampleCols,
            patientCols,
            removeForCenterConsortiumSamples,
            removeForMergedConsortiumSamples,
            consortiumReleaseSynId,
            current_release_staging,
            CENTER_MAPPING_DF)

    assert not clinicalDf['SAMPLE_ID'].duplicated().any()

    store_maf_files(
        syn,
        genieVersion,
        centerMafFileViewSynId,
        consortiumReleaseSynId,
        clinicalDf[['SAMPLE_ID', 'CENTER']],
        CENTER_MAPPING_DF,
        keepForMergedConsortiumSamples,
        keepForCenterConsortiumSamples,
        remove_mafInBed_variants,
        current_release_staging)

    cnaSamples = store_cna_files(
        syn,
        centerMafFileViewSynId,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        CENTER_MAPPING_DF,
        genieVersion,
        consortiumReleaseSynId,
        current_release_staging)

    data_gene_matrix = store_data_gene_matrix(
        syn,
        genieVersion,
        clinicalDf,
        cnaSamples,
        consortiumReleaseSynId,
        current_release_staging)

    genePanelEntities = store_gene_panel_files(
        syn,
        fileviewSynId,
        genieVersion,
        data_gene_matrix,
        consortiumReleaseSynId,
        current_release_staging)

    store_fusion_files(
        syn,
        consortiumReleaseSynId,
        genieVersion,
        fusionSynId,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        current_release_staging,
        CENTER_MAPPING_DF)

    store_seg_files(
        syn,
        genieVersion,
        segSynId,
        consortiumReleaseSynId,
        keepForCenterConsortiumSamples,
        keepForMergedConsortiumSamples,
        CENTER_MAPPING_DF,
        current_release_staging)

    store_bed_files(
        syn,
        genieVersion,
        bedDf,
        clinicalDf['SEQ_ASSAY_ID'].unique(),
        CENTER_MAPPING_DF,
        current_release_staging,
        consortiumReleaseSynId)

    return(genePanelEntities)


def update_process_trackingdf(
        syn, process_trackerdb_synid, center, process_type, start=True):
    '''
    Updates the processing tracking database

    Args:
        syn: synapse object
        process_trackerdb_synid: Synapse id of process tracking db
        center: GENIE center (ie. SAGE)
        process_type: processing type (ie. dbToStage)
        start: Start or end of processing.  Default is True for start
    '''
    logger.info("UPDATE PROCESS TRACKING TABLE")
    column = 'timeStartProcessing' if start else 'timeEndProcessing'
    process_tracker = syn.tableQuery(
        "SELECT {col} FROM {tableid} where center = '{center}' "
        "and processingType = '{process_type}'".format(
            col=column,
            tableid=process_trackerdb_synid,
            center=center,
            process_type=process_type))
    process_trackerdf = process_tracker.asDataFrame()
    process_trackerdf[column][0] = str(int(time.time()*1000))
    syn.store(synapseclient.Table(process_trackerdb_synid, process_trackerdf))


# def command_reviseMetadataFiles(syn, args, databaseSynIdMappingDf):
#     '''
#     Command to call metadata files with args

#     Args:
#         syn: Synapse object
#         args: Argument list
#         databaseSynIdMappingDf: database to synapse id mapping df
#     '''
#     reviseMetadataFiles(
#         syn, args.staging, databaseSynIdMappingDf, args.genieVersion)


def reviseMetadataFiles(
        syn, staging, consortiumId, genieVersion=None):
    '''
    Rewrite metadata files with the correct GENIE version

    Args:
        syn: Synapse object
        staging: staging flag
        databaseSynIdMappingDf: database to synapse id mapping df
        genieVersion: GENIE version, Default to None
    '''
    allFiles = syn.getChildren(consortiumId)
    metadataEnts = [
        syn.get(i['id'],
                downloadLocation=GENIE_RELEASE_DIR,
                ifcollision="overwrite.local")
        for i in allFiles if 'meta' in i['name']]

    for metaEnt in metadataEnts:
        with open(metaEnt.path, "r+") as meta:
            metaText = meta.read()
            if "meta_study" not in metaEnt.path:
                version = ''
            else:
                version = re.search(".+GENIE.+v(.+)", metaText).group(1)
            # Fix this line
            genieVersion = version if genieVersion is None else genieVersion
            dataFileVersion = re.search(".+data_(.+)[.]txt", metaText)
            if dataFileVersion is None:
                dataFileVersion = re.search(".+data_(.+)[.]seg", metaText)
            if dataFileVersion is not None:
                dataFileVersion = dataFileVersion.group(1).split("_")[-1]

            if version != genieVersion:
                metaText = metaText.replace(
                    "GENIE Cohort v{}".format(version),
                    "GENIE Cohort v{}".format(genieVersion))

                metaText = metaText.replace(
                    "GENIE v{}".format(version),
                    "GENIE v{}".format(genieVersion))

                if dataFileVersion is not None:
                    metaText = metaText.replace(dataFileVersion, genieVersion)
                    metaText = metaText.replace(dataFileVersion, genieVersion)
                meta.seek(0)
                meta.write(metaText)
                meta.truncate()
                storeFile(
                    syn,
                    metaEnt.path,
                    parent=consortiumId,
                    genieVersion=genieVersion,
                    staging=staging)


def createLinkVersion(
        syn,
        genieVersion,
        caseListEntities,
        genePanelEntities,
        databaseSynIdMappingDf):
    versioning = genieVersion.split(".")
    main = versioning[0]
    # second = ".".join(versioning[1:])
    releaseSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'release'].values[0]
    releases = synu.walk(syn, releaseSynId)
    mainReleaseFolders = next(releases)[1]
    releaseFolderSynId = [
        synId for folderName, synId in mainReleaseFolders
        if folderName == "Release {}".format(main)]

    if len(releaseFolderSynId) > 0:
        secondRelease = synu.walk(syn, releaseFolderSynId[0])
        secondReleaseFolders = next(secondRelease)[1]
        secondReleaseFolderSynIdList = [
            synId for folderName, synId in secondReleaseFolders
            if folderName == genieVersion]

        if len(secondReleaseFolderSynIdList) > 0:
            secondReleaseFolderSynId = secondReleaseFolderSynIdList[0]
        else:
            secondReleaseFolderSynId = syn.store(
                synapseclient.Folder(genieVersion,
                                     parent=releaseFolderSynId[0]))
    else:
        mainReleaseFolderId = syn.store(
            synapseclient.Folder("Release %s" % main,
                                 parent=releaseSynId)).id
        secondReleaseFolderSynId = syn.store(
            synapseclient.Folder(genieVersion,
                                 parent=mainReleaseFolderId)).id

    caselistId = findCaseListId(syn, secondReleaseFolderSynId)
    consortiumSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'consortium'].values[0]
    consortiumRelease = syn.getChildren(consortiumSynId)
    # data_clinical.txt MUST be pulled in because the clinical file is
    # needed in the consortium_to_public.py
    for ents in consortiumRelease:
        if ents['type'] != "org.sagebionetworks.repo.model.Folder" \
                and not ents['name'].startswith("data_gene_panel"):
            syn.store(synapseclient.Link(
                ents['id'],
                parent=secondReleaseFolderSynId,
                targetVersion=ents['versionNumber']))

    releaseFiles = syn.getChildren(secondReleaseFolderSynId)
    clinEnt = [
        ents['id']
        for ents in releaseFiles
        if ents['name'] == "data_clinical.txt"][0]
    # Set private permission for the data_clinical.txt Link
    syn.setPermissions(clinEnt, principalId=3346558, accessType=[])
    syn.setPermissions(clinEnt, principalId=3326313, accessType=[])
    # caselistEnts = syn.getChildren("syn9689663")
    for ents in caseListEntities:
        syn.store(synapseclient.Link(
            ents.id,
            parent=caselistId,
            targetVersion=ents.versionNumber))

    # Store gene panels
    for ents in genePanelEntities:
        syn.store(synapseclient.Link(
            ents.id,
            parent=secondReleaseFolderSynId,
            targetVersion=ents.versionNumber))


def main(genie_version,
         processing_date,
         cbioportal_path,
         oncotree_link=None,
         consortium_release_cutoff=184,
         pemfile=None,
         test=False,
         staging=False,
         debug=False,
         skip_mutationsincis=False):

    syn = process.synLogin(pemfile, debug=debug)
    genie_user = os.environ['GENIE_USER']
    if pemfile is not None:
        genie_pass = process.get_password(pemfile)
    else:
        genie_pass = None

    assert not (test and staging), \
        "You can only specify --test or --staging, not both"

    if test:
        databaseSynIdMappingId = 'syn11600968'
        genie_version = "TESTING"
    elif staging:
        skip_mutationsincis = True
        databaseSynIdMappingId = 'syn12094210'
    else:
        databaseSynIdMappingId = 'syn10967259'
    # Database/folder syn id mapping
    databaseSynIdMapping = syn.tableQuery(
        'select * from {}'.format(databaseSynIdMappingId))
    databaseSynIdMappingDf = databaseSynIdMapping.asDataFrame()
    # databaseSynIdMappingDf.index = databaseSynIdMappingDf.Database
    # del databaseSynIdMappingDf['Database']
    # databaseSynIdMappingDf.to_dict()

    if oncotree_link is None:
        oncoLink = databaseSynIdMappingDf['Id'][
            databaseSynIdMappingDf['Database'] == 'oncotreeLink'].values[0]
        oncoLinkEnt = syn.get(oncoLink)
        oncotree_link = oncoLinkEnt.externalURL

    # Check if you can connect to oncotree link,
    # if not then don't run validation / processing
    process.checkUrl(oncotree_link)

    consortiumSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'consortium'].values[0]
    processTrackerSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'processTracker'].values[0]
    # get syn id of case list folder in consortium release
    caseListSynId = findCaseListId(syn, consortiumSynId)

    if not staging:
        update_process_trackingdf(
            syn, processTrackerSynId, 'SAGE', 'dbToStage', start=True)

    syn.table_query_timeout = 50000
    centerMappingSynId = databaseSynIdMappingDf['Id'][
        databaseSynIdMappingDf['Database'] == 'centerMapping'].values[0]
    # Only release files where release is true
    center_mapping = syn.tableQuery(
        'SELECT * FROM {} where release is true'.format(centerMappingSynId))
    center_mappingdf = center_mapping.asDataFrame()
    processingDate = datetime.datetime.strptime(processing_date, '%b-%Y')

    cbioValidatorPath = os.path.join(
        cbioportal_path, "core/src/main/scripts/importer/validateData.py")
    assert os.path.exists(cbioValidatorPath),\
        "Please specify correct cbioportalPath"

    logger.info("STAGING TO CONSORTIUM")
    genePanelEntities = stagingToCbio(
        syn,
        processingDate,
        genie_version,
        center_mappingdf,
        databaseSynIdMappingDf,
        oncotree_url=oncotree_link,
        consortiumReleaseCutOff=consortium_release_cutoff,
        current_release_staging=staging,
        skipMutationsInCis=skip_mutationsincis,
        test=test,
        genie_user=genie_user,
        genie_pass=genie_pass)

    # Create case lists files
    logger.info("CREATE CASE LIST FILES")
    # Remove old caselists first
    if not os.path.exists(CASE_LIST_PATH):
        os.mkdir(CASE_LIST_PATH)
    caselists = os.listdir(CASE_LIST_PATH)
    for caselist in caselists:
        os.remove(os.path.join(CASE_LIST_PATH, caselist))
    clinical_path = os.path.join(
        GENIE_RELEASE_DIR,
        'data_clinical_{}.txt'.format(genie_version))
    gene_matrix_path = os.path.join(
        GENIE_RELEASE_DIR,
        "data_gene_matrix_{}.txt".format(genie_version))
    create_case_lists.main(
        clinical_path,
        gene_matrix_path,
        CASE_LIST_PATH,
        "genie_private")
    caseListFiles = os.listdir(CASE_LIST_PATH)
    caseListEntities = []
    for casePath in caseListFiles:
        casePath = os.path.join(CASE_LIST_PATH, casePath)
        caseListEntities.append(storeFile(
            syn,
            casePath,
            parent=caseListSynId,
            staging=staging,
            caseLists=True,
            genieVersion=genie_version))

    logger.info("REMOVING UNNECESSARY FILES")
    genie_files = os.listdir(GENIE_RELEASE_DIR)
    for genieFile in genie_files:
        if genie_version not in genieFile and \
             "meta" not in genieFile and "case_lists" not in genieFile:
            os.remove(os.path.join(GENIE_RELEASE_DIR, genieFile))
    os.remove(clinical_path)

    logger.info("REVISE METADATA FILES")
    reviseMetadataFiles(syn, staging, consortiumSynId, genie_version)
    logger.info("CBIO VALIDATION")
    '''
    Must be exit 0 because the validator sometimes fails,
    but we still want to capture the output
    '''
    command = [cbioValidatorPath, '-s', GENIE_RELEASE_DIR, '-n', '; exit 0']
    cbioOutput = subprocess.check_output(" ".join(command), shell=True)
    logger.info(cbioOutput.decode("utf-8"))
    cbio_validator_log = \
        "cbioValidatorLogsConsortium_{}.txt".format(genie_version)
    if not test and not staging:
        with open(cbio_validator_log, "w") as cbioLog:
            cbioLog.write(cbioOutput.decode("utf-8"))
        syn.store(synapseclient.File(
            cbio_validator_log, parentId="syn10155804"))
        os.remove(cbio_validator_log)
    logger.info("REMOVING OLD FILES")

    process.rmFiles(CASE_LIST_PATH)
    private_cna_meta_path = \
        '%s/genie_private_meta_cna_hg19_seg.txt' % GENIE_RELEASE_DIR
    if os.path.exists(private_cna_meta_path):
        os.unlink(private_cna_meta_path)
    logger.info("CREATING LINK VERSION")
    createLinkVersion(
        syn, genie_version, caseListEntities,
        genePanelEntities, databaseSynIdMappingDf)

    if not staging:
        update_process_trackingdf(
            syn, processTrackerSynId, 'SAGE', 'dbToStage', start=False)

    logger.info("COMPLETED DATABASE TO STAGING")

    if not test:
        logger.info("DASHBOARD UPDATE")
        dashboard_table_updater.run_dashboard(
            syn,
            databaseSynIdMappingDf,
            genie_version,
            staging=staging)
        dashboard_markdown_html_commands = [
            'Rscript',
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'dashboard_markdown_generator.R'),
            genie_version]

        if genie_user is not None and genie_pass is not None:
            dashboard_markdown_html_commands.extend(
                ['--syn_user', genie_user, '--syn_pass', genie_pass])
        if staging:
            dashboard_markdown_html_commands.append('--staging')
        subprocess.check_call(dashboard_markdown_html_commands)
        logger.info("DASHBOARD UPDATE COMPLETE")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Release GENIE consortium files')

    parser.add_argument(
        "processingDate",
        type=str,
        metavar="Jan-2017",
        help="The processing date of GENIE in Month-Year format"
             " (ie. Apr-2017)")

    parser.add_argument(
        "cbioportalPath",
        type=str,
        metavar="/path/to/cbioportal",
        help="Make sure you clone the cbioportal github: "
             "git clone https://github.com/cBioPortal/cbioportal.git")

    parser.add_argument(
        "genieVersion",
        type=str,
        metavar="1.0.1",
        help="GENIE release version")

    parser.add_argument(
        "--oncotreeLink",
        type=str,
        help="Link to oncotree code")

    parser.add_argument(
        "--consortiumReleaseCutOff",
        type=int,
        metavar=184,
        default=184,
        help="Consortium release cut off time in days")

    parser.add_argument(
        "--staging",
        action='store_true',
        help="Store into staging folder")

    parser.add_argument(
        "--skipMutationsInCis",
        action='store_true',
        help="Skip running mutation in cis script")

    parser.add_argument(
        "--test",
        action='store_true',
        help="Run test")

    parser.add_argument(
        "--pemFile",
        type=str,
        help="Path to PEM file (genie.pem)")

    parser.add_argument(
        "--debug",
        action='store_true',
        help="Synapse debug feature")
    args = parser.parse_args()

    main(genie_version=args.genieVersion,
         processing_date=args.processingDate,
         cbioportal_path=args.cbioportalPath,
         oncotree_link=args.oncotreeLink,
         consortium_release_cutoff=args.consortiumReleaseCutOff,
         pemfile=args.pemFile,
         test=args.test,
         staging=args.staging,
         debug=args.debug,
         skip_mutationsincis=args.skipMutationsInCis)
