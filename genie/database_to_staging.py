#!/usr/local/bin/ python3

# import argparse
# import datetime
import logging
import math
import os
import re
import subprocess
import time

import pandas as pd
import synapseclient
import synapseutils

from . import process_functions
# from . import create_case_lists
# from . import dashboard_table_updater

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

    Returns:
        string: Synapse id of case list
    '''
    releaseEnts = synapseutils.walk(syn, parentId)
    releaseFolders = next(releaseEnts)
    if len(releaseFolders[1]) == 0:
        caselistId = syn.store(
            synapseclient.Folder(name="case_lists", parent=parentId)).id
    else:
        caselistId = releaseFolders[1][0][1]
    return(caselistId)


def storeFile(syn,
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
    #   process_functions.center_anon(filePath, ANONYMIZE_CENTER_DF)
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
    #   process_functions.center_convert_back(filePath, ANONYMIZE_CENTER_DF)
    return(temp)


def redaction_phi(values):
    '''
    Boolean vector of pediatric and PHI rows

    Args:
        values: Dataframe column or list of values to check

    Returns:
        list: redact and bool vectors
        list: pediatric redaction
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
    '''
    Re annotate PHI data

    Args:
        mergedClinical: merged clinical dataframe

    Returns:
        pandas.DataFrame: Re-annotated clinical file
    '''
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
def configureMafRow(rowArray, headers, keepSamples, remove_variants,
                    flagged_variants):
    chrom = str(rowArray[headers.index('Chromosome')])
    start = str(rowArray[headers.index('Start_Position')])
    end = str(rowArray[headers.index('End_Position')])
    ref = str(rowArray[headers.index('Reference_Allele')])
    seq = str(rowArray[headers.index('Tumor_Seq_Allele2')])
    sampleId = str(rowArray[headers.index('Tumor_Sample_Barcode')])
    hgvsp = str(rowArray[headers.index('HGVSp_Short')])
    variant = chrom+' '+start+' '+end+' '+ref+' '+seq+' '+sampleId
    # Add this line for now because merge check uses
    # different primary key from maf
    mergecheck_variant = \
        chrom+' '+start+' '+hgvsp+' '+ref+' '+seq+' '+sampleId
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
        if mergecheck_variant in flagged_variants.tolist():
            rowArray.append('True')
        else:
            rowArray.append('')
        newRow = "\t".join(rowArray)
        newRow += "\n"
        newRow = process_functions.removeStringFloat(newRow)
        return(newRow)
    else:
        return(None)


def runMAFinBED(syn,
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

        # mafText = process_functions.removePandasDfFloat(center_mutation)
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
    removeSeqDateSamples = process_functions.seqDateFilter(
        clinicalDf, processingDate, consortiumReleaseCutOff)
    return(removeSeqDateSamples)


def mutation_in_cis_filter(syn,
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
    sample_filtering = syn.tableQuery(
        "SELECT Tumor_Sample_Barcode FROM {} where Flag = 'TOSS' and "
        "Tumor_Sample_Barcode is not null".format(variant_filtering_synId))

    filtered_samplesdf = sample_filtering.asDataFrame()
    # #Alex script #1 removed patients
    remove_samples = \
        filtered_samplesdf['Tumor_Sample_Barcode'].drop_duplicates()

    # Find variants to flag
    variant_flagging = syn.tableQuery(
        "SELECT * FROM {} where Flag = 'FLAG' and "
        "Tumor_Sample_Barcode is not null".format(variant_filtering_synId))
    flag_variantsdf = variant_flagging.asDataFrame()

    flag_variantsdf['flaggedVariants'] = \
        flag_variantsdf['Chromosome'].astype(str) + ' ' + \
        flag_variantsdf['Start_Position'].astype(str) + ' ' + \
        flag_variantsdf['HGVSp_Short'].astype(str) + ' ' + \
        flag_variantsdf['Reference_Allele'].astype(str) + ' ' + \
        flag_variantsdf['Tumor_Seq_Allele2'].astype(str) + ' ' + \
        flag_variantsdf['Tumor_Sample_Barcode'].astype(str)
    return(remove_samples, flag_variantsdf['flaggedVariants'])


def seq_assay_id_filter(clinicaldf):
    '''
    (Deprecated)
    Remove samples that are part of SEQ_ASSAY_IDs with less
    than 50 samples

    Args:
        clinicalDf: Sample clinical dataframe

    Returns:
        pd.Series: samples to remove
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
        pd.Series: samples to remove
    '''

    logger.info("NO GENE PANEL FILTER")
    has_seq_assay = clinicaldf['SEQ_ASSAY_ID'].isin(beddf['SEQ_ASSAY_ID'])
    remove_samples = clinicaldf['SAMPLE_ID'][~has_seq_assay]
    logger.info("Removing samples with no bed file: {}".format(
        ",".join(remove_samples)))
    return(remove_samples)


def store_gene_panel_files(syn,
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


def return_syn_tablequerydf(syn, query):
    table = syn.tableQuery(query)
    return(table.asDataFrame())


def store_fusion_files(syn,
                       release_synid,
                       genie_version,
                       fusion_synid,
                       keep_for_center_consortium_samples,
                       keep_for_merged_consortium_samples,
                       current_release_staging,
                       center_mappingdf):
    '''
    Create, filter, configure, and store fusion file

    Args:
        syn: Synapse object
        release_synid: Synapse id to store release file
        genie_version: GENIE version (ie. v6.1-consortium)
        fusion_synid: Fusion database synid
        keep_for_center_consortium_samples: Samples to keep for center files
        keep_for_merged_consortium_samples: Samples to keep for merged file
        current_release_staging: Staging flag
        center_mappingdf: Center mapping dataframe
    '''
    logger.info("MERING, FILTERING, STORING FUSION FILES")
    FusionsDf = return_syn_tablequerydf(
        syn,
        'select HUGO_SYMBOL,ENTREZ_GENE_ID,CENTER,TUMOR_SAMPLE_BARCODE,FUSION,'
        'DNA_SUPPORT,RNA_SUPPORT,METHOD,FRAME from {}'.format(fusion_synid))
    # FusionsDf = Fusions.asDataFrame()
    FusionsDf['ENTREZ_GENE_ID'][FusionsDf['ENTREZ_GENE_ID'] == 0] = pd.np.nan

    if not current_release_staging:
        FusionsStagingDf = FusionsDf[FusionsDf['TUMOR_SAMPLE_BARCODE'].isin(
            keep_for_center_consortium_samples)]
        for center in center_mappingdf.center:
            center_fusion = FusionsStagingDf[
                FusionsStagingDf['CENTER'] == center]
            if not center_fusion.empty:
                center_fusion.to_csv(
                    FUSIONS_CENTER_PATH % center,
                    sep="\t", index=False)
                storeFile(
                    syn,
                    FUSIONS_CENTER_PATH % center,
                    genieVersion=genie_version,
                    parent=center_mappingdf['stagingSynId'][
                        center_mappingdf['center'] == center][0],
                    centerStaging=True)

    FusionsDf = FusionsDf[FusionsDf['TUMOR_SAMPLE_BARCODE'].isin(
        keep_for_merged_consortium_samples)]
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
    fusionText = process_functions.removePandasDfFloat(FusionsDf)
    fusions_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_fusions_%s.txt' % genie_version)
    with open(fusions_path, "w") as fusionFile:
        fusionFile.write(fusionText)
    storeFile(
        syn, fusions_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="data_fusions.txt",
        staging=current_release_staging)


def store_maf_files(syn,
                    genie_version,
                    flatfiles_view_synid,
                    release_synid,
                    clinicaldf,
                    center_mappingdf,
                    keep_for_merged_consortium_samples,
                    keep_for_center_consortium_samples,
                    remove_mafinbed_variants,
                    flagged_mutationInCis_variants,
                    current_release_staging):
    '''
    Create, filter, configure, and store maf file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        flatfiles_view_synid: Synapse id of fileview with all the flat files
        release_synid: Synapse id to store release file
        clinicaldf: Clinical dataframe with SAMPLE_ID and CENTER
        center_mappingdf: Center mapping dataframe
        keep_for_merged_consortium_samples: Samples to keep for merged file
        keep_for_center_consortium_samples: Samples to keep for center files
        remove_mafinbed_variants: Variants to remove
        flagged_mutationInCis_variants: Variants to flag
        current_release_staging: Staging flag
    '''

    logger.info("FILTERING, STORING MUTATION FILES")
    centerMafSynIds = syn.tableQuery(
        "select id from {} where name like '%mutation%'".format(
            flatfiles_view_synid))
    centerMafSynIdsDf = centerMafSynIds.asDataFrame()
    mutations_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_mutations_extended_%s.txt' % genie_version)
    with open(mutations_path, 'w') as f:
        pass
    for index, mafSynId in enumerate(centerMafSynIdsDf.id):
        mafEnt = syn.get(mafSynId)
        logger.info(mafEnt.path)
        with open(mafEnt.path, "r") as mafFile:
            header = mafFile.readline()
            headers = header.replace("\n", "").split("\t")
            # Add in mutation in cis flag header
            headers.append("mutationInCis_Flag")
            header = "\t".join(headers) + "\n"
            if index == 0:
                with open(mutations_path, 'a') as f:
                    f.write(header)
                # Create maf file per center for their staging directory
                for center in clinicaldf['CENTER'].unique():
                    with open(MUTATIONS_CENTER_PATH % center, 'w') as f:
                        f.write(header)
        # with open(mafEnt.path,"r") as newMafFile:
        #   newMafFile.readline()
            center = mafEnt.path.split("_")[3]
            # Make sure to only write the centers that release = True
            if center in center_mappingdf.center.tolist():
                for row in mafFile:
                    rowArray = row.replace("\n", "").split("\t")
                    center = rowArray[headers.index('Center')]
                    newMergedRow = configureMafRow(
                        rowArray, headers,
                        keep_for_merged_consortium_samples,
                        remove_mafinbed_variants,
                        flagged_mutationInCis_variants)
                    if newMergedRow is not None:
                        with open(mutations_path, 'a') as f:
                            f.write(newMergedRow)
                    newCenterRow = configureMafRow(
                        rowArray, headers,
                        keep_for_center_consortium_samples,
                        remove_mafinbed_variants,
                        flagged_mutationInCis_variants)
                    if newCenterRow is not None:
                        with open(MUTATIONS_CENTER_PATH % center, 'a') as f:
                            f.write(newCenterRow)
    storeFile(
        syn, mutations_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="data_mutations_extended.txt",
        staging=current_release_staging)

    if not current_release_staging:
        for center in clinicaldf['CENTER'].unique():
            storeFile(
                syn, MUTATIONS_CENTER_PATH % center,
                genieVersion=genie_version,
                parent=center_mappingdf['stagingSynId'][
                    center_mappingdf['center'] == center][0],
                centerStaging=True)


def run_genie_filters(syn,
                      genie_user,
                      genie_pass,
                      genie_version,
                      variant_filtering_synId,
                      clinicaldf,
                      beddf,
                      sample_cols,
                      center_mappingdf,
                      processing_date,
                      skip_mutationsincis,
                      consortium_release_cutoff,
                      test):
    '''
    Run GENIE filters and returns variants and samples to remove

    Args:
        syn: Synapse object
        genie_user: Synapse username
        genie_pass: Synapse password
        genie_version: GENIE version (ie. v6.1-consortium)
        variant_filtering_synId: Synapse id of mutationInCis table
        clinicaldf: Clinical dataframe with SAMPLE_ID and SEQ_ASSAY_ID
        beddf: Bed dataframe
        sample_cols: Clinical sample columns
        center_mappingdf: Center mapping dataframe
        processing_date: Processing date
        skip_mutationsincis: Skip mutation in cis filter
        consortium_release_cutoff: Release cutoff in days
        test: Test flag

    Returns:
        pandas.Series: Variants to remove
        set: samples to remove for release files
        set: samples to remove for center files
        pandas.Series: Variants to flag
    '''

    # ADD CHECKS TO CODE BEFORE UPLOAD.
    # Throw error if things don't go through
    logger.info("RUN GENIE FILTERS")
    # STORING CLINICAL FILES INTO CBIOPORTAL

    ''' FILTERING '''
    logger.info("MAF IN BED FILTER")
    remove_mafinbed_variants = runMAFinBED(
        syn, center_mappingdf, test=test,
        genieVersion=genie_version,
        genie_user=genie_user, genie_pass=genie_pass)

    logger.info("MUTATION IN CIS FILTER")
    remove_mutationincis_samples, flagged_mutationincis_variants = \
        mutation_in_cis_filter(
            syn, skip_mutationsincis, variant_filtering_synId,
            center_mappingdf, genieVersion=genie_version, test=test,
            genie_user=genie_user, genie_pass=genie_pass)
    remove_no_genepanel_samples = no_genepanel_filter(clinicaldf, beddf)

    logger.info("SEQ DATE FILTER")
    remove_seqdate_samples = seq_date_filter(
        clinicaldf, processing_date, consortium_release_cutoff)
    # Only certain samples are removed for the files that go into
    # staging center folder
    remove_center_consortium_samples = set(remove_mutationincis_samples).union(
        set(remove_no_genepanel_samples))
    # Most filteres are applied for the files that go into the merged
    # consortium release
    remove_merged_consortium_samples = set(remove_seqdate_samples)
    # set(remove_seqAssayId_samples)#.union(set(remove_seqDate_samples))
    remove_merged_consortium_samples = \
        remove_merged_consortium_samples.union(
            remove_center_consortium_samples)

    return(remove_mafinbed_variants,
           remove_merged_consortium_samples,
           remove_center_consortium_samples,
           flagged_mutationincis_variants)


def store_clinical_files(syn,
                         genie_version,
                         clinicaldf,
                         oncotree_url,
                         sample_cols,
                         patient_cols,
                         remove_center_consortium_samples,
                         remove_merged_consortium_samples,
                         release_synid,
                         current_release_staging,
                         center_mappingdf):
    '''
    Create, filter, configure, and store clinical file

    Args:
        syn: Synapse object
        genie_version: GENIE version (ie. v6.1-consortium)
        clinicaldf: Clinical dataframe with SAMPLE_ID and SEQ_ASSAY_ID
        oncotree_url: Oncotree URL
        sample_cols: Clinical sample columns
        patient_cols: Clinical patient columns
        remove_center_consortium_samples: Samples to remove for center files
        remove_merged_consortium_samples: Samples to remove for merged file
        release_synid: Synapse id to store release file
        current_release_staging: Staging flag
        center_mappingdf: Center mapping dataframe

    Returns:
        pandas.DataFrame: configured clinical dataframe
        pandas.Series: samples to keep for center files
        pandas.Series: samples to keep for release files
    '''

    logger.info("CONFIGURING CLINICAL FILES")
    logger.info("REMOVING PHI")
    clinicaldf = reAnnotatePHI(clinicaldf)
    logger.info("ADD CANCER TYPES")
    # This removes support for both oncotree urls (only support json)
    oncotree_dict = process_functions.get_oncotree_code_mappings(oncotree_url)
    # Add in unknown key which maps to UNKNOWN everything
    oncotree_dict['UNKNOWN'] = {
        'CANCER_TYPE': 'UNKNOWN',
        'CANCER_TYPE_DETAILED': 'UNKNOWN',
        'ONCOTREE_PRIMARY_NODE': 'UNKNOWN',
        'ONCOTREE_SECONDARY_NODE': 'UNKNOWN'}

    clinicaldf['CANCER_TYPE'] = [
        oncotree_dict[code.upper()]["CANCER_TYPE"]
        if code.upper() in oncotree_dict.keys() else float('nan')
        for code in clinicaldf['ONCOTREE_CODE']]

    clinicaldf['CANCER_TYPE_DETAILED'] = [
        oncotree_dict[code.upper()]["CANCER_TYPE_DETAILED"]
        if code.upper() in oncotree_dict.keys() else float('nan')
        for code in clinicaldf['ONCOTREE_CODE']]

    clinicaldf['ONCOTREE_PRIMARY_NODE'] = [
        oncotree_dict[code.upper()]["ONCOTREE_PRIMARY_NODE"]
        if code.upper() in oncotree_dict.keys() else float('nan')
        for code in clinicaldf['ONCOTREE_CODE']]

    clinicaldf['ONCOTREE_SECONDARY_NODE'] = [
        oncotree_dict[code.upper()]["ONCOTREE_SECONDARY_NODE"]
        if code.upper() in oncotree_dict.keys() else float('nan')
        for code in clinicaldf['ONCOTREE_CODE']]

    # All cancer types that are null should have null oncotree codes
    clinicaldf['ONCOTREE_CODE'][
        clinicaldf['CANCER_TYPE'].isnull()] = float('nan')
    # Suggest using AGE_AT_SEQ_REPORT_DAYS instead so that the
    # descriptions can match
    clinicaldf['AGE_AT_SEQ_REPORT_DAYS'] = clinicaldf['AGE_AT_SEQ_REPORT']
    clinicaldf['AGE_AT_SEQ_REPORT'] = [
        int(math.floor(int(float(age))/365.25))
        if process_functions.checkInt(age) else age
        for age in clinicaldf['AGE_AT_SEQ_REPORT']]
    clinicaldf['AGE_AT_SEQ_REPORT'][
        clinicaldf['AGE_AT_SEQ_REPORT'] == ">32485"] = ">89"
    clinicaldf['AGE_AT_SEQ_REPORT'][
        clinicaldf['AGE_AT_SEQ_REPORT'] == "<6570"] = "<18"

    ############################################################
    # CENTER SPECIFIC CODE FOR RIGHT NOW (REMOVE UHN-555-V1)
    ############################################################
    # clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'] != "UHN-555-V1"]
    # clinicalDf = clinicalDf[clinicalDf['SEQ_ASSAY_ID'] != "PHS-TRISEQ-V1"]

    # clinicalDf = clinicalDf[clinicalDf['CENTER'] != "WAKE"]
    # clinicalDf = clinicalDf[clinicalDf['CENTER'] != "CRUK"]
    ############################################################
    ############################################################

    clinicaldf.drop_duplicates("SAMPLE_ID", inplace=True)

    logger.info("STORING CLINICAL FILES")
    # samples must be removed after reading in the clinical file again
    staging_clinicaldf = clinicaldf[
        ~clinicaldf['SAMPLE_ID'].isin(remove_center_consortium_samples)]
    if not current_release_staging:
        for center in center_mappingdf.center:
            center_clinical = \
                staging_clinicaldf[staging_clinicaldf['CENTER'] == center]
            center_sample = \
                center_clinical[sample_cols].drop_duplicates('SAMPLE_ID')
            center_patient = \
                center_clinical[patient_cols].drop_duplicates('PATIENT_ID')
            center_sample.to_csv(
                SAMPLE_CENTER_PATH % center, sep="\t", index=False)
            center_patient.to_csv(
                PATIENT_CENTER_PATH % center, sep="\t", index=False)
            storeFile(
                syn, SAMPLE_CENTER_PATH % center,
                genieVersion=genie_version,
                parent=center_mappingdf['stagingSynId'][
                    center_mappingdf['center'] == center][0],
                centerStaging=True)
            storeFile(
                syn, PATIENT_CENTER_PATH % center,
                genieVersion=genie_version,
                parent=center_mappingdf['stagingSynId'][
                    center_mappingdf['center'] == center][0],
                centerStaging=True)

    clinicaldf = clinicaldf[~clinicaldf['SAMPLE_ID'].isin(
        remove_merged_consortium_samples)]
    # This must happen here because the seq assay filter
    # must happen after all the other filters
    # logger.info("SEQ ASSAY FILTER")
    # remove_seqAssayId_samples = seq_assay_id_filter(clinicalDf)
    # removeForMergedConsortiumSamples = \
    #     removeForMergedConsortiumSamples.union(set(remove_seqAssayId_samples))
    # clinicalDf = clinicalDf[~clinicalDf['SAMPLE_ID'].isin(
    #    remove_seqAssayId_samples)]

    keep_center_consortium_samples = staging_clinicaldf.SAMPLE_ID
    keep_merged_consortium_samples = clinicaldf.SAMPLE_ID
    # This mapping table is the GENIE clinical code to description
    # mapping to generate the headers of the clinical file
    mapping_table = syn.tableQuery('SELECT * FROM syn9621600')
    mapping = mapping_table.asDataFrame()
    clinical_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_clinical_%s.txt' % genie_version)
    clinical_sample_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_clinical_sample_%s.txt' % genie_version)
    clinical_patient_path = os.path.join(
        GENIE_RELEASE_DIR, 'data_clinical_patient_%s.txt' % genie_version)
    process_functions.addClinicalHeaders(
        clinicaldf, mapping, patient_cols, sample_cols,
        clinical_sample_path, clinical_patient_path)
    storeFile(
        syn, clinical_sample_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="data_clinical_sample.txt",
        staging=current_release_staging)

    storeFile(
        syn, clinical_patient_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="data_clinical_patient.txt",
        staging=current_release_staging)

    clinicaldf.to_csv(clinical_path, sep="\t", index=False)
    storeFile(
        syn, clinical_path,
        parent=release_synid,
        name="data_clinical.txt",
        staging=current_release_staging)

    return(clinicaldf,
           keep_center_consortium_samples,
           keep_merged_consortium_samples)


def store_cna_files(syn,
                    flatfiles_view_synid,
                    keep_for_center_consortium_samples,
                    keep_for_merged_consortium_samples,
                    center_mappingdf,
                    genie_version,
                    release_synid,
                    current_release_staging):
    '''
    Create, filter and store cna file

    Args:
        syn: Synapse object
        flatfiles_view_synid: Synapse id of fileview with all the flat files
        keep_for_center_consortium_samples: Samples to keep for center files
        keep_for_merged_consortium_samples: Samples to keep for merged file
        center_mappingdf: Center mapping dataframe
        genie_version: GENIE version (ie. v6.1-consortium)
        release_synid: Synapse id to store release file
        current_release_staging: Staging flag

    Returns:
        list: CNA samples
    '''
    logger.info("MERING, FILTERING, STORING CNA FILES")
    cna_path = os.path.join(
        GENIE_RELEASE_DIR, "data_CNA_%s.txt" % genie_version)
    center_cna_synids = syn.tableQuery(
        "select id from %s " % flatfiles_view_synid +
        "where name like 'data_CNA%'")
    center_cna_synidsdf = center_cna_synids.asDataFrame()
    # Grab all unique symbols and form cna_template
    all_symbols = set()
    for cna_synid in center_cna_synidsdf['id']:
        cna_ent = syn.get(cna_synid)
        with open(cna_ent.path, "r") as cna_file:
            # Read first line first
            cna_file.readline()
            # Get all hugo symbols
            all_symbols = all_symbols.union(
                set(line.split("\t")[0] for line in cna_file))
    cna_template = pd.DataFrame({"Hugo_Symbol": list(all_symbols)})
    cna_template.sort_values("Hugo_Symbol", inplace=True)
    cna_template.to_csv(cna_path, sep="\t", index=False)
    # Loop through to create finalized CNA file
    with_center_hugo_symbol = pd.Series("Hugo_Symbol")
    with_center_hugo_symbol = with_center_hugo_symbol.append(
        pd.Series(keep_for_center_consortium_samples))

    with_merged_hugo_symbol = pd.Series("Hugo_Symbol")
    with_merged_hugo_symbol = with_merged_hugo_symbol.append(
        pd.Series(keep_for_merged_consortium_samples))

    cna_samples = []

    for cna_synId in center_cna_synidsdf['id']:
        cna_ent = syn.get(cna_synId)
        center = cna_ent.name.replace("data_CNA_", "").replace(".txt", "")
        logger.info(cna_ent.path)
        if center in center_mappingdf.center.tolist():
            center_cna = pd.read_csv(cna_ent.path, sep="\t")
            merged_cna = cna_template.merge(
                center_cna, on="Hugo_Symbol", how="outer")
            merged_cna.sort_values("Hugo_Symbol", inplace=True)

            merged_cna = merged_cna[merged_cna.columns[
                merged_cna.columns.isin(with_center_hugo_symbol)]]

            cna_text = process_functions.removePandasDfFloat(merged_cna)
            # Replace blank with NA's
            cna_text = cna_text.replace(
                "\t\t", "\tNA\t").replace(
                "\t\t", "\tNA\t").replace(
                '\t\n', "\tNA\n")

            # Store center CNA file in staging dir
            with open(CNA_CENTER_PATH % center, "w") as cna_file:
                cna_file.write(cna_text)
            storeFile(
                syn, CNA_CENTER_PATH % center,
                genieVersion=genie_version,
                parent=center_mappingdf['stagingSynId'][
                    center_mappingdf['center'] == center][0],
                centerStaging=True)
            # This is to remove more samples for the final cna file
            merged_cna = merged_cna[merged_cna.columns[
                merged_cna.columns.isin(with_merged_hugo_symbol)]]

            cna_text = process_functions.removePandasDfFloat(merged_cna)
            cna_text = cna_text.replace(
                "\t\t", "\tNA\t").replace(
                "\t\t", "\tNA\t").replace(
                '\t\n', "\tNA\n")

            with open(CNA_CENTER_PATH % center, "w") as cna_file:
                cna_file.write(cna_text)
            # Join CNA file
            cna_samples.extend(merged_cna.columns[1:].tolist())
            linux_join_command = ["join", cna_path, CNA_CENTER_PATH % center]
            output = subprocess.check_output(linux_join_command)
            with open(cna_path, "w") as cna_file:
                cna_file.write(output.decode("utf-8").replace(" ", "\t"))

    storeFile(
        syn, cna_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="data_CNA.txt",
        staging=current_release_staging)

    return(cna_samples)


# SEG
def store_seg_files(syn,
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
        staging_segdf = segdf[segdf['ID'].isin(
            keep_for_center_consortium_samples)]
        for center in center_mappingdf.center:
            center_seg = staging_segdf[staging_segdf['CENTER'] == center]
            if not center_seg.empty:
                del center_seg['CENTER']
                segtext = process_functions.removePandasDfFloat(center_seg)
                with open(SEG_CENTER_PATH % center, "w") as seg_file:
                    seg_file.write(segtext)
                storeFile(
                    syn, SEG_CENTER_PATH % center,
                    genieVersion=genie_version,
                    parent=center_mappingdf['stagingSynId'][
                        center_mappingdf['center'] == center][0],
                    centerStaging=True)
    del segdf['CENTER']
    segdf = segdf[segdf['ID'].isin(keep_for_merged_consortium_samples)]
    segtext = process_functions.removePandasDfFloat(segdf)
    with open(seg_path, "w") as seg_file:
        seg_file.write(segtext)
    storeFile(
        syn, seg_path,
        parent=release_synid,
        genieVersion=genie_version,
        name="genie_private_data_cna_hg19.seg",
        staging=current_release_staging)


def store_data_gene_matrix(syn,
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
        pandas.DataFrame: data gene matrix dataframe
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


def store_bed_files(syn,
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
            bed_seq_df = \
                bed_seq_df[bed_seq_df['Hugo_Symbol'] != bed_seq_df['ID']]
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


def stagingToCbio(syn, processingDate, genieVersion,
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

    Returns:
        list: Gene panel entities
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
        removeForCenterConsortiumSamples, \
        flagged_mutationInCis_variants = \
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
        flagged_mutationInCis_variants,
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


def update_process_trackingdf(syn, process_trackerdb_synid, center,
                              process_type, start=True):
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


def revise_metadata_files(syn, staging, consortiumid, genie_version=None):
    '''
    Rewrite metadata files with the correct GENIE version

    Args:
        syn: Synapse object
        staging: staging flag
        consortiumid: Synapse id of consortium release folder
        genie_version: GENIE version, Default to None
    '''
    release_files = syn.getChildren(consortiumid)
    meta_file_ents = [
        syn.get(i['id'],
                downloadLocation=GENIE_RELEASE_DIR,
                ifcollision="overwrite.local")
        for i in release_files if 'meta' in i['name']]

    for meta_ent in meta_file_ents:
        with open(meta_ent.path, "r+") as meta:
            meta_text = meta.read()
            if "meta_study" not in meta_ent.path:
                version = ''
            else:
                version = re.search(".+GENIE.+v(.+)", meta_text).group(1)
            # Fix this line
            genie_version = version if genie_version is None else genie_version
            version_on_file = re.search(".+data_(.+)[.]txt", meta_text)
            if version_on_file is None:
                version_on_file = re.search(".+data_(.+)[.]seg", meta_text)
            if version_on_file is not None:
                version_on_file = version_on_file.group(1).split("_")[-1]

            if version != genie_version:
                meta_text = meta_text.replace(
                    "GENIE Cohort v{}".format(version),
                    "GENIE Cohort v{}".format(genie_version))

                meta_text = meta_text.replace(
                    "GENIE v{}".format(version),
                    "GENIE v{}".format(genie_version))

                if version_on_file is not None:
                    meta_text = meta_text.replace(
                        version_on_file, genie_version)
                    meta_text = meta_text.replace(
                        version_on_file, genie_version)
                meta.seek(0)
                meta.write(meta_text)
                meta.truncate()
                storeFile(
                    syn,
                    meta_ent.path,
                    parent=consortiumid,
                    genieVersion=genie_version,
                    staging=staging)


def search_and_create_folder(syn, parentid, folder_name):
    '''
    This function will search for an existing directory given a parent id
    It will create the Synapse folder if it doesn't exist

    Args:
        syn: Synapse object
        parentid: Synapse id of a project or folder
        folder_name: Folder being searched for

    Returns:
        string: Synapse folder id
        boolean: Checks if the folder already existed.  True if it did
    '''
    folders = syn.getChildren(parentid, includeTypes=['folder'])
    search_for_folder = \
        filter(lambda folder: folder['name'] == folder_name, folders)
    search_for_folder = list(search_for_folder)
    if len(search_for_folder) == 0:
        folder_ent = synapseclient.Folder(name=folder_name, parent=parentid)
        folder_synid = syn.store(folder_ent)['id']
        already_exists = False
    elif len(search_for_folder) == 1:
        folder_synid = search_for_folder[0]['id']
        already_exists = True
    else:
        raise ValueError("There should not be any duplicated folder names")
    return(folder_synid, already_exists)


def create_link_version(syn,
                        genie_version,
                        case_list_entities,
                        gene_panel_entities,
                        database_synid_mappingdf):
    '''
    Create release links from the actual entity and version

    TODO: Refactor to use fileviews

    Args:
        syn: Synapse object
        genie_version: GENIE version number
        case_list_entities: Case list entities
        gene_panel_entities: Gene panel entities
        database_synid_mappingdf: dataframe containing database to
                                  synapse id mapping
    '''
    # Grab major release numbers (ie 1,2,3 ...)
    major_release = genie_version.split(".")[0]
    all_releases_synid = database_synid_mappingdf['Id'][
        database_synid_mappingdf['Database'] == 'release'].values[0]
    # Create major release folder
    major_release_folder_synid, already_exists = search_and_create_folder(
        syn, all_releases_synid, "Release {}".format(major_release))
    # If the major release folder didn't exist, go ahead and create the
    # release folder
    if not already_exists:
        release_folder_ent = synapseclient.Folder(
            genie_version,
            parent=major_release_folder_synid)
        release_folder_synid = syn.store(release_folder_ent)['id']
    else:
        release_folder_synid, already_exists = search_and_create_folder(
            syn, major_release_folder_synid, genie_version)
    # Search or create case lists folder
    caselist_folder_synid, already_exists = search_and_create_folder(
        syn, release_folder_synid, "case_lists")

    # caselistId = findCaseListId(syn, release_folder_synid)
    consortium_synid = database_synid_mappingdf['Id'][
        database_synid_mappingdf['Database'] == 'consortium'].values[0]
    consortium_release_files = syn.getChildren(consortium_synid)
    # data_clinical.txt MUST be pulled in because the clinical file is
    # needed in the consortium_to_public.py
    for release_file in consortium_release_files:
        if release_file['type'] != "org.sagebionetworks.repo.model.Folder" \
                and not release_file['name'].startswith("data_gene_panel"):
            syn.store(synapseclient.Link(
                release_file['id'],
                parent=release_folder_synid,
                targetVersion=release_file['versionNumber']))

    release_files = syn.getChildren(release_folder_synid)
    clinical_ent = [
        ents['id']
        for ents in release_files
        if ents['name'] == "data_clinical.txt"][0]
    # Set private permission for the data_clinical.txt link
    syn.setPermissions(clinical_ent, principalId=3346558, accessType=[])
    syn.setPermissions(clinical_ent, principalId=3326313, accessType=[])

    for ents in case_list_entities:
        syn.store(synapseclient.Link(
            ents.id,
            parent=caselist_folder_synid,
            targetVersion=ents.versionNumber))

    # Store gene panels
    for ents in gene_panel_entities:
        syn.store(synapseclient.Link(
            ents.id,
            parent=release_folder_synid,
            targetVersion=ents.versionNumber))
