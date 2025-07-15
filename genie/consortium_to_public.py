"""Converts consortium release files to public release files"""

import logging
import os
import sys

import pandas as pd
import synapseclient
import synapseutils
from genie import (
    create_case_lists,
    database_to_staging,
    extract,
    load,
    process_functions,
)
from synapseclient import Entity, File, Folder, Link, Project, Schema

logger = logging.getLogger(__name__)
stdout_handler = logging.StreamHandler(stream=sys.stdout)
stdout_handler.setLevel(logging.INFO)
logger.addHandler(stdout_handler)
from typing import Dict

import synapseutils as synu


# TODO: Add to transform.py
def commonVariantFilter(mafDf):
    """
    This filter returns variants to keep

    Args:
        mafDf: Maf dataframe
    """
    mafDf["FILTER"] = mafDf["FILTER"].fillna("")
    to_keep = ["common_variant" not in i for i in mafDf["FILTER"]]
    mafDf = mafDf[to_keep]
    return mafDf


def _copyRecursive(
    syn: synapseclient.Synapse,
    entity: str,
    destinationId: str,
    mapping: Dict[str, str] = None,
    skipCopyAnnotations: bool = False,
    **kwargs,
) -> Dict[str, str]:
    """
    NOTE: This is a copy of the function found here: https://github.com/Sage-Bionetworks/synapsePythonClient/blob/develop/synapseutils/copy_functions.py#L409
    This was copied because there is a restriction that doesn't allow for copying entities with access requirements

    Recursively copies synapse entites, but does not copy the wikis

    Arguments:
        syn: A Synapse object with user's login
        entity: A synapse entity ID
        destinationId: Synapse ID of a folder/project that the copied entity is being copied to
        mapping: A mapping of the old entities to the new entities
        skipCopyAnnotations: Skips copying the annotations
                                Default is False

    Returns:
        a mapping between the original and copied entity: {'syn1234':'syn33455'}
    """

    version = kwargs.get("version", None)
    setProvenance = kwargs.get("setProvenance", "traceback")
    excludeTypes = kwargs.get("excludeTypes", [])
    updateExisting = kwargs.get("updateExisting", False)
    if mapping is None:
        mapping = dict()
    # Check that passed in excludeTypes is file, table, and link
    if not isinstance(excludeTypes, list):
        raise ValueError("Excluded types must be a list")
    elif not all([i in ["file", "link", "table"] for i in excludeTypes]):
        raise ValueError(
            "Excluded types can only be a list of these values: file, table, and link"
        )

    ent = syn.get(entity, downloadFile=False)
    if ent.id == destinationId:
        raise ValueError("destinationId cannot be the same as entity id")

    if (isinstance(ent, Project) or isinstance(ent, Folder)) and version is not None:
        raise ValueError("Cannot specify version when copying a project of folder")

    if not isinstance(ent, (Project, Folder, File, Link, Schema, Entity)):
        raise ValueError("Not able to copy this type of file")

    permissions = syn.restGET("/entity/{}/permissions".format(ent.id))
    # Don't copy entities without DOWNLOAD permissions
    if not permissions["canDownload"]:
        syn.logger.warning(
            "%s not copied - this file lacks download permission" % ent.id
        )
        return mapping

    # HACK: These lines of code were removed to allow for data with access requirements to be copied
    # https://github.com/Sage-Bionetworks/synapsePythonClient/blob/2909fa778e814f62f6fe6ce2d951ce58c0080a4e/synapseutils/copy_functions.py#L464-L470

    copiedId = None

    if isinstance(ent, Project):
        project = syn.get(destinationId)
        if not isinstance(project, Project):
            raise ValueError(
                "You must give a destinationId of a new project to copy projects"
            )
        copiedId = destinationId
        # Projects include Docker repos, and Docker repos cannot be copied
        # with the Synapse rest API. Entity views currently also aren't
        # supported
        entities = syn.getChildren(
            entity, includeTypes=["folder", "file", "table", "link"]
        )
        for i in entities:
            mapping = _copyRecursive(
                syn,
                i["id"],
                destinationId,
                mapping=mapping,
                skipCopyAnnotations=skipCopyAnnotations,
                **kwargs,
            )

        if not skipCopyAnnotations:
            project.annotations = ent.annotations
            syn.store(project)
    elif isinstance(ent, Folder):
        copiedId = synu.copy_functions._copyFolder(
            syn,
            ent.id,
            destinationId,
            mapping=mapping,
            skipCopyAnnotations=skipCopyAnnotations,
            **kwargs,
        )
    elif isinstance(ent, File) and "file" not in excludeTypes:
        copiedId = synu.copy_functions._copyFile(
            syn,
            ent.id,
            destinationId,
            version=version,
            updateExisting=updateExisting,
            setProvenance=setProvenance,
            skipCopyAnnotations=skipCopyAnnotations,
        )
    elif isinstance(ent, Link) and "link" not in excludeTypes:
        copiedId = synu.copy_functions._copyLink(
            syn, ent.id, destinationId, updateExisting=updateExisting
        )
    elif isinstance(ent, Schema) and "table" not in excludeTypes:
        copiedId = synu.copy_functions._copyTable(
            syn, ent.id, destinationId, updateExisting=updateExisting
        )
    # This is currently done because copyLink returns None sometimes
    if copiedId is not None:
        mapping[ent.id] = copiedId
        syn.logger.info("Copied %s to %s" % (ent.id, copiedId))
    else:
        syn.logger.info("%s not copied" % ent.id)
    return mapping


# TODO: Add to etl.py
def consortiumToPublic(
    syn,
    processingDate,
    genie_version,
    releaseId,
    databaseSynIdMappingDf,
    publicReleaseCutOff=365,
):
    cna_path = os.path.join(database_to_staging.GENIE_RELEASE_DIR, "data_CNA.txt")
    clinical_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR, "data_clinical.txt"
    )
    clinical_sample_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "data_clinical_sample.txt",
    )
    clinicl_patient_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "data_clinical_patient.txt",
    )
    data_gene_panel_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR, "data_gene_matrix.txt"
    )
    mutations_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "data_mutations_extended.txt",
    )
    seg_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR,
        "data_cna_hg19.seg",
    )
    combined_bed_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR, "genie_combined.bed"
    )

    if not os.path.exists(database_to_staging.GENIE_RELEASE_DIR):
        os.mkdir(database_to_staging.GENIE_RELEASE_DIR)
    if not os.path.exists(database_to_staging.CASE_LIST_PATH):
        os.mkdir(database_to_staging.CASE_LIST_PATH)

    # public release preview
    public_release_preview = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "public"
    ].values[0]
    public_release_preview_caselist = database_to_staging.search_or_create_folder(
        syn=syn, parentid=public_release_preview, folder_name="case_lists"
    )

    #######################################################################
    # Sponsored projects filter
    #######################################################################
    # if before release date -> go into staging consortium
    # if after date -> go into public
    # sponsoredReleaseDate = syn.tableQuery('SELECT * FROM syn8545108')
    # sponsoredReleaseDateDf = sponsoredReleaseDate.asDataFrame()
    # sponsoredProjectSamples = syn.tableQuery('SELECT * FROM syn8545106')
    # sponsoredProjectSamplesDf = sponsoredProjectSamples.asDataFrame()
    # sponsoredProjectsDf = sponsoredProjectSamplesDf.merge(
    #     sponsoredReleaseDateDf, left_on="sponsoredProject",
    #     right_on="sponsoredProjects")
    # dates = sponsoredProjectsDf['releaseDate'].apply(
    #     lambda date: datetime.datetime.strptime(date, '%b-%Y'))
    # publicReleaseSamples = sponsoredProjectsDf['genieSampleId'][
    #     dates < processingDate]
    #######################################################################

    # SEQ_DATE filter
    # Jun-2015, given processing date (today) -> public release
    # (processing date - Jun-2015 > 12 months)
    consortiumReleaseWalk = synapseutils.walk(syn, releaseId)

    consortiumRelease = next(consortiumReleaseWalk)
    for filename, synid in consortiumRelease[2]:
        if filename == "data_clinical.txt":
            clinical = syn.get(synid, followLink=True)
        elif filename == "data_gene_matrix.txt":
            gene_matrix = syn.get(synid, followLink=True)
        elif filename == "assay_information.txt":
            assay_info = syn.get(synid, followLink=True)

    clinicalDf = pd.read_csv(clinical.path, sep="\t", comment="#")
    gene_matrixdf = pd.read_csv(gene_matrix.path, sep="\t")

    removeForPublicSamples = process_functions.seqDateFilter(
        clinicalDf, processingDate, publicReleaseCutOff
    )
    logger.info("SAMPLE CLASS FILTER")
    # comment back in when public release filter back on
    # publicReleaseSamples = publicReleaseSamples.append(keepForPublicSamples)
    # Make sure all null oncotree codes are removed
    clinicalDf = clinicalDf[~clinicalDf["ONCOTREE_CODE"].isnull()]
    publicReleaseSamples = clinicalDf.SAMPLE_ID[
        ~clinicalDf.SAMPLE_ID.isin(removeForPublicSamples)
    ]

    existing_seq_dates = clinicalDf.SEQ_DATE[
        clinicalDf.SAMPLE_ID.isin(publicReleaseSamples)
    ]

    logger.info(
        "SEQ_DATES for public release: "
        + ", ".join(set(existing_seq_dates.astype(str)))
    )

    # Clinical release scope filter
    clinical_tier_release_scope_synid = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "clinical_tier_release_scope"
    ][0]
    publicRelease = extract.get_syntabledf(
        syn=syn,
        query_string=f"SELECT * FROM {clinical_tier_release_scope_synid} where releaseScope = 'public'",
    )

    # check if SAMPLE_CLASS is present
    if not process_functions.check_values_in_column(
        publicRelease, "fieldName", "SAMPLE_CLASS"
    ):
        logger.error("Must have SAMPLE_CLASS column in the public release scope.")

    allClin = clinicalDf[clinicalDf["SAMPLE_ID"].isin(publicReleaseSamples)]
    # check if cfDNA samples are present
    if not process_functions.check_values_in_column(allClin, "SAMPLE_CLASS", "cfDNA"):
        logger.error(
            "cfDNA samples should not be filtered out in the clinical dataframe."
        )

    allClin.to_csv(clinical_path, sep="\t", index=False)

    gene_matrixdf = gene_matrixdf[gene_matrixdf["SAMPLE_ID"].isin(publicReleaseSamples)]
    gene_matrixdf.to_csv(data_gene_panel_path, sep="\t", index=False)
    load.store_file(
        syn=syn,
        filepath=data_gene_panel_path,
        parentid=public_release_preview,
        version_comment=genie_version,
        name="data_gene_matrix.txt",
    )
    load.store_file(
        syn=syn,
        filepath=clinical_path,
        parentid=public_release_preview,
        version_comment=genie_version,
        name="data_clinical.txt",
    )

    create_case_lists.main(
        clinical_path,
        assay_info.path,
        database_to_staging.CASE_LIST_PATH,
        "genie_public",
    )

    caseListFiles = os.listdir(database_to_staging.CASE_LIST_PATH)
    caseListEntities = []
    for casePath in caseListFiles:
        casePath = os.path.join(database_to_staging.CASE_LIST_PATH, casePath)
        caseListEntities.append(
            load.store_file(
                syn=syn,
                filepath=casePath,
                parentid=public_release_preview_caselist,
                version_comment=genie_version,
            )
        )

    # Grab mapping table to fill in clinical headers
    clinical_code_to_desc_map_synid = databaseSynIdMappingDf["Id"][
        databaseSynIdMappingDf["Database"] == "clinical_code_to_desc_map"
    ][0]
    mapping = extract.get_syntabledf(
        syn=syn, query_string=f"SELECT * FROM {clinical_code_to_desc_map_synid}"
    )
    genePanelEntities = []
    for entName, entId in consortiumRelease[2]:
        is_deprecated_file = entName in ["data_fusions.txt"]
        # skip files to convert
        if (
            entName.startswith("data_linear")
            or "meta_" in entName
            or entName.endswith(".html")
            or entName
            in [
                "data_clinical_sample.txt",
                "data_gene_matrix.txt",
                "data_clinical_patient.txt",
                "data_guide.pdf",
                "release_notes.pdf",
                "samples_to_retract.csv",
                "non_somatic.csv",
                "snv_as_dnp.csv",
                "snv_as_onp.csv",
                "duplicated_variants.csv",
            ]
            or is_deprecated_file
        ):
            # data_gene_matrix was processed above because it had to be
            # used for generating caselists
            continue
        if entName == "data_clinical.txt":
            patientCols = publicRelease["fieldName"][
                publicRelease["level"] == "patient"
            ].tolist()
            sampleCols = ["PATIENT_ID"]
            sampleCols.extend(
                publicRelease["fieldName"][publicRelease["level"] == "sample"].tolist()
            )
            # clinicalDf is defined on line 127
            clinicalDf = clinicalDf[clinicalDf["SAMPLE_ID"].isin(publicReleaseSamples)]

            # Delete columns that are private scope
            # for private in privateRelease:
            #   del clinicalDf[private]
            process_functions.addClinicalHeaders(
                clinicalDf,
                mapping,
                patientCols,
                sampleCols,
                clinical_sample_path,
                clinicl_patient_path,
            )

            load.store_file(
                syn=syn,
                filepath=clinical_sample_path,
                parentid=public_release_preview,
                version_comment=genie_version,
                name="data_clinical_sample.txt",
            )
            load.store_file(
                syn=syn,
                filepath=clinicl_patient_path,
                parentid=public_release_preview,
                version_comment=genie_version,
                name="data_clinical_patient.txt",
            )

        elif "mutation" in entName:
            mutation = syn.get(entId, followLink=True)
            mutationDf = pd.read_csv(mutation.path, sep="\t", comment="#")
            # mutationDf = commonVariantFilter(mutationDf)
            mutationDf["FILTER"] = "PASS"
            mutationDf = mutationDf[
                mutationDf["Tumor_Sample_Barcode"].isin(publicReleaseSamples)
            ]
            text = process_functions.removeFloat(mutationDf)
            with open(mutations_path, "w") as f:
                f.write(text)
            load.store_file(
                syn=syn,
                filepath=mutations_path,
                parentid=public_release_preview,
                version_comment=genie_version,
                name="data_mutations_extended.txt",
            )
        elif "CNA" in entName:
            cna = syn.get(entId, followLink=True)
            cnaDf = pd.read_csv(cna.path, sep="\t")
            cna_columns = pd.concat([publicReleaseSamples, pd.Series("Hugo_Symbol")])
            # parse out the CNA columns to keep
            cnaDf = cnaDf[cnaDf.columns[cnaDf.columns.isin(cna_columns)]]
            text = process_functions.removeFloat(cnaDf)
            text = (
                text.replace("\t\t", "\tNA\t")
                .replace("\t\t", "\tNA\t")
                .replace("\t\n", "\tNA\n")
            )
            with open(cna_path, "w") as cnaFile:
                cnaFile.write(text)
            load.store_file(
                syn=syn,
                filepath=cna_path,
                parentid=public_release_preview,
                version_comment=genie_version,
                name="data_CNA.txt",
            )
        elif entName.endswith(".seg"):
            seg = syn.get(entId, followLink=True)
            segDf = pd.read_csv(seg.path, sep="\t")
            segDf = segDf[segDf["ID"].isin(publicReleaseSamples)]
            text = process_functions.removeFloat(segDf)
            with open(seg_path, "w") as segFile:
                segFile.write(text)
            load.store_file(
                syn=syn,
                filepath=seg_path,
                parentid=public_release_preview,
                version_comment=genie_version,
                name="data_cna_hg19.seg",
            )
        elif entName == "genomic_information.txt":
            bed = syn.get(entId, followLink=True)
            bedDf = pd.read_csv(bed.path, sep="\t")
            bedDf = bedDf[bedDf.SEQ_ASSAY_ID.isin(allClin.SEQ_ASSAY_ID)]
            bedDf.to_csv(combined_bed_path, sep="\t", index=False)
            load.store_file(
                syn=syn,
                filepath=combined_bed_path,
                parentid=public_release_preview,
                version_comment=genie_version,
                name="genomic_information.txt",
            )
        else:
            ent = syn.get(entId, followLink=True, downloadFile=False)
            copiedId = _copyRecursive(
                syn,
                ent,
                public_release_preview,
                version=ent.versionNumber,
                updateExisting=True,
                setProvenance=None,
                skipCopyAnnotations=True,
            )
            copiedEnt = syn.get(copiedId[ent.id], downloadFile=False)
            # Set version comment
            copiedEnt.versionComment = genie_version
            copiedEnt = syn.store(copiedEnt, forceVersion=False)
            # There was a time when gene panel files had to be renamed
            # with an appended genie version. But... GEN-76
            # No longer appending genie verison to release files
            # So just need to track gene panel entities
            if entName.startswith("data_gene_panel"):
                genePanelEntities.append(copiedEnt)

    return caseListEntities, genePanelEntities
