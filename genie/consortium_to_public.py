"""Converts consortium release files to public release files"""

import logging
import os
import shutil

import synapseclient
import synapseutils
import pandas as pd

from . import process_functions
from . import database_to_staging
from . import create_case_lists

logger = logging.getLogger(__name__)


def storeFile(syn, filePath, parentId, genie_version, name=None):
    """Stores file with genie version as comment

    Args:
        syn: Synapse object
        filePath: Path to file
        parentId: Synapse id of folder

    Returns:
        Stored Entity
    """
    if name is None:
        name = os.path.basename(filePath)
    file_ent = synapseclient.File(
        filePath, name=name, parent=parentId, versionComment=genie_version
    )
    file_ent = syn.store(file_ent)
    return file_ent


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
    fusions_path = os.path.join(
        database_to_staging.GENIE_RELEASE_DIR, "data_fusions.txt"
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
    public_release_preview_caselist = database_to_staging.find_caselistid(
        syn, public_release_preview
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
    remove_sc_samples = database_to_staging.sample_class_filter(clinical_df=clinicalDf)
    removeForPublicSamples = list(set(removeForPublicSamples).union(remove_sc_samples))
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
    # If consortium -> Don't release to public
    clinicalReleaseScope = syn.tableQuery(
        "SELECT * FROM syn8545211 where releaseScope = 'public'"
    )
    publicRelease = clinicalReleaseScope.asDataFrame()

    allClin = clinicalDf[clinicalDf["SAMPLE_ID"].isin(publicReleaseSamples)]
    allClin.to_csv(clinical_path, sep="\t", index=False)

    gene_matrixdf = gene_matrixdf[gene_matrixdf["SAMPLE_ID"].isin(publicReleaseSamples)]
    gene_matrixdf.to_csv(data_gene_panel_path, sep="\t", index=False)
    storeFile(
        syn,
        data_gene_panel_path,
        public_release_preview,
        genie_version,
        name="data_gene_matrix.txt",
    )
    storeFile(
        syn,
        clinical_path,
        public_release_preview,
        genie_version,
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
            storeFile(syn, casePath, public_release_preview_caselist, genie_version)
        )

    # Grab mapping table to fill in clinical headers
    mapping_table = syn.tableQuery("SELECT * FROM syn9621600")
    mapping = mapping_table.asDataFrame()
    genePanelEntities = []
    for entName, entId in consortiumRelease[2]:
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
            ]
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

            storeFile(
                syn,
                clinical_sample_path,
                public_release_preview,
                genie_version,
                name="data_clinical_sample.txt",
            )
            storeFile(
                syn,
                clinicl_patient_path,
                public_release_preview,
                genie_version,
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
            storeFile(
                syn,
                mutations_path,
                public_release_preview,
                genie_version,
                name="data_mutations_extended.txt",
            )

        elif "fusion" in entName:
            fusion = syn.get(entId, followLink=True)
            fusionDf = pd.read_csv(fusion.path, sep="\t")
            fusionDf = fusionDf[
                fusionDf["Tumor_Sample_Barcode"].isin(publicReleaseSamples)
            ]
            fusionDf.to_csv(fusions_path, sep="\t", index=False)
            storeFile(
                syn,
                fusions_path,
                public_release_preview,
                genie_version,
                name="data_fusions.txt",
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
            storeFile(
                syn,
                cna_path,
                public_release_preview,
                genie_version,
                name="data_CNA.txt",
            )
        elif entName.endswith(".seg"):
            seg = syn.get(entId, followLink=True)
            segDf = pd.read_csv(seg.path, sep="\t")
            segDf = segDf[segDf["ID"].isin(publicReleaseSamples)]
            text = process_functions.removeFloat(segDf)
            with open(seg_path, "w") as segFile:
                segFile.write(text)
            storeFile(
                syn,
                seg_path,
                public_release_preview,
                genie_version,
                name="data_cna_hg19.seg",
            )
        elif entName == "genomic_information.txt":
            bed = syn.get(entId, followLink=True)
            bedDf = pd.read_csv(bed.path, sep="\t")
            bedDf = bedDf[bedDf.SEQ_ASSAY_ID.isin(allClin.SEQ_ASSAY_ID)]
            bedDf.to_csv(combined_bed_path, sep="\t", index=False)
            storeFile(
                syn,
                combined_bed_path,
                public_release_preview,
                genie_version,
                name="genomic_information.txt",
            )
        # elif entName.startswith("data_gene_panel"):
        #     genePanel = syn.get(entId, followLink=True)
        #     # Create new gene panel naming and store
        #     fileName = os.path.basename(genePanel.path)
        #     # newFileList = fileName.split("_")
        #     # newFileList[-1] = genie_version + ".txt"
        #     # newFileName = "_".join(newFileList)
        #     genePanelPath = os.path.join(
        #         database_to_staging.GENIE_RELEASE_DIR, fileName
        #     )
        #     shutil.copy(genePanel.path, genePanelPath)
        #     # del newFileList[-1]
        #     # entName = "_".join(newFileList)
        #     # entName = entName + ".txt"
        #     genepanel_ent = storeFile(
        #         syn, genePanelPath, public_release_preview, genie_version, name=entName
        #     )
        #     genePanelEntities.append(genepanel_ent)
        else:
            ent = syn.get(entId, followLink=True, downloadFile=False)
            copiedId = synapseutils.copy(
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


def get_public_to_consortium_synid_mapping(
    syn: synapseclient.Synapse, release_synid: str
) -> dict:
    """
    Gets the mapping between potential public release names and
    the consortium release folder

    Args:
        syn (Synapse): Synapse connection
        release_synid (str): Release folder fileview

    Returns:
        dict: Mapping between potential public release and consortium
              release synapse id
    """
    # This dict contains the mapping between public release name and
    # consortium release folder
    public_to_consortium_map = dict()
    # release_files = synapseutils.walk(syn, releaseSynId)
    # TODO: fix the database to mapping table
    consortium_release_folders = syn.tableQuery(
        f"SELECT name, id FROM {release_synid} WHERE "
        "name NOT LIKE 'Release %' "
        "and name NOT LIKE '%-public' "
        "and name NOT IN ('case_lists', 'potential_artifacts')"
        "ORDER BY name"
    )
    consortium_release_folders_df = consortium_release_folders.asDataFrame()
    # Get major release version
    consortium_release_folders_df["major_release"] = [
        release.split(".")[0] for release in consortium_release_folders_df["name"]
    ]
    # only keep the latest consortium release for the public release
    consortium_release_folders_df.drop_duplicates(
        "major_release", keep="last", inplace=True
    )

    for _, release_info in consortium_release_folders_df.iterrows():
        major_release = release_info["major_release"]
        # add support for potential patch releases
        for num in [0, 1, 2, 3]:
            # This has to exist because the the first three GENIE releases
            # used semantic versioning
            if release_info["major_release"] in ["0", "1", "2"]:
                public_release_name = f"{int(major_release) + 1}.{num}.0"
                public_to_consortium_map[public_release_name] = release_info["id"]
            else:
                public_release_name = f"{major_release}.{num}-public"
                public_to_consortium_map[public_release_name] = release_info["id"]
    return public_to_consortium_map
