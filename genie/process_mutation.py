"""Process mutation files"""
import logging
import os
import shutil
import subprocess
import tempfile

import pandas as pd
import synapseclient  # lgtm [py/import-and-import-from]
from synapseclient import Synapse
from synapseclient.core.exceptions import SynapseTimeoutError

from . import process_functions

logger = logging.getLogger(__name__)

# Some columns are already capitalized, so they aren't included here
MAF_COL_MAPPING = {
    "HUGO_SYMBOL": "Hugo_Symbol",
    "ENTREZ_GENE_ID": "Entrez_Gene_Id",
    "CENTER": "Center",
    "NCBI_BUILD": "NCBI_Build",
    "CHROMOSOME": "Chromosome",
    "START_POSITION": "Start_Position",
    "END_POSITION": "End_Position",
    "STRAND": "Strand",
    "VARIANT_CLASSIFICATION": "Variant_Classification",
    "VARIANT_TYPE": "Variant_Type",
    "REFERENCE_ALLELE": "Reference_Allele",
    "TUMOR_SEQ_ALLELE1": "Tumor_Seq_Allele1",
    "TUMOR_SEQ_ALLELE2": "Tumor_Seq_Allele2",
    "DBSNP_RS": "dbSNP_RS",
    "DBSNP_VAL_STATUS": "dbSNP_Val_Status",
    "TUMOR_SAMPLE_BARCODE": "Tumor_Sample_Barcode",
    "MATCHED_NORM_SAMPLE_BARCODE": "Matched_Norm_Sample_Barcode",
    "MATCH_NORM_SEQ_ALLELE1": "Match_Norm_Seq_Allele1",
    "MATCH_NORM_SEQ_ALLELE2": "Match_Norm_Seq_Allele2",
    "TUMOR_VALIDATION_ALLELE1": "Tumor_Validation_Allele1",
    "TUMOR_VALIDATION_ALLELE2": "Tumor_Validation_Allele2",
    "MATCH_NORM_VALIDATION_ALLELE1": "Match_Norm_Validation_Allele1",
    "MATCH_NORM_VALIDATION_ALLELE2": "Match_Norm_Validation_Allele2",
    "VERIFICATION_STATUS": "Verification_Status",
    "VALIDATION_STATUS": "Validation_Status",
    "MUTATION_STATUS": "Mutation_Status",
    "SEQUENCING_PHASE": "Sequencing_Phase",
    "SEQUENCE_SOURCE": "Sequence_Source",
    "VALIDATION_METHOD": "Validation_Method",
    "SCORE": "Score",
    "BAM_FILE": "BAM_File",
    "SEQUENCER": "Sequencer",
    "T_REF_COUNT": "t_ref_count",
    "T_ALT_COUNT": "t_alt_count",
    "N_REF_COUNT": "n_ref_count",
    "N_ALT_COUNT": "n_alt_count",
    "ALLELE": "Allele",
    "AMINO_ACID_CHANGE": "amino_acid_change",
    "AMINO_ACIDS": "Amino_acids",
    "CDS_POSITION": "CDS_position",
    "CODONS": "Codons",
    "CONSEQUENCE": "Consequence",
    "EXISTING_VARIATION": "Existing_variation",
    "EXON_NUMBER": "Exon_Number",
    "FEATURE": "Feature",
    "FEATURE_TYPE": "Feature_type",
    "GENE": "Gene",
    "HGVSC": "HGVSc",
    "HGVSP": "HGVSp",
    "HGVSP_SHORT": "HGVSp_Short",
    "HOTSPOT": "Hotspot",
    "MA:FIMPACT": "MA:FImpact",
    "MA:LINK.MSA": "MA:link.MSA",
    "MA:LINK.PDB": "MA:link.PDB",
    "MA:LINK.VAR": "MA:link.var",
    "MA:PROTEIN.CHANGE": "MA:protein.change",
    "POLYPHEN": "PolyPhen",
    "PROTEIN_POSITION": "Protein_position",
    "REFSEQ": "RefSeq",
    "TRANSCRIPT": "transcript",
    "TRANSCRIPT_ID": "Transcript_ID",
    "ALL_EFFECTS": "all_effects",
    "CDNA_CHANGE": "cdna_change",
    "CDNA_POSITION": "cDNA_position",
    "N_DEPTH": "n_depth",
    "T_DEPTH": "t_depth",
}

KNOWN_STRING_COLS = [
    "IS_NEW",
    "ALLELE_NUM",
    "Chromosome",
    "CLIN_SIG",
    "MOTIF_NAME",
    "HIGH_INF_POS",
    "MINIMISED",
    "CHROMOSOME",
    "VERIFICATION_STATUS",
    "VALIDATION_STATUS",
    "MUTATION_STATUS",
    "SEQUENCE_SOURCE",
    "SEQUENCER",
    "REPORT_AF",
    "CDNA_CHANGE",
    "AMINO_ACID_CHANGE",
    "TRANSCRIPT",
    "STRAND_VEP",
    "HGNC_ID",
    "PUBMED",
    "PICK",
    "Exon_Number",
]


def _convert_to_str_dtype(column_types, known_string_cols):
    """Sometimes the deteremined dtype is incorrect based off the first
    100 rows, update the incorrect dtypes.
    """
    for str_col in known_string_cols:
        if column_types.get(str_col):
            column_types[str_col] = "object"
    return column_types


def determine_dtype(path: str):
    """Reads in a dataframe partially and determines the dtype of columns"""
    # Change this nrows to 5000 so that it better encapsulates the types
    subset_df = pd.read_csv(path, nrows=5000, sep="\t", comment="#")
    column_types = subset_df.dtypes.to_dict()
    return column_types


def move_and_configure_maf(mutation_path: str, input_files_dir: str) -> str:
    """Moves maf files into processing directory. Maf file's column headers
    are renamed if necessary and .0 are stripped.

    Args:
        mutation_path (str): Mutation file path
        input_files_dir (str): Input file directory

    Returns:
        str: Filepath to moved and configured maf
    """
    filename = os.path.basename(mutation_path)
    new_filepath = os.path.join(input_files_dir, filename)
    column_types = determine_dtype(mutation_path)
    new_column_types = _convert_to_str_dtype(column_types, KNOWN_STRING_COLS)
    mafdf = pd.read_csv(mutation_path, sep="\t", dtype=new_column_types, comment="#")
    # If any column headers need to be remapped, remap
    mafdf = mafdf.rename(columns=MAF_COL_MAPPING)
    # Must remove floating .0 or else processing will fail for genome nexus
    maf_text = process_functions.removePandasDfFloat(mafdf)
    with open(new_filepath, "w") as new_maf_f:
        new_maf_f.write(maf_text)
    return new_filepath


def move_mutation(mutation_path, input_files_dir):
    """Move mutation file into processing directory"""
    # If mutation file is vcf, just copy
    if mutation_path.endswith(".vcf"):
        shutil.copy(mutation_path, input_files_dir)
    else:
        move_and_configure_maf(mutation_path, input_files_dir)


def process_mutation_workflow(
    syn: Synapse,
    center: str,
    validfiles: pd.DataFrame,
    genie_config: dict,
    workdir: str,
) -> str:
    """Process vcf/maf workflow

    Args:
        syn: Synapse connection
        center: Center name
        validfiles: Center validated files
        genie_config: GENIE configuration.
        workdir: Working directory

    Returns:
        Annotated Maf Path

    """
    # Get valid files
    mutation_files = validfiles["fileType"].isin(["maf", "vcf"])
    valid_mutation_files = validfiles["path"][mutation_files].tolist()
    # If there are no valid mutation files, return
    if not valid_mutation_files:
        logger.info("No mutation data")
        return
    # Certificate to use GENIE Genome Nexus
    syn.get(
        "syn22053204",
        ifcollision="overwrite.local",
        downloadLocation=genie_config["genie_annotation_pkg"],
    )
    # Genome Nexus Jar file
    syn.get(
        "syn22084320",
        ifcollision="overwrite.local",
        downloadLocation=genie_config["genie_annotation_pkg"],
    )

    annotated_maf_path = annotate_mutation(
        center=center,
        mutation_files=valid_mutation_files,
        genie_annotation_pkg=genie_config["genie_annotation_pkg"],
        workdir=workdir,
    )

    maf_tableid = genie_config["vcf2maf"]
    flatfiles_synid = genie_config["centerMaf"]
    # Split into narrow maf and store into db / flat file
    split_and_store_maf(
        syn=syn,
        center=center,
        maf_tableid=maf_tableid,
        annotated_maf_path=annotated_maf_path,
        flatfiles_synid=flatfiles_synid,
        workdir=workdir,
    )

    return annotated_maf_path


def annotate_mutation(
    center: str, mutation_files: list, genie_annotation_pkg: str, workdir: str
) -> str:
    """Process vcf/maf files

    Args:
        center: Center name
        mutation_files: list of mutation files
        genie_annotation_pkg: Path to GENIE annotation package

    Returns:
        Path to final maf
    """
    input_files_dir = tempfile.mkdtemp(dir=workdir)
    output_files_dir = tempfile.mkdtemp(dir=workdir)

    for mutation_file in mutation_files:
        move_mutation(mutation_file, input_files_dir)

    merged_maf_path = os.path.join(
        output_files_dir, f"data_mutations_extended_{center}.txt"
    )
    annotater_cmd = [
        "bash",
        os.path.join(genie_annotation_pkg, "annotation_suite_wrapper.sh"),
        f"-i={input_files_dir}",
        f"-o={output_files_dir}",
        f"-m={merged_maf_path}",
        f"-c={center}",
        "-s=WXS",
        f"-p={genie_annotation_pkg}",
    ]

    subprocess.check_call(annotater_cmd)

    return merged_maf_path


def append_or_createdf(dataframe: pd.DataFrame, filepath: str):
    """Creates a file with the dataframe or appends to a existing file.

    Args:
        df: pandas.dataframe to write out
        filepath: Filepath to append or create

    """
    if not os.path.exists(filepath) or os.stat(filepath).st_size == 0:
        dataframe.to_csv(filepath, sep="\t", index=False)
    else:
        dataframe.to_csv(filepath, sep="\t", mode="a", index=False, header=None)


def store_full_maf(syn: Synapse, filepath: str, parentid: str):
    """Stores full maf file

    Args:
        syn: Synapse connection
        filepath: Path to file
        parentid: Synapse container id

    """
    syn.store(synapseclient.File(filepath, parentId=parentid))


def store_narrow_maf(syn: Synapse, filepath: str, maf_tableid: str):
    """
    Stores the narrow maf in Synapse Table

    Args:
        syn: Synapse connection
        filepath: Path to maf file
        maf_tableid: database synid

    """
    logger.info(f"STORING {filepath}")
    # database = syn.get(maf_tableid)
    try:
        update_table = synapseclient.Table(maf_tableid, filepath, separator="\t")
        syn.store(update_table)
    except SynapseTimeoutError:
        # This error occurs because of waiting for table to index.
        # Don't worry about this.
        pass


def format_maf(mafdf: pd.DataFrame, center: str) -> pd.DataFrame:
    """Format maf file, shortens the maf file length

    Args:
        mafdf: mutation dataframe
        center: Center name

    Returns:
        Formatted mutation dataframe"""
    mafdf["Center"] = center
    # Leaving here for safe guarding.
    mafdf["Tumor_Sample_Barcode"] = [
        process_functions.checkGenieId(i, center) for i in mafdf["Tumor_Sample_Barcode"]
    ]
    mafdf["Sequence_Source"] = float("nan")
    mafdf["Sequencer"] = float("nan")
    mafdf["Validation_Status"][
        mafdf["Validation_Status"].isin(["Unknown", "unknown"])
    ] = ""

    return mafdf


def split_and_store_maf(
    syn: "Synapse",
    center: str,
    maf_tableid: str,
    annotated_maf_path: str,
    flatfiles_synid: str,
    workdir: str,
):
    """Separates annotated maf file into narrow and full maf and stores them

    Args:
        syn: Synapse connection
        center: Center
        maf_tableid: Mutation table synapse id
        annotated_maf_path: Annotated maf
        flatfiles_synid: GENIE flat files folder

    """
    narrow_maf_cols = [
        col["name"]
        for col in syn.getTableColumns(maf_tableid)
        if col["name"] != "inBED"
    ]
    full_maf_path = os.path.join(
        workdir, center, "staging", f"data_mutations_extended_{center}.txt"
    )
    narrow_maf_path = os.path.join(
        workdir, center, "staging", f"data_mutations_extended_{center}_MAF_narrow.txt"
    )
    maf_chunks = pd.read_csv(
        annotated_maf_path, sep="\t", chunksize=100000, comment="#"
    )

    for maf_chunk in maf_chunks:
        maf_chunk = format_maf(maf_chunk, center)
        append_or_createdf(maf_chunk, full_maf_path)
        narrow_maf_chunk = maf_chunk[narrow_maf_cols]
        append_or_createdf(narrow_maf_chunk, narrow_maf_path)

    store_narrow_maf(syn, narrow_maf_path, maf_tableid)
    # Store MAF flat file into synapse
    store_full_maf(syn, full_maf_path, flatfiles_synid)
