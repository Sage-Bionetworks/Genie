"""Process mutation files"""
import os
import shutil
import subprocess
import tempfile

import pandas as pd
import synapseclient
try:
    from synapesclient.exceptions import SynapseTimeoutError
except ModuleNotFoundError:
    from synapseclient.core.exceptions import SynapseTimeoutError

from . import process_functions

WORKDIR = os.path.expanduser("~/.synapseCache")


def process_mutation_workflow(syn, center, mutation_files,
                              genie_annotation_pkg,
                              maf_tableid, flatfiles_synid):
    """Process vcf/maf workflow"""
    annotated_maf_path = annotate_mutation(
        center=center,
        mutation_files=mutation_files,
        genie_annotation_pkg=genie_annotation_pkg)

    # Split into narrow maf and store into db / flat file
    split_and_store_maf(syn=syn,
                        center=center,
                        maf_tableid=maf_tableid,
                        annotated_maf_path=annotated_maf_path,
                        flatfiles_synid=flatfiles_synid)

    return annotated_maf_path


def annotate_mutation(center: str, mutation_files: list,
                      genie_annotation_pkg: str) -> str:
    """Process vcf/maf files

    Args:
        center: Center name
        mutation_files: list of mutation files
        genie_annotation_pkg: Path to GENIE annotation package

    Returns:
        Path to final maf
    """
    input_files_dir = tempfile.mkdtemp(dir=WORKDIR)
    output_files_dir = tempfile.mkdtemp(dir=WORKDIR)

    for mutation_file in mutation_files:
        shutil.copyfile(mutation_file, input_files_dir)

    annotater_cmd = ['bash', os.path.join(genie_annotation_pkg,
                                          'annotation_suite_wrapper.sh'),
                     f'-i={input_files_dir}',
                     f'-o={output_files_dir}',
                     f'-m=data_mutations_extended_{center}.txt',
                     f'-c={center}',
                     '-s=WXS',
                     f'-p={genie_annotation_pkg}']

    subprocess.check_call(annotater_cmd)

    return os.path.join(output_files_dir,
                        "annotated", f"data_mutations_extended_{center}.txt")


def append_or_createdf(dataframe: 'DataFrame', filepath: str):
    """Creates a file with the dataframe or appends to a existing file.

    Args:
        df: pandas.dataframe to write out
        filepath: Filepath to append or create

    """
    if os.stat(filepath).st_size == 0:
        dataframe.to_csv(filepath, sep="\t", index=False)
    else:
        dataframe.to_csv(filepath, sep="\t", mode='a', index=False,
                         header=None)
    # write_or_append = "wb" if maf else "ab"
    # with open(filepath, write_or_append) as maf_file:
    #     maf_text = process_functions.removeStringFloat(maf_text)
    #     maf_file.write(maf_text.encode("utf-8"))


def store_full_maf(syn: 'Synapse', filepath: str, parentid: str):
    """Stores full maf file"""
    syn.store(synapseclient.File(filepath, parentId=parentid))


def store_narrow_maf(syn: 'Synapse', filepath: str, maf_tableid: str):
    '''
    Stores the processed maf
    There is a isNarrow option, but note that the number of rows
    of the maf file DOES NOT change in this function

    Args:
        filePath: Path to maf file
        mafSynId: database synid
        centerMafSynid: center flat file folder synid
        isNarrow: Is the file a narrow maf. Defaul to False.
    '''
    logger.info('STORING %s' % filepath)
    database = syn.get(maf_tableid)
    try:
        update_table = synapseclient.Table(database.id, filepath,
                                           separator="\t")
        syn.store(update_table)
    except SynapseTimeoutError:
        # This error occurs because of waiting for table to index.
        # Don't worry about this.
        pass


def format_maf(mafdf: 'DataFrame', center: str) -> 'DataFrame':
    """Format maf file, shortens the maf file length"""
    mafdf['Center'] = center
    mafdf['Tumor_Sample_Barcode'] = [
        process_functions.checkGenieId(i, center)
        for i in mafdf['Tumor_Sample_Barcode']
    ]

    mafdf['Sequence_Source'] = float('nan')
    mafdf['Sequencer'] = float('nan')
    mafdf['Validation_Status'][
        mafdf['Validation_Status'].isin(["Unknown", "unknown"])
    ] = ''

    return mafdf


def split_and_store_maf(syn: 'Synapse', center: str, maf_tableid: str,
                        annotated_maf_path: str, flatfiles_synid: str):
    """Separates annotated maf file into narrow and full maf and stores them

    Args:
        syn: Synapse connection
        center: Center
        maf_tableid: Mutation table synapse id
        annotated_maf_path: Annotated maf
        flatfiles_synid: GENIE flat files folder

    """
    narrow_maf_cols = [col['name']
                       for col in syn.getTableColumns(maf_tableid)
                       if col['name'] != 'inBED']
    full_maf_path = os.path.join(
        WORKDIR, center, "staging",
        f"data_mutations_extended_{center}_MAF.txt"
    )
    narrow_maf_path = os.path.join(
        WORKDIR, center, "staging",
        f"data_mutations_extended_{center}_MAF_narrow.txt"
    )
    maf_chunks = pd.read_csv(annotated_maf_path, chunksize=100000)

    for maf_chunk in maf_chunks:
        maf_chunk = format_maf(maf_chunk, center)
        append_or_createdf(maf_chunk, full_maf_path)
        narrow_maf_chunk = maf_chunk[narrow_maf_cols]
        append_or_createdf(narrow_maf_chunk, narrow_maf_path)

    store_narrow_maf(syn, narrow_maf_path, maf_tableid)
    # Store MAF flat file into synapse
    store_full_maf(syn, full_maf_path, flatfiles_synid)
