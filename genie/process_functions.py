"""Processing functions that are used in the GENIE pipeline"""
import datetime
import json
import logging
import os
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import tempfile
import time
from typing import List

import ast
from Crypto.Cipher import PKCS1_OAEP
from Crypto.PublicKey import RSA
import pandas as pd
import synapseclient  # lgtm [py/import-and-import-from]
from synapseclient import Synapse

# try:
#   from urllib.request import urlopen
# except ImportError:
#   from urllib2 import urlopen
# Ignore SettingWithCopyWarning warning
pd.options.mode.chained_assignment = None

logger = logging.getLogger(__name__)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# try:
#   from urlparse import urlparse
# except ImportError:
#   from urllib.parse import urlparse

# Create merged dictionary of remapped genes
# def createDict(invalidated_genes):
#     toRemapRemove = {}
#     invalidated_genes = pd.Series(invalidated_genes)
#     for i in invalidated_genes[invalidated_genes!=True]:
#         toRemapRemove.update(i)
#     return(toRemapRemove)

# Validate genes
# def hgncRestCall(path):
#     """
#     This function does the rest call to the genenames website

#     :params path:     The gene symbol url path to add to the base uri

#     :returns:         If the symbol exists, returns True and the corrected symbol, otherwise returns False and None.
#     """
#     headers = {'Accept': 'application/json',}

#     uri = 'http://rest.genenames.org'

#     target = urlparse(uri+path)
#     method = 'GET'
#     body = ''
#     h = http.Http()
#     response, content = h.request(target.geturl(),
#                                   method,
#                                   body,
#                                   headers)
#     if response['status'] == '200':
#         data = json.loads(content)
#         if len(data['response']['docs']) == 0:
#             return(False, [None])
#         else:
#             mapped = [symbol['symbol'] for symbol in data['response']['docs']]
#             return(True, mapped)
#     else:
#         return(False, [None])

# Validation of gene names
# def validateSymbol(gene, returnMapping=False):
#     """
#     This function does validation of symbols

#     :params gene:               Gene symbol
#     :params returnMapping:      Return mapping of old gene to new gene

#     :returns:                   Check if the provided gene name is a correct symbol and print out genes
#                                 that need to be remapped or can't be mapped to anything
#     """
#     path = '/fetch/symbol/%s' %  gene
#     verified, symbol = hgncRestCall(path)
#     if not verified:
#         path = '/fetch/prev_symbol/%s' %  gene
#         verified, symbol = hgncRestCall(path)
#     if not verified:
#         path = '/fetch/alias_symbol/%s' %  gene
#         verified, symbol = hgncRestCall(path)
#     if gene in symbol:
#         return(True)
#     else:
#         if symbol[0] is None:
#             #logger.error("%s cannot be remapped. Please correct." % gene)
#             logger.warning("%s cannot be remapped. These rows will have an empty gene symbol" % gene)
#         else:
#             #if "MLL4", then the HUGO symbol should be KMT2D and KMT2B
#             if len(symbol) > 1:
#                 #logger.error("%s can be mapped to different symbols: %s. Please correct." % (gene, ", ".join(symbol)))
#                 logger.warning("%s can be mapped to different symbols: %s. Please correct or it will be removed." % (gene, ", ".join(symbol)))
#             else:
#                 logger.info("%s will be remapped to %s" % (gene, symbol[0]))
#                 if returnMapping:
#                     return({gene: symbol[0]})
#                 else:
#                     return(True)
#         if returnMapping:
#             return({gene: float('nan')})
#         else:
#             return(False)
#         #return(True)

# # Remap and remove genes in dataframes given a column
# def remapGenes(invalidated_genes, DF, col,isBedFile=False):
#     nonmapped = []
#     gene_dict = createDict(invalidated_genes)
#     for key in gene_dict:
#         value = gene_dict[key]
#         if isBedFile:
#         #     nonmapped.append(key)
#         #     DF = DF[DF[col] != key]
#         # else:
#             DF[col][DF[col] == key] = value
#     return(DF, nonmapped)


def retry_get_url(url):
    """
    Implement retry logic when getting urls.
    Timesout at 3 seconds, retries 5 times.

    Args:
        url:  Http or https url

    Returns:
        requests.get()
    """
    s = requests.Session()
    retries = Retry(total=5, backoff_factor=1)
    s.mount("http://", HTTPAdapter(max_retries=retries))
    s.mount("https://", HTTPAdapter(max_retries=retries))
    response = s.get(url, timeout=3)
    return response


def checkUrl(url):
    """
    Check if URL link is live

    Args:
        url: web URL
    """
    temp = retry_get_url(url)
    assert temp.status_code == 200, "%s site is down" % url


def checkColExist(DF, key):
    """
    This function checks if the column exists in a dataframe

    Args:
        DF: pandas dataframe
        key: Expected column header name

    Returns:
        bool:  True if column exists
    """
    result = False if DF.get(key) is None else True
    return result


def validate_genie_identifier(
    identifiers: pd.Series, center: str, filename: str, col: str
) -> str:
    """Validate GENIE sample and patient ids.

    Args:
        identifiers (pd.Series): Array of GENIE identifiers
        center (str): GENIE center name
        filename (str): name of file
        col (str): Column with identifiers

    return:
        str: Errors
    """
    total_error = ""
    if not all(identifiers.str.startswith(f"GENIE-{center}")):
        total_error = total_error + (
            f"{filename}: {col} must start with GENIE-{center}\n"
        )
    if any(identifiers.str.len() >= 50):
        total_error = total_error + (
            f"{filename}: {col} must have less than 50 characters.\n"
        )
    return total_error


# def get_oncotree_codes(oncotree_url):
#     '''
#     Deprecated
#     Gets the oncotree data from the specified url.

#     Args:
#         oncotree_url: Oncotree url

#     Returns:
#         dataframe: oncotree codes
#     '''
#     # PATTERN = re.compile('([A-Za-z\' ,-/]*) \\(([A-Za-z_]*)\\)')
#     PATTERN = re.compile('.*[(](.*)[)]')
#     # with requests.get(oncotree_url) as oncotreeUrl:
#     #   oncotree = oncotreeUrl.text.split("\n")
#     oncotreeUrl = retry_get_url(oncotree_url)
#     oncotree = oncotreeUrl.text.split("\n")

#     # oncotree = urlopen(oncotree_url).read().split('\n')
#     allCodes = []
#     for row in oncotree[:-1]:
#         data = row.split("\t")
#         allCodes.extend([
#             PATTERN.match(i).group(1)
#             for i in data[0:5] if PATTERN.match(i)])
#     codes = pd.DataFrame({"ONCOTREE_CODE": list(set(allCodes))})
#     return(codes)


def lookup_dataframe_value(df, col, query):
    """
    Look up dataframe value given query and column

    Args:
        df: dataframe
        col: column with value to return
        query: Query for specific column

    Returns:
        value
    """
    query = df.query(query)
    query_val = query[col].iloc[0]
    return query_val


def get_syntabledf(syn, query_string):
    """
    Get dataframe from table query

    Args:
        syn: Synapse object
        query_string: Table query

    Returns:
        pandas dataframe with query results
    """
    table = syn.tableQuery(query_string)
    tabledf = table.asDataFrame()
    return tabledf


def get_synid_database_mappingdf(syn, project_id):
    """
    Get database to synapse id mapping dataframe

    Args:
        syn: Synapse object
        project_id: Synapse Project ID with a 'dbMapping' annotation.

    Returns:
        database to synapse id mapping dataframe
    """

    project = syn.get(project_id)
    database_mapping_synid = project.annotations["dbMapping"][0]
    database_map_query = "SELECT * FROM {}".format(database_mapping_synid)
    mappingdf = get_syntabledf(syn, database_map_query)
    return mappingdf


def getDatabaseSynId(syn, tableName, project_id=None, databaseToSynIdMappingDf=None):
    """
    Get database synapse id from database to synapse id mapping table

    Args:
        syn: Synapse object
        project_id: Synapse Project ID with a database mapping table.
        tableName: Name of synapse table
        databaseToSynIdMappingDf: Avoid calling rest call to download table
                                  if the mapping table is already downloaded

    Returns:
        str:  Synapse id of wanted database
    """
    if databaseToSynIdMappingDf is None:
        databaseToSynIdMappingDf = get_synid_database_mappingdf(
            syn, project_id=project_id
        )

    synId = lookup_dataframe_value(
        databaseToSynIdMappingDf, "Id", 'Database == "{}"'.format(tableName)
    )
    return synId


def rmFiles(folderPath, recursive=True):
    """
    Convenience function to remove all files in dir

    Args:
        folderPath: Path to folder
        recursive:  Removes all files recursively
    """
    for dirPath, dirNames, filePaths in os.walk(folderPath):
        for filePath in filePaths:
            os.unlink(os.path.join(dirPath, filePath))
        if not recursive:
            break


def removeStringFloat(string):
    """
    remove string float in tsv file

    Args:
        string: tsv file in string format

    Return:
        string: string with float removed
    """
    string = string.replace(".0\t", "\t")
    string = string.replace(".0\n", "\n")
    return string


def removePandasDfFloat(df, header=True):
    """
    Remove decimal for integers due to pandas

    Args:
        df:  Pandas dataframe

    Return:
        str: tsv in text
    """
    if header:
        text = df.to_csv(sep="\t", index=False)
    else:
        text = df.to_csv(sep="\t", index=False, header=None)

    text = removeStringFloat(text)
    return text


def removeFloat(df):
    """
    Need to remove this function
    as it calls another function
    """
    # text = df.to_csv(sep="\t",index=False)
    # text = text.replace(".0\t","\t")
    # text = text.replace(".0\n","\n")
    text = removePandasDfFloat(df)
    return text


def checkGenieId(ID, center):
    """
    Checks if GENIE ID is labelled correctly
    and reformats the GENIE ID

    Args:
        ID: string
        center: GENIE center

    Return:
        str: Formatted GENIE ID string
    """
    if str(ID).startswith("%s-" % center):
        return "GENIE-%s" % str(ID)
    elif not str(ID).startswith("GENIE-%s-" % center):
        return "GENIE-%s-%s" % (center, str(ID))
    else:
        return str(ID)


def storeFile(
    syn,
    fileName,
    parentId,
    center,
    fileFormat,
    dataSubType,
    platform=None,
    cBioFileFormat=None,
    used=None,
):
    """
    # Storing Files along with annotations
    """
    logger.info("STORING FILES")
    fileEnt = synapseclient.File(fileName, parent=parentId)
    fileEnt.center = center
    fileEnt.species = "Human"
    fileEnt.consortium = "GENIE"
    fileEnt.dataType = "genomicVariants"
    fileEnt.fundingAgency = "AACR"
    fileEnt.assay = "targetGeneSeq"
    fileEnt.fileFormat = fileFormat
    fileEnt.dataSubType = dataSubType
    fileEnt.fileStage = "staging"
    fileEnt.platform = platform
    if platform is not None:
        fileEnt.platform = platform
    if cBioFileFormat is not None:
        fileEnt.cBioFileFormat = cBioFileFormat
    ent = syn.store(fileEnt, used=used)
    return ent


def seqDateFilter(clinicalDf, processingDate, days):
    """
    SEQ_DATE filter
    SEQ_DATE - Clinical data (6 and 12 as parameters)
    Jan-2017 , given processing date (today) ->
        staging release (processing date - Jan-2017 < 6 months)
    July-2016 , given processing date (today) ->
        consortium release (processing date - July-2016 between
        6 months - 12 months)

    """
    copyClinicalDf = clinicalDf.copy()
    # copyClinicalDf['SEQ_DATE'][copyClinicalDf['SEQ_DATE'].astype(str) == '999'] = "Jan-1988"
    # copyClinicalDf['SEQ_DATE'][copyClinicalDf['SEQ_DATE'].astype(str) == '999.0'] = "Jan-1988"
    if not isinstance(processingDate, datetime.datetime):
        processingDate = datetime.datetime.strptime(processingDate, "%b-%Y")
    # Remove this null statement after clinical files have been re-validated
    # copyClinicalDf['SEQ_DATE'][copyClinicalDf['SEQ_DATE'].isnull()] = "Jan-1900"
    copyClinicalDf["SEQ_DATE"][copyClinicalDf["SEQ_DATE"] == "Release"] = "Jan-1900"
    # clinicalDf['SEQ_DATE'][clinicalDf['SEQ_DATE'] == '999'] = "Jan-1988"
    dates = copyClinicalDf["SEQ_DATE"].apply(
        lambda date: datetime.datetime.strptime(date, "%b-%Y")
    )
    keep = processingDate - dates > datetime.timedelta(days)
    keepSamples = copyClinicalDf["SAMPLE_ID"][~keep]
    # copyClinicalDf.SEQ_DATE[keep].unique()
    return keepSamples


def addClinicalHeaders(
    clinicalDf, mapping, patientCols, sampleCols, samplePath, patientPath
):
    """
    Add clinical file headers

    Args:
        clinicalDf: clinical dataframe
        mapping: mapping dataframe, maps clinical columns to
                 labels and descriptions
        patientCols: list of patient columns
        sampleCols: list of sample columns
        samplePath: clinical sample path
        patientPath: clinical patient path
    """
    patientLabels = [
        str(mapping["labels"][mapping["cbio"] == i].values[0]) for i in patientCols
    ]
    sampleLabels = [
        str(mapping["labels"][mapping["cbio"] == i].values[0]) for i in sampleCols
    ]
    patientDesc = [
        str(mapping["description"][mapping["cbio"] == i].values[0]) for i in patientCols
    ]
    sampleDesc = [
        str(mapping["description"][mapping["cbio"] == i].values[0]) for i in sampleCols
    ]
    patientType = [
        str(mapping["colType"][mapping["cbio"] == i].values[0]) for i in patientCols
    ]
    sampleType = [
        str(mapping["colType"][mapping["cbio"] == i].values[0]) for i in sampleCols
    ]

    with open(patientPath, "w+") as patientFile:
        patientFile.write("#%s\n" % "\t".join(patientLabels))
        patientFile.write("#%s\n" % "\t".join(patientDesc))
        patientFile.write("#%s\n" % "\t".join(patientType))
        patientFile.write("#%s\n" % "\t".join(["1"] * len(patientLabels)))
        text = removeFloat(clinicalDf[patientCols].drop_duplicates("PATIENT_ID"))
        patientFile.write(text)
    with open(samplePath, "w+") as sampleFile:
        sampleFile.write("#%s\n" % "\t".join(sampleLabels))
        sampleFile.write("#%s\n" % "\t".join(sampleDesc))
        sampleFile.write("#%s\n" % "\t".join(sampleType))
        sampleFile.write("#%s\n" % "\t".join(["1"] * len(sampleLabels)))
        text = removeFloat(clinicalDf[sampleCols].drop_duplicates("SAMPLE_ID"))
        sampleFile.write(text)


########################################################################
# CENTER ANONYMIZING
########################################################################


# def center_anon(filePath, anonymizeCenterDf):
#     '''
#     Deprecated
#     '''
#     with open(filePath, "r") as datafile:
#         text = datafile.read()
#     for center in anonymizeCenterDf['center']:
#         newCenter = anonymizeCenterDf['newCenter'][
#             anonymizeCenterDf['center'] == center].values[0]
#         text = re.sub("\t%s\t" % center, "\t%s\t" % newCenter, text)
#         text = re.sub("GENIE-%s-" % center, "GENIE-%s-" % newCenter, text)
#     with open(filePath, "w") as datafile:
#         datafile.write(text)


# def center_convert_back(filePath, anonymizeCenterDf):
#     '''
#     Deprecated
#     '''
#     with open(filePath, "r") as datafile:
#         text = datafile.read()
#     for center in anonymizeCenterDf['center']:
#         newCenter = anonymizeCenterDf['newCenter'][
#             anonymizeCenterDf['center'] == center].values[0]
#         text = re.sub("\t%s\t" % newCenter, "\t%s\t" % center, text)
#         text = re.sub("GENIE-%s-" % newCenter, "GENIE-%s-" % center, text)
#     with open(filePath, "w") as datafile:
#         datafile.write(text)

##############################################################################
# UPDATING DATABASE
##############################################################################


def _check_valid_df(df, col):
    """
    Checking if variable is a pandas dataframe and column specified exist

    Args:
        df: Pandas dataframe
        col: Column name
    """
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Must pass in pandas dataframe")
    if df.get(col) is None:
        raise ValueError("'{}' column must exist in dataframe".format(col))


def _get_left_diff_df(left, right, checkby):
    """
    Subset the dataframe based on 'checkby' by taking values in the left df
    that arent in the right df

    Args:
        left: Dataframe
        right: Dataframe
        checkby: Column of values to compare

    Return:
        Dataframe: Subset of dataframe from left that don't exist in the right
    """
    _check_valid_df(left, checkby)
    _check_valid_df(right, checkby)
    diffdf = left[~left[checkby].isin(right[checkby])]
    return diffdf


def _get_left_union_df(left, right, checkby):
    """
    Subset the dataframe based on 'checkby' by taking the union of
    values in the left df with the right df

    Args:
        left: Dataframe
        right: Dataframe
        checkby: Column of values to compare

    Return:
        Dataframe: Subset of dataframe from left that also exist in the right
    """
    _check_valid_df(left, checkby)
    _check_valid_df(right, checkby)
    uniondf = left[left[checkby].isin(right[checkby])]
    return uniondf


def _append_rows(new_datasetdf, databasedf, checkby):
    """
    Compares the dataset from the database and determines which rows to
    append from the dataset

    Args:
        new_datasetdf: Input data dataframe
        databasedf: Existing data dataframe
        checkby: Column of values to compare

    Return:
        Dataframe: Dataframe of rows to append
    """
    databasedf.fillna("", inplace=True)
    new_datasetdf.fillna("", inplace=True)

    appenddf = _get_left_diff_df(new_datasetdf, databasedf, checkby)
    if not appenddf.empty:
        logger.info("Adding Rows")
    else:
        logger.info("No new rows")
    del appenddf[checkby]
    appenddf.reset_index(drop=True, inplace=True)
    return appenddf


def _delete_rows(new_datasetdf, databasedf, checkby):
    """
    Compares the dataset from the database and determines which rows to
    delete from the dataset

    Args:
        new_datasetdf: Input data dataframe
        databasedf: Existing data dataframe
        checkby: Column of values to compare

    Return:
        Dataframe: Dataframe of rows to delete
    """

    databasedf.fillna("", inplace=True)
    new_datasetdf.fillna("", inplace=True)
    # If the new dataset is empty, delete everything in the database
    deletedf = _get_left_diff_df(databasedf, new_datasetdf, checkby)
    if not deletedf.empty:
        logger.info("Deleting Rows")
        delete_rowid_version = pd.DataFrame(
            [[rowid.split("_")[0], rowid.split("_")[1]] for rowid in deletedf.index]
        )
        delete_rowid_version.reset_index(drop=True, inplace=True)
    else:
        delete_rowid_version = pd.DataFrame()
        logger.info("No deleted rows")

    # del deletedf[checkby]
    return delete_rowid_version


def _create_update_rowsdf(updating_databasedf, updatesetdf, rowids, differentrows):
    """
    Create the update dataset dataframe

    Args:
        updating_databasedf: Update database dataframe
        updatesetdf:  Update dataset dataframe
        rowids: rowids of the database (Synapse ROW_ID, ROW_VERSION)
        differentrows: vector of booleans for rows that need to be updated
                       True for update, False for not

    Returns:
        dataframe: Update dataframe
    """
    if sum(differentrows) > 0:
        updating_databasedf.loc[differentrows] = updatesetdf.loc[differentrows]
        toupdatedf = updating_databasedf.loc[differentrows]
        logger.info("Updating rows")
        rowid_version = pd.DataFrame(
            [
                [rowid.split("_")[0], rowid.split("_")[1]]
                for rowid, row in zip(rowids, differentrows)
                if row
            ]
        )
        toupdatedf["ROW_ID"] = rowid_version[0].values
        toupdatedf["ROW_VERSION"] = rowid_version[1].values
        toupdatedf.reset_index(drop=True, inplace=True)
    else:
        toupdatedf = pd.DataFrame()
        logger.info("No updated rows")
    return toupdatedf


def _update_rows(new_datasetdf, databasedf, checkby):
    """
    Compares the dataset from the database and determines which rows to
    update from the dataset

    Args:
        new_datasetdf: Input data dataframe
        databasedf: Existing data dataframe
        checkby: Column of values to compare

    Return:
        Dataframe: Dataframe of rows to update
    """
    # initial_database = databasedf.copy()
    databasedf.fillna("", inplace=True)
    new_datasetdf.fillna("", inplace=True)
    updatesetdf = _get_left_union_df(new_datasetdf, databasedf, checkby)
    updating_databasedf = _get_left_union_df(databasedf, new_datasetdf, checkby)

    # If you input the exact same dataframe theres nothing to update
    # must save row version and ids for later
    rowids = updating_databasedf.index.values
    # Set index values to be 'checkby' values
    updatesetdf.index = updatesetdf[checkby]
    updating_databasedf.index = updating_databasedf[checkby]
    del updatesetdf[checkby]
    del updating_databasedf[checkby]

    # Remove duplicated index values
    updatesetdf = updatesetdf[~updatesetdf.index.duplicated()]
    # Reorder dataset index
    updatesetdf = updatesetdf.loc[updating_databasedf.index]
    # Index comparison
    differences = updatesetdf != updating_databasedf
    differentrows = differences.apply(sum, axis=1) > 0

    toupdatedf = _create_update_rowsdf(
        updating_databasedf, updatesetdf, rowids, differentrows
    )

    return toupdatedf


def updateData(
    syn: Synapse,
    databaseSynId: str,
    newData: pd.DataFrame,
    filterBy: str,
    filterByColumn: str = "CENTER",
    col: List[str] = None,
    toDelete: bool = False,
) -> None:
    """Update Synapse table given a new dataframe

    Args:
        syn (Synapse): Synapse connection
        databaseSynId (str): Synapse Id of Synapse Table
        newData (pd.DataFrame): New data in a dataframe
        filterBy (str): Value to filter new data by
        filterByColumn (str, optional): Column to filter values by. Defaults to "CENTER".
        col (List[str], optional): List of columns to ingest. Defaults to None.
        toDelete (bool, optional): Delete rows given the primary key. Defaults to False.
    """
    databaseEnt = syn.get(databaseSynId)
    database = syn.tableQuery(
        f"SELECT * FROM {databaseSynId} where {filterByColumn} ='{filterBy}'"
    )
    database = database.asDataFrame()
    db_cols = set(database.columns)
    if col is not None:
        new_data_cols = set(col)
        # Make sure columns from file exists in database columns
        use_cols = db_cols.intersection(new_data_cols)
        # No need to fail, because there is bound to be at least one
        # column that will exist in the database
        database = database[list(use_cols)]
    else:
        newData = newData[database.columns]
    updateDatabase(
        syn, database, newData, databaseSynId, databaseEnt.primaryKey, toDelete
    )


def updateDatabase(
    syn, database, new_dataset, database_synid, primary_key_cols, to_delete=False
):
    """
    Updates synapse tables by a row identifier with another
    dataset that has the same number and order of columns

    Args:
        syn: Synapse object
        database: The synapse table (pandas dataframe)
        new_dataset: New dataset (pandas dataframe)
        databaseSynId: Synapse Id of the database table
        uniqueKeyCols: Column(s) that make up the unique key
        toDelete: Delete rows, Defaults to False

    Returns:
        Nothing
    """
    primary_key = "UNIQUE_KEY"
    database = database.fillna("")
    orig_database_cols = database.columns
    col_order = ["ROW_ID", "ROW_VERSION"]
    col_order.extend(orig_database_cols.tolist())
    new_dataset = new_dataset.fillna("")
    # Columns must be in the same order
    new_dataset = new_dataset[orig_database_cols]
    database[primary_key_cols] = database[primary_key_cols].applymap(str)
    database[primary_key] = database[primary_key_cols].apply(
        lambda x: " ".join(x), axis=1
    )

    new_dataset[primary_key_cols] = new_dataset[primary_key_cols].applymap(str)
    new_dataset[primary_key] = new_dataset[primary_key_cols].apply(
        lambda x: " ".join(x), axis=1
    )

    allupdates = pd.DataFrame(columns=col_order)
    to_append_rows = _append_rows(new_dataset, database, primary_key)
    to_update_rows = _update_rows(new_dataset, database, primary_key)
    if to_delete:
        to_delete_rows = _delete_rows(new_dataset, database, primary_key)
    else:
        to_delete_rows = pd.DataFrame()
    allupdates = pd.concat([allupdates, to_append_rows, to_update_rows], sort=False)
    storedatabase = False
    update_all_file = tempfile.NamedTemporaryFile(dir=SCRIPT_DIR, delete=False)

    with open(update_all_file.name, "w") as updatefile:
        # Must write out the headers in case there are no appends or updates
        updatefile.write(",".join(col_order) + "\n")
        if not allupdates.empty:
            """
            This is done because of pandas typing.
            An integer column with one NA/blank value
            will be cast as a double.
            """
            updatefile.write(
                allupdates[col_order]
                .to_csv(index=False, header=None)
                .replace(".0,", ",")
                .replace(".0\n", "\n")
            )
            storedatabase = True
        if not to_delete_rows.empty:
            updatefile.write(
                to_delete_rows.to_csv(index=False, header=None)
                .replace(".0,", ",")
                .replace(".0\n", "\n")
            )
            storedatabase = True
    if storedatabase:
        syn.store(synapseclient.Table(database_synid, update_all_file.name))
    # Delete the update file
    os.unlink(update_all_file.name)


def checkInt(element):
    """
    Check if an item can become an integer

    Args:
        element: Any variable and type

    Returns:
        boolean: True/False
    """
    try:
        element = float(element)
        return element.is_integer()
    except (ValueError, TypeError):
        return False


def check_col_and_values(
    df, col, possible_values, filename, na_allowed=False, required=False, sep=None
):
    """
    This function checks if the column exists then checks if the values in the
    column have the correct values

    Args:
        df: Input dataframe
        col: Expected column name
        possible_values: list of possible values
        filename: Name of file
        required: If the column is required.  Default is False

    Returns:
        tuple: warning, error
    """
    warning = ""
    error = ""
    have_column = checkColExist(df, col)
    if not have_column:
        if required:
            error = "{filename}: Must have {col} column.\n".format(
                filename=filename, col=col
            )
        else:
            warning = (
                "{filename}: Doesn't have {col} column. "
                "This column will be added\n".format(filename=filename, col=col)
            )
    else:
        if na_allowed:
            check_values = df[col].dropna()
        else:
            check_values = df[col]
        if sep:
            final = []
            for value in check_values:
                final.extend(value.split(sep))
            check_values = pd.Series(final)
        if not check_values.isin(possible_values).all():
            error = "{filename}: Please double check your {col} column.  This column must only be these values: {possible_vals}\n".format(
                filename=filename,
                col=col,
                possible_vals=", ".join(
                    [
                        # This is done because of pandas typing.
                        # An integer column with one NA/blank value
                        # will be cast as a double.
                        str(value).replace(".0", "")
                        for value in possible_values
                    ]
                ),
            )
    return (warning, error)


def extract_oncotree_code_mappings_from_oncotree_json(
    oncotree_json, primary, secondary
):
    oncotree_code_to_info = {}
    data = oncotree_json["children"]
    for node in data:
        # if not node['code']:
        #     sys.stderr.write('Encountered oncotree node without '
        #                       'oncotree code : ' + node + '\n')
        #     continue
        if data[node]["level"] == 1:
            primary = node
            secondary = ""
        elif data[node]["level"] == 2:
            secondary = node
        cancer_type = data[node]["mainType"]
        cancer_type_detailed = data[node]["name"]
        if not cancer_type_detailed:
            cancer_type_detailed = ""
        oncotree_code_to_info[node.upper()] = {
            "CANCER_TYPE": cancer_type,
            "CANCER_TYPE_DETAILED": cancer_type_detailed,
            "ONCOTREE_PRIMARY_NODE": primary,
            "ONCOTREE_SECONDARY_NODE": secondary,
        }

        if len(data[node]["children"]) > 0:
            recurseDict = extract_oncotree_code_mappings_from_oncotree_json(
                data[node], primary, secondary
            )
            oncotree_code_to_info.update(recurseDict)
    return oncotree_code_to_info


def get_oncotree_code_mappings(oncotree_tumortype_api_endpoint_url):
    """
    CREATE ONCOTREE DICTIONARY MAPPING TO PRIMARY, SECONDARY,
    CANCER TYPE, AND CANCER DESCRIPTION
    """
    # oncotree_raw_response = urlopen(oncotree_tumortype_api_endpoint_url).text
    # with requests.get(oncotree_tumortype_api_endpoint_url) as oncotreeUrl:
    oncotreeUrl = retry_get_url(oncotree_tumortype_api_endpoint_url)
    oncotree_raw_response = oncotreeUrl.text
    oncotree_response = json.loads(oncotree_raw_response)
    oncotree_response = oncotree_response["TISSUE"]
    return extract_oncotree_code_mappings_from_oncotree_json(oncotree_response, "", "")


# Get mapping code #Add USE DESCRIPTION sampletypedetailed -> public
def getCODE(mapping, key, useDescription=False):
    if useDescription:
        value = mapping["DESCRIPTION"][mapping["CODE"] == key].values
    else:
        value = mapping["CBIO_LABEL"][mapping["CODE"] == key].values
    if len(value) > 0:
        return value[0]
    else:
        return ""


def getPrimary(code, oncotreeDict, primary):
    if code != "":
        for level in oncotreeDict:
            if sum(oncotreeDict[level] == code) > 0:
                toAdd = primary[oncotreeDict[level] == code].values[0]
                break
            else:
                toAdd = code
    else:
        toAdd = "NOT_ANNOTATED"
    return toAdd


# def create_key():
#   from Crypto.PublicKey import RSA
#   from Crypto import Random
#   from Crypto.Cipher import PKCS1_OAEP
#
#   random_generator = Random.new().read
#   generate public and private keys
#   key = RSA.generate(1024, random_generator)
#   key_cryptor = PKCS1_OAEP.new(key)
#   encrypted = key_cryptor.encrypt(geniePassword)
#   #message to encrypt is in the above line 'encrypt this message'
#   decrypted = key_cryptor.decrypt(encrypted)
#   with open("genie.pem","wb") as geniepem_f:
#       geniepem_f.write(key.exportKey(format='PEM'))


def read_key(pemfile_path):
    """
    Obtain key from pemfile

    Args:
        pemfile_path:  Path to pemfile

    Returns:
        RSA key
    """
    with open(pemfile_path, "r") as pemfile_f:
        key = RSA.importKey(pemfile_f.read())
    return key


def decrypt_message(message, key):
    """
    Decrypt message with a pem key from
    func read_key

    Args:
        message: Encrypted message
        key: read_key returned key

    Returns:
        Decrypted message
    """
    # Use module Crypto.Cipher.PKCS1_OAEP instead
    decryptor = PKCS1_OAEP.new(key)
    decrypted = decryptor.decrypt(ast.literal_eval(str(message)))
    return decrypted.decode("utf-8")


def get_password(pemfile_path):
    """
    Get password using pemfile

    Args:
        pemfile_path: Path to pem file

    Return:
        Password
    """
    if not os.path.exists(pemfile_path):
        raise ValueError(
            "Path to pemFile must be specified if there " "is no cached credentials"
        )
    key = read_key(pemfile_path)
    genie_pass = decrypt_message(os.environ["GENIE_PASS"], key)
    return genie_pass


def synLogin(pemfile_path, debug=False):
    """
    Use pem file to log into synapse if credentials aren't cached

    Args:
        pemfile_path: Path to pem file
        debug: Synapse debug feature.  Defaults to False

    Returns:
        Synapse object logged in
    """
    try:
        syn = synapseclient.Synapse(debug=debug)
        # Get auth token via scheduled job secrets
        if os.getenv("SCHEDULED_JOB_SECRETS") is not None:
            secrets = json.loads(os.getenv("SCHEDULED_JOB_SECRETS"))
            auth_token = secrets["SYNAPSE_AUTH_TOKEN"]
        else:
            auth_token = None
        syn.login(authToken=auth_token)
    except Exception:
        # TODO: deprecate this feature soon
        genie_pass = get_password(pemfile_path)
        syn = synapseclient.Synapse(debug=debug)
        syn.login(os.environ["GENIE_USER"], genie_pass)
    return syn


def get_gdc_data_dictionary(filetype):
    """
    Use the GDC API to get the values allowed for columns of
    different filetypes (ie. disease_type in the case file)

    Args:
        filetype: GDC file type (ie. case, read_group)

    Return:
        json:  Dictionary of allowed columns for the filetype and
               allowed values for those columns
    """
    gdc_dict = retry_get_url(
        "https://api.gdc.cancer.gov/v0/submission/_dictionary/{filetype}".format(
            filetype=filetype
        )
    )
    gdc_response = json.loads(gdc_dict.text)
    return gdc_response


def _create_schema(syn, table_name, parentid, columns=None, annotations=None):
    """Creates Table Schema

    Args:
        syn: Synapse object
        table_name: Name of table
        parentid: Project synapse id
        columns: Columns of Table
        annotations: Dictionary of annotations to add

    Returns:
        Schema
    """
    schema = synapseclient.Schema(
        name=table_name, columns=columns, parent=parentid, annotations=annotations
    )
    new_schema = syn.store(schema)
    return new_schema


def _update_database_mapping(
    syn, database_synid_mappingdf, database_mapping_synid, fileformat, new_tableid
):
    """Updates database to synapse id mapping table

    Args:
        syn: Synapse object
        database_synid_mappingdf: Database to synapse id mapping dataframe
        database_mapping_synid: Database to synapse id table id
        fileformat: File format updated
        new_tableid: New file format table id

    Returns:
        Updated Table object
    """
    fileformat_ind = database_synid_mappingdf["Database"] == fileformat
    # Store in the new database synid
    database_synid_mappingdf["Id"][fileformat_ind] = new_tableid
    # Only update the one row
    to_update_row = database_synid_mappingdf[fileformat_ind]

    syn.store(synapseclient.Table(database_mapping_synid, to_update_row))
    return database_synid_mappingdf


# TODO: deprecate once move function is out of the cli into the
# client master branch
def _move_entity(syn, ent, parentid, name=None):
    """Moves an entity (works like linux mv)

    Args:
        syn: Synapse object
        ent: Synapse Entity
        parentid: Synapse Project id
        name: New Entity name if a new name is desired

    Returns:
        Moved Entity
    """
    ent.parentId = parentid
    if name is not None:
        ent.name = name
    moved_ent = syn.store(ent)
    return moved_ent


def get_dbmapping(syn: Synapse, projectid: str) -> dict:
    """Gets database mapping information

    Args:
        syn: Synapse connection
        projectid: Project id where new data lives

    Returns:
        {'synid': database mapping syn id,
         'df': database mapping pd.DataFrame}

    """
    project_ent = syn.get(projectid)
    dbmapping_synid = project_ent.annotations.get("dbMapping", "")[0]
    database_mapping = syn.tableQuery(f"select * from {dbmapping_synid}")
    database_mappingdf = database_mapping.asDataFrame()
    return {"synid": dbmapping_synid, "df": database_mappingdf}


def create_new_fileformat_table(
    syn: Synapse,
    file_format: str,
    newdb_name: str,
    projectid: str,
    archive_projectid: str,
) -> dict:
    """Creates new database table based on old database table and archives
    old database table

    Args:
        syn: Synapse object
        file_format: File format to update
        newdb_name: Name of new database table
        projectid: Project id where new database should live
        archive_projectid: Project id where old database should be moved

    Returns:
        {"newdb_ent": New database synapseclient.Table,
         "newdb_mappingdf": new databse pd.DataFrame,
         "moved_ent": old database synpaseclient.Table}
    """
    db_info = get_dbmapping(syn, projectid)
    database_mappingdf = db_info["df"]
    dbmapping_synid = db_info["synid"]

    olddb_synid = getDatabaseSynId(
        syn, file_format, databaseToSynIdMappingDf=database_mappingdf
    )
    olddb_ent = syn.get(olddb_synid)
    olddb_columns = list(syn.getTableColumns(olddb_synid))

    newdb_ent = _create_schema(
        syn,
        table_name=newdb_name,
        columns=olddb_columns,
        parentid=projectid,
        annotations=olddb_ent.annotations,
    )

    newdb_mappingdf = _update_database_mapping(
        syn, database_mappingdf, dbmapping_synid, file_format, newdb_ent.id
    )
    # Automatically rename the archived entity with ARCHIVED
    # This will attempt to resolve any issues if the table already exists at
    # location
    new_table_name = f"ARCHIVED {time.time()}-{olddb_ent.name}"
    moved_ent = _move_entity(syn, olddb_ent, archive_projectid, name=new_table_name)
    return {
        "newdb_ent": newdb_ent,
        "newdb_mappingdf": newdb_mappingdf,
        "moved_ent": moved_ent,
    }


def get_database_mapping_config(syn: Synapse, synid: str) -> dict:
    """Gets Synapse database to Table mapping in dict

    Args:
        syn: Synapse connection
        synid: Synapse id of database mapping table

    Returns:
        dict: {'databasename': 'synid'}
    """
    config = syn.tableQuery(f"SELECT * FROM {synid}")
    configdf = config.asDataFrame()
    configdf.index = configdf["Database"]
    config_dict = configdf.to_dict()
    return config_dict["Id"]


def get_genie_config(
    syn: Synapse,
    project_id: str,
) -> dict:
    """Get configurations needed for the GENIE codebase

    Args:
        syn (Synapse): Synapse connection
        project_id (str): Synapse project id

    Returns:
        dict: _description_
    """
    # Get the Synapse Project where data is stored
    # Should have annotations to find the table lookup
    project = syn.get(project_id)

    # Get project GENIE configurations
    database_to_synid_mapping_synid = project.annotations.get("dbMapping", "")
    genie_config = get_database_mapping_config(
        syn=syn, synid=database_to_synid_mapping_synid[0]
    )
    # Fill in GENIE center configurations
    center_mapping_id = genie_config["centerMapping"]
    center_mapping_df = get_syntabledf(
        syn=syn, query_string=f"SELECT * FROM {center_mapping_id} where release is true"
    )
    center_mapping_df.index = center_mapping_df.center
    # Add center configurations including input/staging synapse ids
    genie_config["center_config"] = center_mapping_df.to_dict("index")

    genie_config["ethnicity_mapping"] = "syn7434242"
    genie_config["race_mapping"] = "syn7434236"
    genie_config["sex_mapping"] = "syn7434222"
    genie_config["sampletype_mapping"] = "syn7434273"

    return genie_config


def _get_oncotreelink(syn, genie_config, oncotree_link=None):
    """
    Gets oncotree link unless a link is specified by the user

    Args:
        syn: Synapse object
        databasetosynid_mappingdf: database to synid mapping
        oncotree_link: link to oncotree. Default is None

    Returns:
        oncotree link

    """
    if oncotree_link is None:
        onco_link_ent = syn.get(genie_config["oncotreeLink"])
        oncotree_link = onco_link_ent.externalURL
    return oncotree_link
