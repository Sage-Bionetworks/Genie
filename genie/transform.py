"""This module contains all the transformation functions used throughout the GENIE
package"""
import warnings

import pandas as pd
from pandas.api.types import is_float_dtype


def _col_name_to_titlecase(string: str) -> str:
    """Convert strings to titlecase. Supports strings separated by _.

    Args:
        string (str): A string

    Returns:
        str: A string converted to title case

    """
    # This is a mapping of strings that are abbreviations
    abbrev_map = {"Dna_": "DNA_", "Rna_": "RNA_", "Sv_": "SV_", "Ncbi_": "NCBI_"}
    # The reason I split the string by _ and utilize the capitalize function
    # instead of using something like .title(), so outlined here:
    # https://stackoverflow.com/questions/1549641/how-can-i-capitalize-the-first-letter-of-each-word-in-a-string
    converted_str = "_".join([each.capitalize() for each in string.split("_")])
    for titlecase, abbrev in abbrev_map.items():
        converted_str = converted_str.replace(titlecase, abbrev)
    return converted_str


def _convert_col_with_nas_to_str(df: pd.DataFrame, col: str) -> list:
    """This converts a column into str while preserving NAs"""
    new_vals = [str(val) if pd.notna(val) else val for val in df[col]]
    return new_vals


def _convert_float_col_with_nas_to_int(df: pd.DataFrame, col: str) -> list:
    """This converts int column that was turned into a float col because
    pandas does that with int values that have NAs back into an int col
    with NAs intact"""
    if is_float_dtype(df[col]) and df[col].isnull().values.any():
        new_vals = df[col].astype(pd.Int64Dtype()).tolist()
        return new_vals
    else:
        return df[col].tolist()


def _convert_df_with_mixed_dtypes(read_csv_params: dict) -> pd.DataFrame:
    """This checks if a dataframe read in normally comes out with mixed data types (which happens
    when low_memory = True because read_csv parses in chunks and guesses dtypes by chunk) and
    converts a dataframe with mixed datatypes to one datatype.

    Args:
        read_csv_params (dict): of input params and values to pandas's read_csv function.
            needs to include filepath to dataset to be read in

    Returns:
        pd.DataFrame : The dataset read in
    """
    warnings.simplefilter("error", pd.errors.DtypeWarning)
    try:
        df = pd.read_csv(**read_csv_params, low_memory=True)
    except pd.errors.DtypeWarning:
        # setting engine to c as that is the only engine that works with low_memory=False
        df = pd.read_csv(**read_csv_params, low_memory=False, engine="c")
    warnings.resetwarnings()
    return df
