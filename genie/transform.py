"""This module contains all the transformation functions used throughout the GENIE
package"""


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
