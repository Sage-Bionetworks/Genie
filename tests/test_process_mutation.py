"""Test process mutation functions"""
import pandas as pd

from genie import process_mutation


def test_format_maf():
    maf_dict = {}
    maf_dict['Center'] = ["foo", "dsdf", "sdf"]
    maf_dict['Tumor_Sample_Barcode'] = ["GENIE-SAGE-1-3", "1-2", "3-2"]
    maf_dict['Sequence_Source'] = ["3", "e", "sd"]
    maf_dict['Sequencer'] = ["dsf", "sdf", "d"]
    maf_dict['Validation_Status'] = ["Unknown", "unknown", "f"]
    mafdf = pd.DataFrame(maf_dict)

    formatted_mafdf = process_mutation.format_maf(mafdf, center="SAGE")

    expected_maf_dict = {}
    expected_maf_dict['Center'] = ["SAGE", "SAGE", "SAGE"]
    expected_maf_dict['Tumor_Sample_Barcode'] = [
        "GENIE-SAGE-1-3", "GENIE-SAGE-1-2", "GENIE-SAGE-3-2"
    ]
    expected_maf_dict['Sequence_Source'] = [float('nan'), float('nan'),
                                            float('nan')]
    expected_maf_dict['Sequencer'] = [float('nan'), float('nan'),
                                      float('nan')]
    expected_maf_dict['Validation_Status'] = ['', '', "f"]
    expected_mafdf = pd.DataFrame(expected_maf_dict)
    assert expected_mafdf.equals(formatted_mafdf[expected_mafdf.columns])
