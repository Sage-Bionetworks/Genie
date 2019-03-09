import synapseclient
import pandas as pd
import mock
from genie.assay import Assayinfo
import pytest

syn = synapseclient.Synapse()
assay_info = Assayinfo(syn, "SAGE")

def test_filetype():
    assert assay_info._fileType == "assayinfo"

def test_incorrect_validatefilename():
    with pytest.raises(AssertionError):
        assay_info.validateFilename(['foo'])

def test_correct_validatefilename():
    assert assay_info.validateFilename(["assay_information.yaml"]) == "assayinfo"

def test_valid__validate():
    assay_info_dict = {'SEQ_ASSAY_ID':['SAGE-1','SAGE-3'],
    'is_paired_end':[True,False],
    'library_strategy':['ChIP-Seq','ChIP-Seq'],
    'library_selection':['PCR','PCR'],
    'platform':['Illumina','Illumina'],
    'instrument_model':['Illumina HiSeq 4000','Illumina HiSeq 4000'],
    'variant_classifications':['Frame_Shift_Ins','Frame_Shift_Ins'],
    'target_capture_kit':['foo','doo'],
    'read_length':[22,333],
    'number_of_genes':[5,20],
    'gene_padding':[10, None]}
    assay_info_df = pd.DataFrame(assay_info_dict)
    error, warning = assay_info._validate(assay_info_df)
    assert error == ''
    assert warning == '' 

def test_missingcols__validate():
    assay_info_df = pd.DataFrame()
    error, warning = assay_info._validate(assay_info_df)
    expected_errors = ('Assay_information.yaml: Must have SEQ_ASSAY_ID column.\n'
                       'Assay_information.yaml: Must have is_paired_end column.\n'
                       'Assay_information.yaml: Must have library_selection column.\n'
                       'Assay_information.yaml: Must have library_strategy column.\n'
                       'Assay_information.yaml: Must have platform column.\n'
                       'Assay_information.yaml: Must have instrument_model column.\n'
                       'Assay_information.yaml: Must have target_capture_kit column.\n'
                       'Assay_information.yaml: Must have read_length column.\n'
                       'Assay_information.yaml: Must have number_of_genes column.\n')
    assert error == expected_errors
    expected_warnings = ("Assay_information.yaml: Doesn't have variant_classifications column. This column will be added\n"
                         "Assay_information.yaml: gene_padding is by default 10 if not specified.\n")
    assert warning == expected_warnings

def test_fillcols__process():
    '''
    Standardization of SEQ_ASSAY_ID
    Add in CENTER, gene_padding, and variant_classifications if missing
    '''

    assay_info_dict = {'SEQ_ASSAY_ID':['SAGE-Foo_1']}
    assay_info_df = pd.DataFrame(assay_info_dict)
    processed_assay_df = assay_info._process(assay_info_df)
    expected_assay_df = pd.DataFrame({'SEQ_ASSAY_ID':['SAGE-FOO-1'],
                                      'gene_padding':[10],
                                      'variant_classifications':[float('nan')],
                                      'CENTER':['SAGE']})
    assert expected_assay_df.equals(processed_assay_df[expected_assay_df.columns])

def test_default10__process():
    '''
    gene_padding default 10
    '''
    
    assay_info_dict = {'SEQ_ASSAY_ID':['SAGE-1','SAGE-2'],
                       'gene_padding':[20,float('nan')],
                       'variant_classifications':['test','test']}
    assay_info_df = pd.DataFrame(assay_info_dict)
    processed_assay_df = assay_info._process(assay_info_df)
    expected_assay_df = pd.DataFrame({'SEQ_ASSAY_ID':['SAGE-1','SAGE-2'],
                                      'gene_padding':[20, 10],
                                      'variant_classifications':['test','test'],
                                      'CENTER':['SAGE','SAGE']})
    assert expected_assay_df.equals(processed_assay_df[expected_assay_df.columns])