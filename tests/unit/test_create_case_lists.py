# import pytest
from genie import create_case_lists
import os

study_id = "test"

clinical_file_map = {'': ['FOOBAR', 'NEED']}
clinical_file_map['Testing2'] = ['test1']
clinical_file_map['Test1 Now, Please/foo'] = ['wow']


case_list_files = create_case_lists.write_case_list_files(
    clinical_file_map, "./", study_id)
sequenced_case_list_files = create_case_lists.write_case_list_sequenced(
    ['test1', 'test2'], "./", study_id)
case_list_cna_path = create_case_lists.write_case_list_cna(
    ['test1', 'test2'], "./", study_id)
case_list_cnaseq_path = create_case_lists.write_case_list_cnaseq(
    ['test1', 'test2'], "./", study_id)


def test_filenames_write_case_list_files():
    first = os.path.basename(case_list_files[2])
    assert first == "cases_Test1_Now_Please_foo.txt"
    second = os.path.basename(case_list_files[0])
    assert second == "cases_no_oncotree_code.txt"
    third = os.path.basename(case_list_files[1])
    assert third == "cases_Testing2.txt"


def test_textfix_write_case_list_files():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_Test1_Now_Please_foo\n'
        'case_list_name: Tumor Type: Test1 Now, Please/foo\n'
        'case_list_description: All tumors with cancer type '
        'Test1 Now, Please/foo\n'
        'case_list_ids: wow')
    with open(case_list_files[2], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_files[2])


def test_nocode_write_case_list_files():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_no_oncotree_code\n'
        'case_list_name: Tumor Type: NA\n'
        'case_list_description: All tumors with cancer type NA\n'
        'case_list_ids: FOOBAR\tNEED')
    with open(case_list_files[0], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_files[0])


def test_nochange_write_case_list_files():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_Testing2\n'
        'case_list_name: Tumor Type: Testing2\n'
        'case_list_description: All tumors with cancer type Testing2\n'
        'case_list_ids: test1')
    with open(case_list_files[1], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_files[1])


def test_filenames_write_case_list_sequenced():
    first = os.path.basename(sequenced_case_list_files[0])
    assert first == "cases_sequenced.txt"
    second = os.path.basename(sequenced_case_list_files[1])
    assert second == "cases_all.txt"


def test_sequencetext_write_case_list_sequenced():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_sequenced\n'
        'case_list_name: Sequenced Tumors\n'
        'case_list_description: All sequenced samples\n'
        'case_list_ids: test1\ttest2')
    with open(sequenced_case_list_files[0], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(sequenced_case_list_files[0])


def test_all_write_case_list_sequenced():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_all\n'
        'case_list_name: All samples\n'
        'case_list_description: All samples\n'
        'case_list_ids: test1\ttest2')
    with open(sequenced_case_list_files[1], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(sequenced_case_list_files[1])


def test_filename_write_case_list_cna():
    assert os.path.basename(case_list_cna_path) == "cases_cna.txt"


def test_filename_write_case_list_cnaseq():
    assert os.path.basename(case_list_cnaseq_path) == "cases_cnaseq.txt"


def test_cnatext_write_case_list_cna():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_cna\n'
        'case_list_name: Samples with CNA\n'
        'case_list_description: Samples with CNA\n'
        'case_list_ids: test1\ttest2')
    with open(case_list_cna_path, 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_cna_path)


def test_cnaseq_write_case_list_cnaseq():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_cnaseq\n'
        'case_list_name: Samples with CNA and mutation\n'
        'case_list_description: Samples with CNA and mutation\n'
        'case_list_ids: test1\ttest2')
    with open(case_list_cnaseq_path, 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_cnaseq_path)
