# import pytest
from genie import create_case_lists
import os

study_id = "test"

clinical_file_map = {
    '': ['FOOBAR', 'NEED'],
    'UNKNOWN': ['test1'],
    "Test1 Now, Please/foo": ["wow"]}

case_list_files = create_case_lists.write_case_list_files(
    clinical_file_map, "./", study_id)
#sequenced_case_list_files = create_case_lists.write_case_list_sequenced(
#    ['test1', 'test2'], "./", study_id)


def test_filenames_write_case_list_files():
    first = os.path.basename(case_list_files[0])
    assert first == "cases_Test1_Now_Please_foo.txt"
    second = os.path.basename(case_list_files[1])
    assert second == "cases_no_oncotree_code.txt"
    third = os.path.basename(case_list_files[2])
    assert third == "cases_UNKNOWN.txt"


def test_textfix_write_case_list_files():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_Test1_Now_Please_foo\n'
        'case_list_name: Tumor Type: Test1 Now, Please/foo\n'
        'case_list_description: All tumors with cancer type '
        'Test1 Now, Please/foo\n'
        'case_list_ids: wow')
    with open(case_list_files[0], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_files[0])


def test_nocode_write_case_list_files():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_no_oncotree_code\n'
        'case_list_name: Tumor Type: NA\n'
        'case_list_description: All tumors with cancer type NA\n'
        'case_list_ids: FOOBAR\tNEED')
    with open(case_list_files[1], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_files[1])


def test_nochange_write_case_list_files():
    expected_text = (
        'cancer_study_identifier: test\n'
        'stable_id: test_UNKNOWN\n'
        'case_list_name: Tumor Type: UNKNOWN\n'
        'case_list_description: All tumors with cancer type UNKNOWN\n'
        'case_list_ids: test1')
    with open(case_list_files[2], 'r') as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_files[2])
