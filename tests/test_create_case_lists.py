import os

import pytest

from genie import create_case_lists

study_id = "test"

clinical_file_map = {"": ["FOOBAR", "NEED"]}
clinical_file_map["Testing2"] = ["test1"]
clinical_file_map["Test1 Now, Please/foo"] = ["wow"]

sequenced_case_list_files = create_case_lists.write_case_list_sequenced(
    ["test1", "test2"], "./", study_id
)
case_list_cna_path = create_case_lists.write_case_list_cna(
    ["test1", "test2"], "./", study_id
)
case_list_cnaseq_path = create_case_lists.write_case_list_cnaseq(
    ["test1", "test2"], "./", study_id
)


def test_filenames_write_case_list_files():
    case_list_files = create_case_lists.write_case_list_files(
        clinical_file_map, "./", study_id
    )
    required_files = [
        "cases_Test1_Now_Please_foo.txt",
        "cases_no_oncotree_code.txt",
        "cases_Testing2.txt",
    ]
    basenames = [os.path.basename(case_file) for case_file in case_list_files]

    assert all([req_file in basenames for req_file in required_files])


expected_change_text = (
    "cancer_study_identifier: test\n"
    "stable_id: test_Test1_Now_Please_foo\n"
    "case_list_name: Tumor Type: Test1 Now, Please/foo\n"
    "case_list_description: All tumors with cancer type "
    "Test1 Now, Please/foo\n"
    "case_list_ids: wow"
)

expected_same_text = (
    "cancer_study_identifier: test\n"
    "stable_id: test_Testing2\n"
    "case_list_name: Tumor Type: Testing2\n"
    "case_list_description: All tumors with cancer type Testing2\n"
    "case_list_ids: test1"
)

expected_nocode_text = (
    "cancer_study_identifier: test\n"
    "stable_id: test_no_oncotree_code\n"
    "case_list_name: Tumor Type: NA\n"
    "case_list_description: All tumors with cancer type NA\n"
    "case_list_ids: FOOBAR\tNEED"
)


@pytest.fixture(
    params=[
        # tuple with (input, expectedOutput)
        ("", clinical_file_map[""], expected_nocode_text),
        ("Testing2", clinical_file_map["Testing2"], expected_same_text),
        (
            "Test1 Now, Please/foo",
            clinical_file_map["Test1 Now, Please/foo"],
            expected_change_text,
        ),
    ]
)
def oncotree_write_params(request):
    return request.param


def test__write_single_oncotree_case_list(oncotree_write_params):
    (cancer_type, ids, expected_text) = oncotree_write_params
    caselist_path = create_case_lists._write_single_oncotree_case_list(
        cancer_type, ids, study_id, "./"
    )
    with open(caselist_path, "r") as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(caselist_path)


def test_filenames_write_case_list_sequenced():
    first = os.path.basename(sequenced_case_list_files[0])
    assert first == "cases_sequenced.txt"
    second = os.path.basename(sequenced_case_list_files[1])
    assert second == "cases_all.txt"


def test_sequencetext_write_case_list_sequenced():
    expected_text = (
        "cancer_study_identifier: test\n"
        "stable_id: test_sequenced\n"
        "case_list_name: Sequenced Tumors\n"
        "case_list_description: All sequenced samples\n"
        "case_list_ids: test1\ttest2"
    )
    with open(sequenced_case_list_files[0], "r") as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(sequenced_case_list_files[0])


def test_all_write_case_list_sequenced():
    expected_text = (
        "cancer_study_identifier: test\n"
        "stable_id: test_all\n"
        "case_list_name: All samples\n"
        "case_list_description: All samples\n"
        "case_list_ids: test1\ttest2"
    )
    with open(sequenced_case_list_files[1], "r") as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(sequenced_case_list_files[1])


def test_filename_write_case_list_cna():
    assert os.path.basename(case_list_cna_path) == "cases_cna.txt"


def test_filename_write_case_list_cnaseq():
    assert os.path.basename(case_list_cnaseq_path) == "cases_cnaseq.txt"


def test_cnatext_write_case_list_cna():
    expected_text = (
        "cancer_study_identifier: test\n"
        "stable_id: test_cna\n"
        "case_list_name: Samples with CNA\n"
        "case_list_description: Samples with CNA\n"
        "case_list_ids: test1\ttest2"
    )
    with open(case_list_cna_path, "r") as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_cna_path)


def test_cnaseq_write_case_list_cnaseq():
    expected_text = (
        "cancer_study_identifier: test\n"
        "stable_id: test_cnaseq\n"
        "case_list_name: Samples with CNA and mutation\n"
        "case_list_description: Samples with CNA and mutation\n"
        "case_list_ids: test1\ttest2"
    )
    with open(case_list_cnaseq_path, "r") as case_list:
        caselist_text = case_list.read()
    assert caselist_text == expected_text
    os.remove(case_list_cnaseq_path)
