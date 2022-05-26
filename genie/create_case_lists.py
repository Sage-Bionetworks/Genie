"""
Creates case lists per cancer type
"""
from collections import defaultdict
import csv
import os

CASE_LIST_TEXT_TEMPLATE = (
    "cancer_study_identifier: {study_id}\n"
    "stable_id: {stable_id}\n"
    "case_list_name: {case_list_name}\n"
    "case_list_description: {case_list_description}\n"
    "case_list_ids: {case_list_ids}"
)


def create_case_lists_map(clinical_file_name):
    """
    Creates the case list dictionary

    Args:
        clinical_file_name: clinical file path

    Returns:
        dict: key = cancer_type
              value = list of sample ids
        dict: key = seq_assay_id
              value = list of sample ids
        list: Clinical samples
    """
    with open(clinical_file_name, "rU") as clinical_file:
        seq_assay_map = defaultdict(list)
        clinical_file_map = defaultdict(list)
        clin_samples = []
        reader = csv.DictReader(clinical_file, dialect="excel-tab")
        for row in reader:
            clinical_file_map[row["CANCER_TYPE"]].append(row["SAMPLE_ID"])
            seq_assay_map[row["SEQ_ASSAY_ID"]].append(row["SAMPLE_ID"])
            clin_samples.append(row["SAMPLE_ID"])
    return clinical_file_map, seq_assay_map, clin_samples


def _write_single_oncotree_case_list(cancer_type, ids, study_id, output_directory):
    """
    Writes one oncotree case list. Python verisons below
    3.6 will sort the dictionary keys which causes tests to fail

    Args:
        cancer_type: Oncotree code cancer type
        ids: GENIE sample ids
        study_id: cBioPortal study id
        output_directory: case list output directory

    Returns:
        case list file path
    """
    cancer_type = "NA" if cancer_type == "" else cancer_type
    cancer_type_no_spaces = (
        cancer_type.replace(" ", "_").replace(",", "").replace("/", "_")
    )
    cancer_type_no_spaces = (
        "no_oncotree_code" if cancer_type_no_spaces == "NA" else cancer_type_no_spaces
    )
    case_list_text = CASE_LIST_TEXT_TEMPLATE.format(
        study_id=study_id,
        stable_id=study_id + "_" + cancer_type_no_spaces,
        case_list_name="Tumor Type: " + cancer_type,
        case_list_description="All tumors with cancer type " + cancer_type,
        case_list_ids="\t".join(ids),
    )
    case_list_path = os.path.abspath(
        os.path.join(output_directory, "cases_" + cancer_type_no_spaces + ".txt")
    )
    with open(case_list_path, "w") as case_list_file:
        case_list_file.write(case_list_text)
    return case_list_path


def write_case_list_files(clinical_file_map, output_directory, study_id):
    """
    Writes the cancer_type case list file to case_lists directory

    Args:
        clinical_file_map: cancer type to sample id mapping from
                           create_case_lists_map
        output_directory: Directory to write case lists
        study_id: cBioPortal study id

    Returns:
        list: oncotree code case list files
    """
    case_list_files = []
    for cancer_type, ids in clinical_file_map.items():
        case_list_path = _write_single_oncotree_case_list(
            cancer_type, ids, study_id, output_directory
        )
        case_list_files.append(case_list_path)
    return case_list_files


def create_sequenced_samples(seq_assay_map, assay_info_file_name):
    """
    Gets samples sequenced

    Args:
        seq_assay_map: dictionary containing lists of samples per seq_assay_id
        assay_info_file_name: Assay information name

    Returns:
        list of cna samples
    """
    with open(assay_info_file_name, "r") as assay_info_file:
        reader = csv.DictReader(assay_info_file, dialect="excel-tab")
        cna_samples = []
        for row in reader:
            if "cna" in row["alteration_types"]:
                cna_samples.extend(seq_assay_map[row["SEQ_ASSAY_ID"]])
    return cna_samples


def write_case_list_sequenced(clinical_samples, output_directory, study_id):
    """
    Writes the genie sequenced and all samples. Since all samples
    are sequenced, _all and _sequenced are equal

    Args:
        clinical_samples: List of clinical samples
        output_directory: Directory to write case lists
        study_id: cBioPortal study id

    Returns:
        list: case list sequenced and all
    """
    caselist_files = []
    case_list_ids = "\t".join(clinical_samples)
    case_sequenced_path = os.path.abspath(
        os.path.join(output_directory, "cases_sequenced.txt")
    )
    with open(case_sequenced_path, "w") as case_list_sequenced_file:
        case_list_file_text = CASE_LIST_TEXT_TEMPLATE.format(
            study_id=study_id,
            stable_id=study_id + "_sequenced",
            case_list_name="Sequenced Tumors",
            case_list_description="All sequenced samples",
            case_list_ids=case_list_ids,
        )
        case_list_sequenced_file.write(case_list_file_text)
    cases_all_path = os.path.abspath(os.path.join(output_directory, "cases_all.txt"))
    with open(cases_all_path, "w") as case_list_all_file:
        case_list_file_text = CASE_LIST_TEXT_TEMPLATE.format(
            study_id=study_id,
            stable_id=study_id + "_all",
            case_list_name="All samples",
            case_list_description="All samples",
            case_list_ids=case_list_ids,
        )
        case_list_all_file.write(case_list_file_text)
    caselist_files.extend([case_sequenced_path, cases_all_path])
    return caselist_files


def write_case_list_cna(cna_samples, output_directory, study_id):
    """
    Writes the cna sequenced samples

    Args:
        cna_samples: List of cna samples
        output_directory: Directory to write case lists
        study_id: cBioPortal study id

    Returns:
        cna caselist path
    """
    case_list_ids = "\t".join(cna_samples)
    cna_caselist_path = os.path.abspath(os.path.join(output_directory, "cases_cna.txt"))
    with open(cna_caselist_path, "w") as case_list_file:
        case_list_file_text = CASE_LIST_TEXT_TEMPLATE.format(
            study_id=study_id,
            stable_id=study_id + "_cna",
            case_list_name="Samples with CNA",
            case_list_description="Samples with CNA",
            case_list_ids=case_list_ids,
        )
        case_list_file.write(case_list_file_text)
    return cna_caselist_path


def write_case_list_cnaseq(cna_samples, output_directory, study_id):
    """
    writes both cna and mutation samples (Just _cna file for now)

    Args:
        cna_samples: List of cna samples
        output_directory: Directory to write case lists
        study_id: cBioPortal study id

    Returns:
        cnaseq path
    """
    case_list_ids = "\t".join(cna_samples)
    cnaseq_caselist_path = os.path.abspath(
        os.path.join(output_directory, "cases_cnaseq.txt")
    )
    with open(cnaseq_caselist_path, "w") as case_list_file:
        case_list_file_text = CASE_LIST_TEXT_TEMPLATE.format(
            study_id=study_id,
            stable_id=study_id + "_cnaseq",
            case_list_name="Samples with CNA and mutation",
            case_list_description="Samples with CNA and mutation",
            case_list_ids=case_list_ids,
        )
        case_list_file.write(case_list_file_text)
    return cnaseq_caselist_path


def main(clinical_file_name, assay_info_file_name, output_directory, study_id):
    """Gets clinical file and gene matrix file and processes it
    to obtain case list files

    Args:
        clinical_file_name: Clinical file path
        assay_info_file_name: Assay information name
        output_directory: Output directory of case list files
        study_id: cBioPortal study id
    """
    case_lists_map, seq_assay_map, clin_samples = create_case_lists_map(
        clinical_file_name
    )
    write_case_list_files(case_lists_map, output_directory, study_id)
    # create_sequenced_samples used to get the samples, but since the removal
    # of WES samples, we can no longer rely on the gene matrix file to grab
    # all sequenced samples, must use assay information file
    cna_samples = create_sequenced_samples(seq_assay_map, assay_info_file_name)
    write_case_list_sequenced(clin_samples, output_directory, study_id)
    write_case_list_cna(cna_samples, output_directory, study_id)
    write_case_list_cnaseq(cna_samples, output_directory, study_id)
