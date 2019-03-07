#! /usr/bin/env python

# ---------------------------------------------------------------
#  Script to create case lists per cancer type
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# imports
import os
import sys
import csv
import argparse

from future.utils import iteritems
# or
from six import iteritems

# ---------------------------------------------------------------
# functions

# ---------------------------------------------------------------
# creates the case list dictionary
# key = cancer_type
# value = list of sids
def create_case_lists_map(clinical_file_name):
    clinical_file = open(clinical_file_name,'rU')
    clinical_file_map = {}
    reader = csv.DictReader(clinical_file,dialect='excel-tab')
    for row in reader:
        if row['CANCER_TYPE'] not in clinical_file_map:
            clinical_file_map[row['CANCER_TYPE']] = [row['SAMPLE_ID']]
        else:
            clinical_file_map[row['CANCER_TYPE']].append(row['SAMPLE_ID'])
    clinical_file.close()
    return(clinical_file_map)

# ---------------------------------------------------------------
# writes the file to case_lists directory inside the directory
def write_case_list_files(clinical_file_map, output_directory, study_id):
    for cancer_type,ids in iteritems(clinical_file_map):
        cancer_type_no_spaces = cancer_type.replace(' ','_').replace(',','').replace("/","_")
        case_list_file = open(os.path.abspath(output_directory + '/' + 'case_list_' + cancer_type_no_spaces + '.txt'),'w')
        stable_id = study_id + '_' + cancer_type_no_spaces
        case_list_name = 'Tumor Type: ' + cancer_type
        case_list_description = 'All tumors with cancer type ' + cancer_type
        case_list_ids = '\t'.join(ids)
        case_list_file.write('cancer_study_identifier: ' + study_id + '\n' +
                                'stable_id: ' + stable_id + '\n' + 
                                'case_list_name: ' + case_list_name + '\n' + 
                                'case_list_description: ' + case_list_description + '\n' +
                                'case_list_ids: ' + case_list_ids)
        case_list_file.close()

def create_sequenced_samples(gene_matrix_file_name):
    gene_matrix_file = open(gene_matrix_file_name,'r')
    reader = csv.DictReader(gene_matrix_file,dialect='excel-tab')
    clinical_samples = []
    cna_samples = []
    for row in reader:
        if row['cna'] != "NA":
            cna_samples.append(row['SAMPLE_ID'])
        clinical_samples.append(row['SAMPLE_ID'])

    return(clinical_samples, cna_samples)
# ---------------------------------------------------------------
# writes the genie sequenced samples
def write_case_list_sequenced(clinical_samples, output_directory, study_id):
    case_list_file = open(os.path.abspath(output_directory + '/' + 'case_list_sequenced.txt'),'w')
    stable_id = study_id + '_sequenced'
    case_list_name = 'Sequenced Tumors'
    case_list_description = 'All sequenced samples (%s)'
    case_list_ids = '\t'.join(clinical_samples)
    case_list_file.write('cancer_study_identifier: ' + study_id + '\n' +
                            'stable_id: ' + stable_id + '\n' + 
                            'case_list_name: ' + case_list_name + '\n' + 
                            'case_list_description: ' + case_list_description + '\n' +
                            'case_list_ids: ' + case_list_ids)
    case_list_file.close()

# ---------------------------------------------------------------
# writes the cna sequenced samples
def write_case_list_cna(cna_samples, output_directory, study_id):
    case_list_file = open(os.path.abspath(output_directory + '/' + 'case_list_cna.txt'),'w')
    stable_id = study_id + '_cna'
    case_list_name = 'Samples with CNA'
    case_list_description = 'Samples with CNA (%s)'
    case_list_ids = '\t'.join(cna_samples)
    case_list_file.write('cancer_study_identifier: ' + study_id + '\n' +
                            'stable_id: ' + stable_id + '\n' + 
                            'case_list_name: ' + case_list_name + '\n' + 
                            'case_list_description: ' + case_list_description + '\n' +
                            'case_list_ids: ' + case_list_ids)
    case_list_file.close()

# ---------------------------------------------------------------
# writes both cna and mutation samples (Just _cna file for now)
def write_case_list_cnaseq(cna_samples, output_directory, study_id):
    case_list_file = open(os.path.abspath(output_directory + '/' + 'case_list_cnaseq.txt'),'w')
    stable_id = study_id + '_cnaseq'
    case_list_name = 'Samples with CNA and mutation'
    case_list_description = 'Samples with CNA and mutation (%s)'
    case_list_ids = '\t'.join(cna_samples)
    case_list_file.write('cancer_study_identifier: ' + study_id + '\n' +
                            'stable_id: ' + stable_id + '\n' + 
                            'case_list_name: ' + case_list_name + '\n' + 
                            'case_list_description: ' + case_list_description + '\n' +
                            'case_list_ids: ' + case_list_ids)
    case_list_file.close()
# ---------------------------------------------------------------
# gets clin file and processes it 
def create_case_lists(clinical_file_name, gene_matrix_file_name, output_directory, study_id):
    case_lists_map = create_case_lists_map(clinical_file_name)
    write_case_list_files(case_lists_map, output_directory, study_id)
    clinical_samples, cna_samples = create_sequenced_samples(gene_matrix_file_name)
    write_case_list_sequenced(clinical_samples, output_directory, study_id)
    write_case_list_cna(cna_samples, output_directory, study_id)
    write_case_list_cnaseq(cna_samples, output_directory, study_id)

# ---------------------------------------------------------------
# the main
def main():
    parser = argparse.ArgumentParser(description='Creation of case lists')

    parser.add_argument("clinical_file_name", type=str,
                        help="Clinical file path")
    parser.add_argument("gene_matrix_file_name", type=str,
                        help="CNA file path")
    parser.add_argument("output_dir", type=str,
                        help="Output directory")
    parser.add_argument("study_id", type=str,
                        help="Output directory")
    args = parser.parse_args()

    clinical_file_name = args.clinical_file_name
    gene_matrix_file_name = args.gene_matrix_file_name
    output_directory = args.output_dir
    study_id = args.study_id

    create_case_lists(clinical_file_name, gene_matrix_file_name, output_directory, study_id)

# ---------------------------------------------------------------
# do a main
if __name__ == '__main__':
    main()
