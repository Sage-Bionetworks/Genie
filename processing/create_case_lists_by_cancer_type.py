#! /usr/bin/env python

# ---------------------------------------------------------------
#  Script to create case lists per cancer type
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# imports
import os
import sys
import getopt
import csv

# ---------------------------------------------------------------
# globals
ERROR_FILE = sys.stderr
OUTPUT_FILE = sys.stdout

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

	return clinical_file_map

# ---------------------------------------------------------------
# writes the file to case_lists directory inside the directory
def write_case_list_files(clinical_file_map, output_directory, study_id):
	for cancer_type,ids in clinical_file_map.iteritems():
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

# ---------------------------------------------------------------
# gets clin file and processes it 
def create_case_lists(clinical_file_name, output_directory, study_id):
    case_lists_map = create_case_lists_map(clinical_file_name)
    write_case_list_files(case_lists_map, output_directory, study_id)

# ---------------------------------------------------------------
# displays usage of program
def usage():
	print >> OUTPUT_FILE, 'create_case_lists_by_cancer_type.py --clinical-file <path/to/clinical/file> --output-directory <path/to/output/directory> study-id <cancer_study_identifier>'

# ---------------------------------------------------------------
# the main
def main():
	# parse command line
    try:
        opts,args = getopt.getopt(sys.argv[1:],'',['clinical-file=', 'output-directory=', 'study-id='])
    except getopt.error,msg:
        print >> ERROR_FILE,msg
        usage()
        sys.exit(2)

    clinical_file_name = ''
    output_directory = ''
    study_id = ''

	# process options
    for o, a in opts:
        if o == '--clinical-file':
            clinical_file_name = a
        elif o == '--output-directory':
        	output_directory = a
        elif o == '--study-id':
        	study_id = a
	
    if clinical_file_name == '' or output_directory == '' or study_id == '':
        usage()
        sys.exit(2)

	# check existence of file
    if not os.path.exists(os.path.abspath(clinical_file_name)):
        print >> ERROR_FILE, 'clinical file cannot be found: ' + clinical_file_name
        sys.exit(2)
    if not os.path.isdir(os.path.abspath(output_directory)):
    	print >> ERROR_FILE, 'directory cannot be found or is not a directory: ' + output_directory
        sys.exit(2)

    create_case_lists(clinical_file_name, output_directory, study_id)

# ---------------------------------------------------------------
# do a main
if __name__ == '__main__':
	main()
