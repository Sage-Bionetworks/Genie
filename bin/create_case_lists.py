import argparse

from genie import create_case_lists


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creation of case lists') #pylint: disable=invalid-name

    parser.add_argument("clinical_file_name",
                        type=str,
                        help="Clinical file path")
    parser.add_argument("assay_info_file_name",
                        type=str,
                        help="gene matrix file path")
    parser.add_argument("output_dir",
                        type=str,
                        help="Output directory")
    parser.add_argument("study_id",
                        type=str,
                        help="Output directory")
    args = parser.parse_args() #pylint: disable=invalid-name

    create_case_lists.main(args.clinical_file_name,
                           args.assay_info_file_name,
                           args.output_dir,
                           args.study_id)
