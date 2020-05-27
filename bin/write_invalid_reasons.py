"""write invalid reasons cli"""
import argparse

from genie import process_functions, write_invalid_reasons


def main():
    parser = argparse.ArgumentParser(description='Write invalid reasons')
    parser.add_argument("--pemFile", type=str,
                        help="Path to PEM file (genie.pem)")
    parser.add_argument("--debug", action='store_true',
                        help="Synapse Debug Feature")
    args = parser.parse_args()

    syn = process_functions.synLogin(args.pemFile, debug=args.debug)
    center_mapping = syn.tableQuery(
        'SELECT * FROM syn10061452 where inputSynId is not null '
        'and release is true'
    )
    center_mappingdf = center_mapping.asDataFrame()
    error_tracker_synid = "syn10153306"
    write_invalid_reasons.write(syn, center_mappingdf, error_tracker_synid)


if __name__ == "__main__":
    main()
