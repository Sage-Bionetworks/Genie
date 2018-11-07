#!/usr/bin/env python
import logging
logger = logging.getLogger("genie")
import synapseclient
import process_functions
import subprocess
import argparse
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))



def main():
    parser = argparse.ArgumentParser(description='validate/processing genie sites')
    parser.add_argument("process", choices=['vcf','maf','main','mafSP'], help='Process vcf, maf or the rest of the files')
    parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")
    parser.add_argument("--center", type=str, help="Center")
    parser.add_argument("--onlyValidate", action='store_true', help = "Only validate the files, don't process")
    parser.add_argument("--createNewMafDatabase", action='store_true', help = "Creates a new maf database")
    parser.add_argument("--testing", action='store_true', help = "Testing the infrastructure!")
    # DEFAULT ARGUMENTS
    parser.add_argument("--vcf2mafPath", type=str, help="Path to vcf2maf", default="~/vcf2maf-1.6.14")
    parser.add_argument("--vepPath", type=str, help="Path to VEP", default="~/vep")
    parser.add_argument("--vepData", type=str, help="Path to VEP data", default="~/.vep")
    #parser.add_argument("--oncotreeLink", type=str, default='http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2017_06_21', help="Link to oncotree code")

    args = parser.parse_args()
    syn = process_functions.synLogin(args)
    
    #Database/folder syn id mapping
    CENTER_MAPPING = syn.tableQuery('SELECT * FROM %s where inputSynId is not null and release is true' % process_functions.getDatabaseSynId(syn, "centerMapping", test=args.testing))
    CENTER_MAPPING_DF = CENTER_MAPPING.asDataFrame()

    #If a center if specified, check if center is part of center mapping table,
    #Then make sure only that center is validated
    if args.center is not None:
        assert args.center in CENTER_MAPPING_DF.center.tolist()
        centers = [args.center]
    else:
        centers = CENTER_MAPPING_DF.center

    inputToDatabaseScript = os.path.join(SCRIPT_DIR, "input_to_database.py")
    for center in centers:
        if args.onlyValidate:
            logger.info("VALIDATING %s" % center)
            extendCommands = ["--onlyValidate"]
            args.process = "vcf"
            logPath = os.path.join(SCRIPT_DIR, "%s_validation_log.txt" % center)
        else:
            logger.info("PROCESSING %s" % center)
            extendCommands = ['--deleteOld', '--vcf2mafPath', os.path.expanduser(args.vcf2mafPath), '--vepPath', os.path.expanduser(args.vepPath), '--vepData', os.path.expanduser(args.vepData)]
                               #'--oncotreeLink', args.oncotreeLink]
            logPath = os.path.join(SCRIPT_DIR, "%s_%s_log.txt" % (center, args.process))

        command = ["python", inputToDatabaseScript, center, args.process, "--pemFile", os.path.expanduser(args.pemFile)] if args.pemFile is not None else ["python", inputToDatabaseScript, center, args.process]
        if args.testing:
            extendCommands.append("--testing")
        if args.createNewMafDatabase:
            #Only need to create the mafdatabase once
            extendCommands.append('--createNewMafDatabase')
            args.createNewMafDatabase = False
        command.extend(extendCommands)
        logger.info("COMMAND EXECUTED: %s" % " ".join(command)) 
        output = subprocess.check_output(command,stderr= subprocess.STDOUT)
        #logger.info(output)
        with open(logPath, "w") as centerLog:
            centerLog.write(output.decode("utf-8"))
        syn.store(synapseclient.File(logPath, parentId="syn10155804"))
        os.remove(logPath)
    #Only write out invalid reasons if the center isnt specified and if only validate
    if args.center is None and args.onlyValidate:
        logger.info("WRITING INVALID REASONS TO CENTER STAGING DIRS")
        writeInvalidReasonsScript = os.path.join(SCRIPT_DIR, 'writeInvalidReasons.py')
        writeInvalidReasons = ['python', writeInvalidReasonsScript, '--pemFile', os.path.expanduser(args.pemFile)] if args.pemFile is not None else ['python', writeInvalidReasonsScript]
        output = subprocess.check_call(writeInvalidReasons)

if __name__ == "__main__":
    main()
   