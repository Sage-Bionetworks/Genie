import logging
import os
import subprocess

import pandas as pd

from .maf import maf
from . import process_functions

logger = logging.getLogger(__name__)


def contains_whitespace(row):
    """
    Helper function for validateVCF.
    No whitespace is allowed in VCF files

    Args:
        row: Each pandas dataframe row of a vcf

    Returns:
        Sum of the the amount of whitespace in a string
    """
    return(sum([" " in i for i in row if isinstance(i, str)]))


class vcf(maf):
    _fileType = "vcf"

    _process_kwargs = [
        "validVCF", "processing", "path_to_GENIE", "databaseToSynIdMappingDf",
        "vcf2mafPath", "veppath", "vepdata", "reference"]

    def _validateFilename(self, filePath):
        basename = os.path.basename(filePath[0])
        startswith_genie = basename.startswith("GENIE-{}-".format(self.center))
        endswith_vcf = basename.endswith(".vcf")
        assert startswith_genie and endswith_vcf

    def _get_dataframe(self, filePathList):
        headers = None
        filepath = filePathList[0]
        with open(filepath, "r") as vcffile:
            for row in vcffile:
                if row.startswith("#CHROM"):
                    headers = \
                        row.replace("\n", "").replace("\r", "").split("\t")
        if headers is not None:
            vcf = pd.read_csv(
                filepath, sep="\t", comment="#", header=None, names=headers)
        else:
            raise ValueError("Your vcf must start with the header #CHROM")
        return(vcf)

    def process_helper(self, vcffiles, path_to_GENIE, mafSynId, centerMafSynId,
                       vcf2mafPath, veppath, vepdata,
                       reference=None):
        logger.info('VCF2MAF %s' % self.center)
        centerInputFolder = os.path.join(path_to_GENIE, self.center, "input")
        centerStagingFolder = \
            os.path.join(path_to_GENIE, self.center, "staging")
        mafFiles = []
        for path in vcffiles:
            vcfName = os.path.basename(path)
            logger.info(vcfName)
            newVCFPath = os.path.join(centerInputFolder, vcfName)
            # remove chr from each row
            command = ["sed", "'s/^chr//'", path, ">", newVCFPath]
            subprocess.check_call(" ".join(command), shell=True)
            # Empty spaces must be replaced with a period
            command = ["sed", '-i', "'s/\t\t/\t.\t/g'", newVCFPath]
            subprocess.check_call(" ".join(command), shell=True)
            # All INFO/HGVS values have a whitespace,
            # which is not allowed in VCF specs. Replace that with a comma
            command = ['sed', '-i', "'s/ p\./,p./'", newVCFPath]
            subprocess.check_call(" ".join(command), shell=True)
            # Remove empty white space lines (VCF2MAF fails)
            command = ['sed', '-i', "'/^\s*$/d'", newVCFPath]
            subprocess.check_call(" ".join(command), shell=True)
            # Strips out windows indentations \r
            command = ['dos2unix', newVCFPath]
            subprocess.check_call(command)
            vcfCols = ["#CHROM", "POS", "ID", "REF", "ALT",
                       "QUAL", "FILTER", "INFO", "FORMAT"]
            with open(newVCFPath, "r") as f:
                for line in f:
                    if line.startswith("#CHROM"):
                        cols = line

            cols = cols.replace("\n", "")
            cols = cols.replace("\r", "")
            cols = cols.split("\t")

            samples = [i for i in cols if i not in vcfCols]

            # tumorName = vcfName.replace(".vcf","")

            if len(samples) == 1:
                tumor = samples[0]
                normal = "NORMAL"
            elif len(samples) == 2:
                # Assumes that Tumor is always first, normal is second
                tumor = samples[0]
                normal = samples[1]
            else:
                tumor = "TUMOR"
                normal = "NORMAL"
            # If the tumor name isn't TUMOR,
            # set the sample id to be the tumor name
            if tumor == "TUMOR":
                tumorName = vcfName.replace(".vcf", "")
            else:
                tumorName = tumor
            newMAFPath = newVCFPath + ".maf"
            if os.path.isfile(newMAFPath):
                mafFiles.append(newMAFPath)
            else:
                command = ['perl', os.path.join(vcf2mafPath, 'vcf2maf.pl'),
                           '--input-vcf', newVCFPath,
                           '--output-maf', newMAFPath,
                           '--vep-path', veppath,
                           '--vep-data', vepdata,
                           '--vep-forks', '8',
                           '--tumor-id', tumorName,
                           '--normal-id', normal,
                           '--vcf-tumor-id', tumor,
                           '--custom-enst', os.path.join(vcf2mafPath, 'data/isoform_overrides_uniprot')]
                if reference is not None:
                    command.extend(["--ref-fasta", reference])
                subprocess.check_call(command)
                if (os.path.isfile(newMAFPath)):
                    mafFiles.append(newMAFPath)

        logger.info("MERGING MAFS")
        # maf = pd.DataFrame()
        newMafPath = os.path.join(
            centerStagingFolder,
            "data_mutations_extended_{}_VCF.txt".format(self.center))
        narrowMafPath = os.path.join(
            centerStagingFolder,
            "data_mutations_extended_{}_VCF_narrow.txt".format(self.center))
        open(newMafPath, "w").close()
        open(narrowMafPath, "w").close()
        narrowMafColumns = [
            col['name'] for col in self.syn.getTableColumns(mafSynId)
            if col['name'] != 'inBED']

        narrowMafPaths = [narrowMafPath]
        for mafFile in mafFiles:
            mafDf = pd.read_csv(mafFile, sep="\t", comment="#")
            mafDf = self.formatMAF(mafDf)
            self.createFinalMaf(mafDf, newMafPath)
            narrowMafDf = mafDf[narrowMafColumns]
            maf_statinfo = os.stat(narrowMafPath)
            # This is done because of the table update size limit.
            # Can only update 1 GB at a time
            if maf_statinfo.st_size / 1000000 > 900:
                narrowMafPath = narrowMafPath.replace(
                    ".txt", str(len(narrowMafPaths)) + ".txt")
                narrowMafPaths.append(narrowMafPath)
            self.createFinalMaf(narrowMafDf, narrowMafPath)

        if len(mafFiles) > 0:
            # These functions have to be next to each other,
            # because no modifications can happen
            # Store MAF flat file into synapse
            # First because of 1 GB limit
            self.storeProcessedMaf(newMafPath, mafSynId, centerMafSynId)
            # Store Narrow MAF into db
            # Need to figure out how to fix this 1 GB limit
            for narrow_path in narrowMafPaths:
                self.storeProcessedMaf(
                    narrow_path, mafSynId, centerMafSynId, isNarrow=True)

        return(newMafPath)

    def process_steps(self, filePath, processing, databaseToSynIdMappingDf,
                      vcf2mafPath, veppath, vepdata, validVCF, path_to_GENIE,
                      reference):
        mutationFiles = []
        if processing == self._fileType:
            mafProcessing = "mafSP" if self._fileType == "mafSP" else 'vcf2maf'
            mafSynId = databaseToSynIdMappingDf.Id[
                databaseToSynIdMappingDf['Database'] == mafProcessing][0]
            centerMafSynId = databaseToSynIdMappingDf.Id[
                databaseToSynIdMappingDf['Database'] == "centerMaf"][0]
            logger.info(validVCF)
            vcfFilePath = self.process_helper(
                validVCF, path_to_GENIE, mafSynId, centerMafSynId,
                vcf2mafPath, veppath, vepdata, reference=reference)
            mutationFiles = [vcfFilePath]
            logger.info("UPDATED DATABASE WITH: %s" % ", ".join(mutationFiles))
        else:
            logger.info(
                "Please run with `--process {filetype}` parameter "
                "if you want to reannotate the {filetype} files".format(
                    filetype=self._fileType))
        return(mutationFiles)

    def _validate(self, vcf):
        '''
        Validates the content of a vcf file

        Args:
            vcf: pandas dataframe containing vcf content

        Returns:
            total_error - error messages
            warning - warning messages
        '''
        REQUIRED_HEADERS = pd.Series([
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])

        total_error = ""
        warning = ""
        if not all(REQUIRED_HEADERS.isin(vcf.columns)):
            total_error += (
                "Your vcf file must have these headers: "
                "CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.\n")

        if len(vcf.columns) > 8:
            if "FORMAT" not in vcf.columns:
                total_error += (
                    "Your vcf file must have FORMAT header "
                    "if genotype columns exist.\n")

        # Require that they report variants mapped to
        # either GRCh37 or hg19 without
        # the chr-prefix. variants on chrM are not supported
        haveColumn = process_functions.checkColExist(vcf, "#CHROM")
        if haveColumn:
            nochr = ["chr" in i for i in vcf['#CHROM'] if isinstance(i, str)]
            if sum(nochr) > 0:
                warning += (
                    "Your vcf file should not have the chr prefix "
                    "in front of chromosomes.\n")
            if sum(vcf['#CHROM'].isin(["chrM"])) > 0:
                total_error += \
                    "Your vcf file must not have variants on chrM.\n"

        # No white spaces
        temp = vcf.apply(lambda x: contains_whitespace(x), axis=1)
        if sum(temp) > 0:
            warning += (
                "Your vcf file should not have any "
                "white spaces in any of the columns.\n")
        # I can also recommend a `bcftools query` command that
        # will parse a VCF in a detailed way,
        # and output with warnings or errors if the format is not adhered too
        return(total_error, warning)
