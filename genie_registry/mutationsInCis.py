import logging
import os

import pandas as pd

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions

logger = logging.getLogger(__name__)

# def updateMutationInCisData(syn, databaseSynId, newData, center, col, toDelete=False):
#     databaseEnt = syn.get(databaseSynId)
#     database = syn.tableQuery("SELECT * FROM %s where Center ='%s'" % (databaseSynId, center))
#     database = database.asDataFrame()[col]
#     process_functions.updateDatabase(syn, database, newData, databaseSynId, databaseEnt.primaryKey, toDelete)


class mutationsInCis(FileTypeFormat):

    _fileType = "mutationsInCis"

    _validation_kwargs = []

    def _get_dataframe(self, filePathList):
        """
        Mutation In Cis is a csv file
        """
        filePath = filePathList[0]
        df = pd.read_csv(filePath, comment="#")
        return df

    # VALIDATE FILENAME
    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "mutationsInCis_filtered_samples.csv"

    # PROCESS
    def process_steps(self, mutationInCis, newPath, databaseSynId):
        process_functions.updateData(
            self.syn, databaseSynId, mutationInCis, self.center, filterByColumn="Center"
        )
        mutationInCis.to_csv(newPath, sep="\t", index=False)
        return newPath

    def _validate(self, mutationInCisDf):
        mutationInCisSynId = self.genie_config["mutationsInCis"]
        # Pull down the correct database
        existingMergeCheck = self.syn.tableQuery(
            "select * from {} where Center = '{}'".format(
                mutationInCisSynId, self.center
            )
        )
        existingMergeCheckDf = existingMergeCheck.asDataFrame()

        total_error = ""
        warning = ""
        required_headers = pd.Series(
            [
                "Flag",
                "Center",
                "Tumor_Sample_Barcode",
                "Hugo_Symbol",
                "HGVSp_Short",
                "Variant_Classification",
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
                "t_alt_count_num",
                "t_depth",
            ]
        )
        primaryKeys = [
            "Tumor_Sample_Barcode",
            "HGVSp_Short",
            "Start_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
        ]
        if not all(required_headers.isin(mutationInCisDf.columns)):
            missing_headers = required_headers[
                ~required_headers.isin(mutationInCisDf.columns)
            ]
            total_error += (
                "Mutations In Cis Filter File: "
                "Must at least have these headers: %s.\n" % ",".join(missing_headers)
            )
        else:
            new = mutationInCisDf[primaryKeys].fillna("")
            existing = existingMergeCheckDf[primaryKeys].fillna("")

            existing["primaryAll"] = [
                " ".join(values.astype(str)) for i, values in existing.iterrows()
            ]
            new["primaryAll"] = [
                " ".join(values.astype(str)) for i, values in new.iterrows()
            ]
            if not all(new.primaryAll.isin(existing.primaryAll)):
                total_error += (
                    "Mutations In Cis Filter File: "
                    "All variants must come from the original "
                    "mutationInCis_filtered_samples.csv file in "
                    "each institution's staging folder.\n"
                )

        if process_functions.checkColExist(mutationInCisDf, "Tumor_Sample_Barcode"):
            error = process_functions.validate_genie_identifier(
                identifiers=mutationInCisDf["Tumor_Sample_Barcode"],
                center=self.center,
                filename="Mutations In Cis Filter File",
                col="TUMOR_SAMPLE_BARCODE",
            )
            total_error += error

        return total_error, warning
