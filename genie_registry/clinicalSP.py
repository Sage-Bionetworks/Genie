import logging
import os

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions

logger = logging.getLogger(__name__)


class clinicalSP(FileTypeFormat):

    _fileType = "clinicalSP"

    # VALIDATE FILENAME
    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "nonGENIE_data_clinical.txt"

    def _process(self, clinicalSPDf):
        clinicalSPDf["SAMPLE_ID"] = [
            process_functions.checkGenieId(sample, self.center)
            for sample in clinicalSPDf["SAMPLE_ID"]
        ]
        clinicalSPDf["CENTER"] = self.center
        clinicalSPDf["PATIENT_ID"] = [
            process_functions.checkGenieId(sample, self.center)
            for sample in clinicalSPDf["PATIENT_ID"]
        ]
        return clinicalSPDf

    def process_steps(self, clinicalSPDf, newPath, databaseSynId):
        clinicalSPDf = self._process(clinicalSPDf)
        process_functions.updateData(self.syn, databaseSynId, clinicalSPDf, self.center)
        clinicalSPDf.to_csv(newPath, sep="\t", index=False)
        return newPath

    def _validate(self, clinicalDF):
        clinicalDF.columns = [col.upper() for col in clinicalDF.columns]
        # clinicalDF = clinicalDF.fillna("")
        total_error = ""
        warning = ""

        # CHECK: SAMPLE_ID
        haveColumn = process_functions.checkColExist(clinicalDF, "SAMPLE_ID")
        if not haveColumn:
            total_error += (
                "nonGENIE_data_clinical.txt: " "File must have SAMPLE_ID column.\n"
            )
        else:
            if sum(clinicalDF["SAMPLE_ID"].isnull()) > 0:
                total_error += (
                    "nonGENIE_data_clinical.txt: "
                    "There can't be any blank values for SAMPLE_ID\n"
                )

        # CHECK: SEQ_ASSAY_ID
        haveColumn = process_functions.checkColExist(clinicalDF, "SEQ_ASSAY_ID")
        if haveColumn:
            if sum(clinicalDF["SEQ_ASSAY_ID"].isnull()) > 0:
                warning += (
                    "nonGENIE_data_clinical.txt: Please double check "
                    "your SEQ_ASSAY_ID columns, there are empty rows.\n"
                )
        else:
            total_error += (
                "nonGENIE_data_clinical.txt: " "File must have SEQ_ASSAY_ID column.\n"
            )

        # CHECK: PATIENT_ID
        haveColumn = process_functions.checkColExist(clinicalDF, "PATIENT_ID")
        if not haveColumn:
            total_error += (
                "nonGENIE_data_clinical.txt: " "File must have PATIENT_ID column.\n"
            )
        else:
            if sum(clinicalDF["PATIENT_ID"].isnull()) > 0:
                total_error += (
                    "nonGENIE_data_clinical.txt: "
                    "There can't be any blank values for PATIENT_ID\n"
                )
        return (total_error, warning)
