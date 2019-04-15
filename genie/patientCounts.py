from __future__ import absolute_import
from genie import FileTypeFormat, process_functions
import logging
import os
import pandas as pd
logger = logging.getLogger(__name__)


class patientCounts(FileTypeFormat):

    _fileType = "patientCounts"

    _process_kwargs = ["newPath", "oncotreeLink", "databaseSynId"]

    _validation_kwargs = ["oncotreeLink"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "patient_counts.txt"

    def _process(self, patientCountsDf, oncotreeLink):
        patientCountsDf['CENTER'] = self.center

        oncotree_mapping_dict = \
            process_functions.get_oncotree_code_mappings(oncotreeLink)
        patientCountsDf['PRIMARY_CODE'] = [
            oncotree_mapping_dict[i.upper()]['ONCOTREE_PRIMARY_NODE']
            for i in patientCountsDf.ONCOTREE_CODE]
        return(patientCountsDf)

    def process_steps(
            self, patientCountsDf, newPath, oncotreeLink, databaseSynId):
        patientCountsDf = self._process(patientCountsDf, oncotreeLink)
        process_functions.updateData(
            self.syn, databaseSynId, patientCountsDf, self.center)
        patientCountsDf.to_csv(newPath, sep="\t",index=False)
        return(newPath)

    def _validate(self, patCountsDf, oncotreeLink):
        total_error = ""
        warning = ""
        # oncotree_mapping = process_functions.get_oncotree_codes(oncotreeLink)
        # if oncotree_mapping.empty:
        oncotree_mapping = pd.DataFrame()
        oncotree_mapping_dict = \
            process_functions.get_oncotree_code_mappings(oncotreeLink)
        oncotree_mapping['ONCOTREE_CODE'] = oncotree_mapping_dict.keys()
        haveColumn = process_functions.checkColExist(patCountsDf, "ONCOTREE_CODE")
        if haveColumn:
            if sum(patCountsDf['ONCOTREE_CODE'].duplicated()) > 0:
                total_error += "Patient Counts: Must not have any duplicated ONCOTREE CODES.\n"
            if not all(patCountsDf['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])):
                unmapped_oncotrees = patCountsDf['ONCOTREE_CODE'][~patCountsDf['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])]
                total_error += "Patient Counts: Please double check that all your ONCOTREE CODES exist in the mapping. You have %d codes that don't map. These are the codes that don't map: %s\n" % (len(unmapped_oncotrees),",".join(set(unmapped_oncotrees)))
        else:
            total_error += "Patient Counts: File must have ONCOTREE_CODE column.\n"

        haveColumn = process_functions.checkColExist(patCountsDf, "NUM_PATIENTS_PD1_PDL1")
        if haveColumn:
            if not all([isinstance(i, int) for i in patCountsDf['NUM_PATIENTS_PD1_PDL1']]):
                total_error += "Patient Counts: Must not have any null values, and must be all integers.\n"
        else:
            total_error += "Patient Counts: File must have NUM_PATIENTS_PD1_PDL1 column.\n"

        return(total_error, warning)
