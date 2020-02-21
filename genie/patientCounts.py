import logging
import os

import pandas as pd

from .example_filetype_format import FileTypeFormat
from . import process_functions

logger = logging.getLogger(__name__)


class patientCounts(FileTypeFormat):
    '''
    Deprecated class
    '''
    _fileType = "patientCounts"

    _process_kwargs = ["newPath", "oncotree_link", "databaseSynId"]

    _validation_kwargs = ["oncotree_link"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "patient_counts.txt"

    def _process(self, patientCountsDf, oncotree_link):
        patientCountsDf['CENTER'] = self.center
        oncotree_mapping_dict = \
            process_functions.get_oncotree_code_mappings(oncotree_link)
        patientCountsDf['PRIMARY_CODE'] = [
            oncotree_mapping_dict[i.upper()]['ONCOTREE_PRIMARY_NODE']
            for i in patientCountsDf.ONCOTREE_CODE]

        return(patientCountsDf)

    def process_steps(
            self, patientCountsDf, newPath, oncotree_link, databaseSynId):
        patientCountsDf = self._process(patientCountsDf, oncotree_link)
        process_functions.updateData(
            self.syn, databaseSynId, patientCountsDf, self.center)
        patientCountsDf.to_csv(newPath, sep="\t", index=False)
        return(newPath)

    def _validate(self, patCountsDf, oncotree_link):
        total_error = ""
        warning = ""
        # oncotree_mapping = process_functions.get_oncotree_codes(oncotree_link)
        # if oncotree_mapping.empty:
        oncotree_mapping = pd.DataFrame()
        oncotree_mapping_dict = \
            process_functions.get_oncotree_code_mappings(oncotree_link)
        oncotree_mapping['ONCOTREE_CODE'] = oncotree_mapping_dict.keys()

        haveColumn = \
            process_functions.checkColExist(patCountsDf, "ONCOTREE_CODE")

        if haveColumn:
            if sum(patCountsDf['ONCOTREE_CODE'].duplicated()) > 0:
                total_error += (
                    "Patient Counts: "
                    "Must not have any duplicated ONCOTREE CODES.\n")
            if not all(patCountsDf['ONCOTREE_CODE'].isin(
                    oncotree_mapping['ONCOTREE_CODE'])):
                unmapped_oncotrees = patCountsDf['ONCOTREE_CODE'][
                    ~patCountsDf['ONCOTREE_CODE'].isin(
                        oncotree_mapping['ONCOTREE_CODE'])]
                total_error += (
                    "Patient Counts: Please double check that all your "
                    "ONCOTREE CODES exist in the mapping. You have {} codes "
                    "that don't map. These are the codes that "
                    "don't map: {}\n".format(
                        len(unmapped_oncotrees),
                        ",".join(set(unmapped_oncotrees))))
        else:
            total_error += (
                "Patient Counts: File must have ONCOTREE_CODE column.\n")

        haveColumn = process_functions.checkColExist(
            patCountsDf, "NUM_PATIENTS_PD1_PDL1")

        if haveColumn:
            if not all([isinstance(i, int)
                        for i in patCountsDf['NUM_PATIENTS_PD1_PDL1']]):
                total_error += (
                    "Patient Counts: Must not have any null values, "
                    "and must be all integers.\n")
        else:
            total_error += (
                "Patient Counts: File must have "
                "NUM_PATIENTS_PD1_PDL1 column.\n")
        return(total_error, warning)
