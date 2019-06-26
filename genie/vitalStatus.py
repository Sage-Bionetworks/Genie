import datetime
import logging
import os

import pandas as pd

from .example_filetype_format import FileTypeFormat
from . import process_functions

logger = logging.getLogger(__name__)


class vitalStatus(FileTypeFormat):

    _fileType = "vitalStatus"

    def _validateFilename(self, filePath):
        '''
        Validate filename: vital_status.txt
        '''
        assert os.path.basename(filePath[0]) == "vital_status.txt"

    def _validate(self, vitalStatusDf):
        total_error = ""
        warning = ""

        # PATIENT ID
        haveColumn = \
            process_functions.checkColExist(vitalStatusDf, "PATIENT_ID")
        if haveColumn:
            if vitalStatusDf.PATIENT_ID.isnull().any():
                total_error += (
                    "Vital status file: Please double check your PATIENT_ID"
                    " column. No null values allowed.\n")
        else:
            total_error += "Vital status file: Must have PATIENT_ID column.\n"

        # YEAR DEATH
        haveColumn = \
            process_functions.checkColExist(vitalStatusDf, "YEAR_DEATH")
        if haveColumn:
            notNullYears = \
                vitalStatusDf.YEAR_DEATH[~vitalStatusDf.YEAR_DEATH.isnull()]
            try:
                notNullYears.apply(
                    lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
            except Exception:
                total_error += (
                    "Vital status file: Please double check your YEAR_DEATH"
                    " column, it must be an integer in YYYY format or an "
                    "empty string.\n")
        else:
            total_error += "Vital status file: Must have YEAR_DEATH column.\n"

        # YEAR CONTACT
        haveColumn = \
            process_functions.checkColExist(vitalStatusDf, "YEAR_CONTACT")
        if haveColumn:
            notNullYears = vitalStatusDf.YEAR_CONTACT[
                ~vitalStatusDf.YEAR_CONTACT.isnull()]
            try:
                notNullYears.apply(
                    lambda x: datetime.datetime.strptime(str(int(x)), '%Y'))
            except Exception:
                total_error += (
                    "Vital status file: Please double check your "
                    "YEAR_CONTACT column, it must be an integer in YYYY "
                    "format or an empty string.\n")
        else:
            total_error += \
                "Vital status file: Must have YEAR_CONTACT column.\n"

        # INT CONTACT
        haveColumn = \
            process_functions.checkColExist(vitalStatusDf, "INT_CONTACT")
        if haveColumn:
            check_ints = [
                process_functions.checkInt(i)
                for i in vitalStatusDf.INT_CONTACT
                if not pd.isnull(i) and i not in ['>32485', '<6570']]
            if not all(check_ints):
                total_error += (
                    "Vital status file: Please double check your "
                    "INT_CONTACT column, it must be an integer, an empty "
                    "string, >32485, or <6570.\n")
        else:
            total_error += "Vital status file: Must have INT_CONTACT column.\n"

        # INT DOD
        haveColumn = process_functions.checkColExist(vitalStatusDf, "INT_DOD")
        if haveColumn:
            check_dod_ints = [
                process_functions.checkInt(i)
                for i in vitalStatusDf.INT_DOD
                if not pd.isnull(i) and i not in ['>32485', '<6570']]
            if not all(check_dod_ints):
                total_error += (
                    "Vital status file: Please double check your "
                    "INT_DOD column, it must be an integer, an empty "
                    "string, >32485, or <6570.\n")
        else:
            total_error += "Vital status file: Must have INT_DOD column.\n"

        haveColumn = process_functions.checkColExist(vitalStatusDf, "DEAD")
        if haveColumn:
            check_boolean = [
                isinstance(i, bool)
                for i in vitalStatusDf.DEAD if not pd.isnull(i)]
            if not all(check_boolean):
                total_error += (
                    "Vital status file: Please double check your DEAD "
                    "column, it must be a boolean value or an empty string.\n")
        else:
            total_error += "Vital status file: Must have DEAD column.\n"

        return(total_error, warning)

    def _process(self, vitalStatusDf):
        vitalStatusDf.PATIENT_ID = [
            process_functions.checkGenieId(patient, self.center)
            for patient in vitalStatusDf.PATIENT_ID]
        vitalStatusDf['CENTER'] = self.center
        return(vitalStatusDf)

    # PROCESS
    def process_steps(self, vitalStatusDf, databaseSynId, newPath):
        vitalStatusDf = self._process(vitalStatusDf)
        process_functions.updateData(
            self.syn, databaseSynId, vitalStatusDf, self.center)
        vitalStatusDf.to_csv(newPath, sep="\t", index=False)
        return(newPath)
