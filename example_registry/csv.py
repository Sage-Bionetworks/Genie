import logging
import os

import synapseclient

from synapsegenie.example_filetype_format import FileTypeFormat
from synapsegenie import process_functions

logger = logging.getLogger(__name__)


class Csv(FileTypeFormat):

    _fileType = "csv"

    _process_kwargs = ["databaseSynId"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]).endswith(".csv")

    def _process(self, df):
        df.columns = [df.upper() for col in df.columns]
        return df

    def process_steps(self, df, newPath, databaseSynId):
        df = self._process(df)
        process_functions.updateData(self.syn, databaseSynId, df, self.center,
                                     toDelete=True)
        df.to_csv(newPath, sep="\t", index=False)
        return newPath

    def _validate(self, df):
        total_error = ""
        warning = ""
        if df.empty:
            total_error += "{}: File must not be empty".format(self._fileType)
        return total_error, warning
