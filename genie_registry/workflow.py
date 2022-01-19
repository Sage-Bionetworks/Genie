import logging
import os

import synapseclient

from genie.example_filetype_format import FileTypeFormat

logger = logging.getLogger(__name__)


class workflow(FileTypeFormat):

    _fileType = "md"

    _process_kwargs = ["databaseSynId"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]).startswith(
            self.center
        ) and os.path.basename(filePath[0]).endswith(".md")

    def process_steps(self, filePath, databaseSynId):
        self.syn.store(synapseclient.File(filePath, parent=databaseSynId))
        return filePath
