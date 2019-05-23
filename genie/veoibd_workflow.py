from __future__ import absolute_import
from genie import FileTypeFormat
import os
import logging
import synapseclient
logger = logging.getLogger(__name__)


class Workflow(FileTypeFormat):

    _fileType = "veoibd_workflow"

    _process_kwargs = ["databaseSynId"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]).startswith("VEOIBD-" + self.center + "-ASSAY") and \
               os.path.basename(filePath[0]).endswith(".md")

    def process_steps(self, filePath, databaseSynId):
        self.syn.store(synapseclient.File(filePath, parent=databaseSynId))
        return(filePath)
