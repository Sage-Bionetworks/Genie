from __future__ import absolute_import
from collections.abc import Sequence
import logging
import os

import pandas
import synapseclient

from .example_filetype_format import FileTypeFormat
from . import process_functions

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Workflow(FileTypeFormat):

    _fileType = "veoibd_workflow"

    _process_kwargs = ["databaseSynId"]

    _required_columns = ['assay_id', 'center']

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]).startswith("VEOIBD-" + self.center + "-ASSAY") and \
               os.path.basename(filePath[0]).endswith(".md")

    def process_steps(self, filePath, databaseSynId):
        logger.debug("Performing process_steps for {}".format(self._fileType))

        logger.debug("Storing file at {}".format(databaseSynId))
        self.syn.store(synapseclient.File(filePath, parent=databaseSynId,
                                          annotations=dict(center=self.center, 
                                                           fileType=self._fileType)),

                       forceVersion=False)
        return(filePath)

    def _get_assay_id(self, filePathList):
        """Get the assay ID from the file name.
        """
        if isinstance(filePathList, Sequence):
            filePathList = filePathList[0]
        
        return os.path.basename(filePathList).split(".md")[0]
