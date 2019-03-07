from __future__ import absolute_import
from genie import example_filetype_format
import os
import logging
import synapseclient
logger = logging.getLogger(__name__)


class workflow(example_filetype_format.FileTypeFormat):

    _fileType = "md"

    _process_kwargs = ["databaseSynId"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]).startswith(self.center) and os.path.basename(filePath[0]).endswith(".md")

    def process_steps(self, filePath, *args, **kwargs):
        logger.info('PROCESSING %s' % filePath)
        databaseSynId = kwargs['databaseSynId']
        self.syn.store(synapseclient.File(filePath, parent=databaseSynId))
        return(filePath)
