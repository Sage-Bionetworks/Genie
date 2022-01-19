from .bed import bed
import os
import logging

logger = logging.getLogger(__name__)


class bedSP(bed):

    _fileType = "bedSP"

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]).startswith(
            "nonGENIE_%s-" % self.center
        ) and os.path.basename(filePath[0]).endswith(".bed")
