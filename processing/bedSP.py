import os
import bed
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class bedSP(bed.bed):

	_fileType = "bedSP"

	_process_kwargs = ["newPath", "parentId", "databaseSynId"]

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]).startswith("nonGENIE_%s-" % self.center) and os.path.basename(filePath[0]).endswith(".bed")