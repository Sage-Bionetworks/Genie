import os
import logging
import maf
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class mafSP(maf.maf):

	_fileType = "mafSP"

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "nonGENIE_data_mutations_extended_%s.txt" % self.center

	def validate_steps(self, filePathList, **kwargs):
		total_error, warning = self.validate_helper(filePathList[0],SP=True)
		return(total_error, warning)
