from __future__ import absolute_import
from genie import maf, process_functions
import os
import logging
import pandas as pd
logger = logging.getLogger(__name__)

class mafSP(maf):

	_fileType = "mafSP"

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "nonGENIE_data_mutations_extended_%s.txt" % self.center

	def storeProcessedMaf(self, filePath, mafSynId, centerMafSynId, isNarrow=False):
		logger.info('STORING %s' % filePath)
		database = self.syn.get(mafSynId)
		mafDataFrame = pd.read_csv(filePath,sep="\t")
		process_functions.updateData(self.syn, mafSynId, mafDataFrame, self.center, filterByColumn="Center", toDelete=True)
		return(filePath)