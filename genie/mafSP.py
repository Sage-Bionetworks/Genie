from __future__ import absolute_import
from genie import maf
from genie import process_functions

import os
import logging
import pandas as pd
#logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def updateData(syn, databaseSynId, newData, center, col, toDelete=False):
	databaseEnt = syn.get(databaseSynId)
	database = syn.tableQuery("SELECT * FROM %s where Center ='%s'" % (databaseSynId, center))
	database = database.asDataFrame()
	process_functions.updateDatabase(syn, database, newData, databaseSynId, databaseEnt.primaryKey, toDelete)
	
class mafSP(maf.maf):

	_fileType = "mafSP"

	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "nonGENIE_data_mutations_extended_%s.txt" % self.center

	def validate_steps(self, filePathList, **kwargs):
		logger.info("VALIDATING %s" % os.path.basename(filePathList[0]))
		mutationDF = pd.read_csv(filePathList[0],sep="\t",comment="#",na_values = ['-1.#IND', '1.#QNAN', '1.#IND', 
								 '-1.#QNAN', '#N/A N/A', '#N/A', 'N/A', '#NA', 'NULL', 'NaN', 
								 '-NaN', 'nan','-nan',''],keep_default_na=False)
		total_error, warning = self.validate_helper(mutationDF,SP=True)
		return(total_error, warning)


	def storeProcessedMaf(self, filePath, mafSynId, centerMafSynId, isNarrow=False):
		logger.info('STORING %s' % filePath)
		database = self.syn.get(mafSynId)
		mafDataFrame = pd.read_csv(filePath,sep="\t")
		updateData(self.syn, mafSynId, mafDataFrame, self.center, database.primaryKey, toDelete=True)
		return(filePath)