from __future__ import absolute_import
from genie import example_filetype_format
from genie import process_functions

import os
import pandas as pd
import logging
from functools import partial
import synapseclient

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def validateSymbol(gene, bedDf, returnMappedDf=True):
	valid=False
	if sum(bedDf['Hugo_Symbol'] == gene) > 0:
		valid=True
	elif sum(bedDf['ID'] == gene) > 0:
		mismatch = bedDf[bedDf['ID'] == gene]
		mismatch.drop_duplicates(inplace=True)
		logger.info("%s will be remapped to %s" % (gene, mismatch['Hugo_Symbol'].values[0]))
		gene = mismatch['Hugo_Symbol'].values[0]
	else:
		logger.warning("%s cannot be remapped and will not be released. The symbol must exist in your seq assay ids (bed files) and must be mappable to a gene." % gene)
		gene = pd.np.nan
	if returnMappedDf:
		return(gene)
	else:
		return(valid)

def makeCNARow(row, symbols):
	totalrow = "%s\n%s" % (",".join(symbols),",".join(row.astype(str)))
	totalrow = totalrow.replace(".0","")
	return(totalrow)

def mergeCNAvalues(x):
	uniqueValues = set(x.tolist())
	if len(uniqueValues) == 1:
		return(x.tolist()[0])
	elif len(uniqueValues) == 2 and 0 in uniqueValues:
		uniqueValues.remove(0)
		return(list(uniqueValues)[0])
	else:
		return(pd.np.nan)

	
def checkIfOneZero(x):
	assert len(set(x.tolist())) == 1, "Can only be one unique value"

class cna(example_filetype_format.FileTypeFormat):
   
	_fileType = "cna"

	_process_kwargs = ["newPath", "databaseSynId",'test','databaseToSynIdMappingDf']

	_validation_kwargs = ['testing','noSymbolCheck']

	# VALIDATE FILENAME
	def _validateFilename(self, filePath):
		assert os.path.basename(filePath[0]) == "data_CNA_%s.txt" % self.center

	def _process(self, cnaDf, test=False):
		checkBy = "TUMOR_SAMPLE_BARCODE"

		cnaDf.rename(columns= {cnaDf.columns[0]:cnaDf.columns[0].upper()}, inplace=True)
		cnaDf.rename(columns= {"HUGO_SYMBOL":"Hugo_Symbol"}, inplace=True)

		columns = [col.upper() for col in cnaDf.columns]
		index = [i for i, col in enumerate(cnaDf.columns) if col.upper() == "ENTREZ_GENE_ID"]
		if len(index) > 0:
			del cnaDf[cnaDf.columns[index][0]]
		#validateSymbol = partial(process_functions.validateSymbol,returnMapping=True)
		#invalidated_genes = self.pool.map(validateSymbol, cna["HUGO_SYMBOL"].drop_duplicates())
		#cna, nonmapped = process_functions.remapGenes(invalidated_genes, cna, "HUGO_SYMBOL",isBedFile=True)
		bedSynId = process_functions.getDatabaseSynId(self.syn, "bed", test=test)
		bed = self.syn.tableQuery("select Hugo_Symbol, ID from %s where CENTER = '%s'" % (bedSynId, self.center))
		bedDf = bed.asDataFrame()
		#originalSymbols = cnaDf['HUGO_SYMBOL'].copy()
		cnaDf['Hugo_Symbol'] = cnaDf['Hugo_Symbol'].apply(lambda x: validateSymbol(x, bedDf))
		order = cnaDf.columns
		# unmappable = cnaDf[cnaDf['HUGO_SYMBOL'].isnull()]
		# unmappableSymbols = originalSymbols[cnaDf['HUGO_SYMBOL'].isnull()]

		cnaDf = cnaDf[~cnaDf['Hugo_Symbol'].isnull()]
		duplicatedGenes = pd.DataFrame()
		for i in cnaDf['Hugo_Symbol'][cnaDf['Hugo_Symbol'].duplicated()].unique():
			dups = cnaDf[cnaDf['Hugo_Symbol'] == i]
			newVal = dups[dups.columns[dups.columns!="HUGO_SYMBOL"]].apply(mergeCNAvalues)
			temp = pd.DataFrame(newVal).transpose()
			temp['Hugo_Symbol'] = i
			duplicatedGenes = duplicatedGenes.append(temp)
		cnaDf.drop_duplicates('Hugo_Symbol',keep=False, inplace=True)
		cnaDf = cnaDf.append(duplicatedGenes)
		cnaDf = cnaDf[order]
		#symbols = cnaDf['HUGO_SYMBOL']
		#del cnaDf['HUGO_SYMBOL']
		cnaDf = cnaDf.fillna('NA')
		newsamples = [process_functions.checkGenieId(i,self.center) if i != "Hugo_Symbol" else i for i in cnaDf.columns]
		#Transpose matrix
		# cnaDf = cnaDf.transpose()
		# data = cnaDf.apply(lambda row: makeCNARow(row, symbols), axis=1)

		#Transpose matrix
		# del unmappable['HUGO_SYMBOL']
		# unmappable = unmappable.transpose()
		# unmappableData = unmappable.apply(lambda row: makeCNARow(row, unmappableSymbols), axis=1)

		# newCNA = pd.DataFrame()
		# newCNA[checkBy] = newsamples
		# newCNA['CNAData'] = data.values
		# newCNA['CENTER'] = self.center
		# newCNA['unmappedData'] = unmappableData.values
		#newCNA = newCNA[~newCNA['CNAData'].isnull()]
		#remove the 0.0, 1.0 and 2.0
		# os.system("sed 's/[.]0//g' %s > %s" % (newPath + "temp", newPath))
		# os.remove(newPath + "temp")
		return(cnaDf)

	def process_steps(self, filePath, **kwargs):
		logger.info('PROCESSING %s' % filePath)
		databaseToSynIdMappingDf = kwargs['databaseToSynIdMappingDf']
		databaseSynId = kwargs['databaseSynId']
		newPath = kwargs['newPath']
		test = kwargs['test']

		cnaDf = pd.read_csv(filePath, sep="\t",comment="#")
		newCNA = self._process(cnaDf, test=test)

		centerMafSynId = databaseToSynIdMappingDf.Id[databaseToSynIdMappingDf['Database'] == "centerMaf"][0]
		if not newCNA.empty:
		# 	cols = newCNA.columns   
		# 	process_functions.updateData(self.syn, databaseSynId, newCNA, self.center, cols, toDelete=True)
		# 	newCNA.to_csv(newPath, sep="\t",index=False)
			newCNA.to_csv(newPath, sep="\t",index=False)
			self.syn.store(synapseclient.File(newPath, parent=centerMafSynId))
		return(newPath)

	def _validate(self, cnvDF, noSymbolCheck, test=False):
		total_error = ""
		warning = ""
		cnvDF.columns = [col.upper() for col in cnvDF.columns]

		if cnvDF.columns[0] != "HUGO_SYMBOL":
			total_error += "Your cnv file's first column must be Hugo_Symbol\n"
		haveColumn = process_functions.checkColExist(cnvDF, "HUGO_SYMBOL")
		if haveColumn:
			keepSymbols = cnvDF["HUGO_SYMBOL"]
			cnvDF.drop("HUGO_SYMBOL", axis=1, inplace=True)

		if sum(cnvDF.apply(lambda x: sum(x.isnull()))) > 0:
			total_error += "Your cnv file must not have any empty values\n"

		if process_functions.checkColExist(cnvDF, "ENTREZ_GENE_ID"):
			del cnvDF['ENTREZ_GENE_ID']

		if not all(cnvDF.applymap(lambda x: isinstance(x, (int, float))).all()):
			total_error += "All values must be numerical values\n"
		else:
			cnvDF['HUGO_SYMBOL'] = keepSymbols
			if haveColumn and not noSymbolCheck:
				#logger.info("VALIDATING %s GENE SYMBOLS" % os.path.basename(filePath))

				bedSynId = process_functions.getDatabaseSynId(self.syn, "bed", test=test)
				bed = self.syn.tableQuery("select Hugo_Symbol, ID from %s where CENTER = '%s'" % (bedSynId, self.center))
				bedDf = bed.asDataFrame()
				cnvDF['remapped'] = cnvDF['HUGO_SYMBOL'].apply(lambda x: validateSymbol(x, bedDf))
				cnvDF = cnvDF[~cnvDF['remapped'].isnull()]

				#Do not allow any duplicated genes after symbols have been remapped
				if sum(cnvDF['remapped'].duplicated()) >0:
					total_error+= "Your CNA file has duplicated Hugo_Symbols (After remapping of genes): %s -> %s.\n" % (",".join(cnvDF['HUGO_SYMBOL'][cnvDF['remapped'].duplicated(keep=False)]), ",".join(cnvDF['remapped'][cnvDF['remapped'].duplicated(keep=False)]))
		return(total_error, warning)

	# VALIDATION
	def validate_steps(self, filePathList, **kwargs):
		"""
		This function validates the CNV (linear or discrete) file to make sure it adhere to the genomic SOP.
		
		:params filePath:     Path to CNV file

		:returns:             Text with all the errors in the CNV file
		"""
		filePath = filePathList[0]
		logger.info("VALIDATING %s" % os.path.basename(filePath))
		test = kwargs['testing']
		noSymbolCheck = kwargs['noSymbolCheck']
		cnvDF = pd.read_csv(filePath,sep="\t",comment="#")
		return(self._validate(cnvDF, noSymbolCheck, test))
