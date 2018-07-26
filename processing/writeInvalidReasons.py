import synapseclient
import os
import process_functions
import argparse

def writeInvalidReasons(x, syn, errorFile):
	ent = syn.get(x['id'],downloadFile=False)
	errors = x['errors']
	showName = "\t" + ent.name + " (%s):\n\n" % ent.id
	errorFile.write(showName)
	errorFile.write(errors + "\n\n\n")

def main():
	parser = argparse.ArgumentParser(description='Write invalid reasons')
	parser.add_argument("--pemFile", type=str, help="Path to PEM file (genie.pem)")
	args = parser.parse_args()
	
	syn = process_functions.synLogin(args)
	CENTER_MAPPING = syn.tableQuery('SELECT * FROM syn10061452 where inputSynId is not null and release is true')
	CENTER_MAPPING_DF = CENTER_MAPPING.asDataFrame()

	for center in CENTER_MAPPING_DF['center']:
		print(center)
		errorTracker = syn.tableQuery("SELECT * FROM syn10153306 where center = '%s'" % center)
		errorTrackerDf = errorTracker.asDataFrame()
		stagingSynId = CENTER_MAPPING_DF['stagingSynId'][CENTER_MAPPING_DF['center'] == center][0]
		with open(center + "_errors.txt", 'w') as errorFile:
			errorTrackerDf.apply(lambda x: writeInvalidReasons(x, syn, errorFile),axis=1)
			if errorTrackerDf.empty:
				errorFile.write("No errors!")
		entFile = synapseclient.File(center + "_errors.txt",parentId = stagingSynId)
		syn.store(entFile)
		os.remove(center + "_errors.txt")

if __name__ == "__main__":
	main()