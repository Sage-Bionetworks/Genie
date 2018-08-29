import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys

from genie.bed import bed
from genie.bedSP import bedSP

def test_processing():
	def createMockTable(dataframe):
		table = mock.create_autospec(synapseclient.table.CsvFileTable)
		table.asDataFrame.return_value= dataframe
		return(table)

	def table_query_results(*args):
		return(table_query_results_map[args])

	symbols = pd.DataFrame(dict(hgnc_symbol=['AAK1','AAED1','AAAS','PINLYP','XRCC1'],
					  	        chromosome_name=['2', '9', '12','19','19'],
					  	        start_position=[69688532, 99401859, 53701240,44080952,44047192],
					  	        end_position=[69901481, 99417585, 53718648,44088116,44084625]))
	#This is the gene positions that all bed dataframe will be processed against
	table_query_results_map = {
	("SELECT * FROM syn11806563",) : createMockTable(symbols),
	}	


	syn = mock.create_autospec(synapseclient.Synapse) 
	syn.tableQuery.side_effect=table_query_results

	bedClass = bed(syn, "SAGE")
	bedSPclass = bedSP(syn,"SAGE")

	expectedBedDf = pd.DataFrame(dict(Chromosome =['2', '9', '12','19','19'],
						  	          Start_Position =[69688533, 99401860, 53701241,44084466,44084466],
						  	          End_Position =[69901480, 99417584, 53718647,44084638,44084638],
						  	          Hugo_Symbol =['AAK1','AAED1','AAAS','XRCC1','PINLYP'],
						  	          includeInPanel =[True, True, True, True,True],
						  	          ID=['AAK1','AAED1','AAAS','XRCC1','foo'],
						  	          SEQ_ASSAY_ID=['SAGE-TEST','SAGE-TEST','SAGE-TEST','SAGE-TEST','SAGE-TEST'],
						  	          Feature_Type=['exon','exon','exon','exon','exon'],
						  	          CENTER=['SAGE','SAGE','SAGE','SAGE','SAGE']))
	
	expectedBedDf.sort_values("ID",inplace=True)
	expectedBedDf.reset_index(drop=True, inplace=True)
	bedDf = pd.DataFrame({0:['2', '9', '12','19','19'],
				  	      1:[69688533, 99401860, 53701241,44084466,44084466],
				  	      2:[69901480, 99417584, 53718647,44084638,44084638],
				  	      3:['AAK1','AAED1','AAAS','XRCC1','foo'],
				  	      4:['d','d','d','d','d']})

	seq_assay_id = "SAGE-Test"
	newPath = "new.bed"
	parentId = "synTest"
	newBedDf = bedClass._process(bedDf, seq_assay_id, newPath, parentId, createPanel=False)
	newBedDf.sort_values("ID",inplace=True)
	newBedDf.reset_index(drop=True, inplace=True)
	assert expectedBedDf.equals(newBedDf[expectedBedDf.columns])

	expectedBedDf = pd.DataFrame(dict(Chromosome =['2', '9', '12','19'],
					  	          Start_Position =[69688432, 1111, 53700240,44080953],
					  	          End_Position =[69689532, 1111, 53719548,44084624],
					  	          Hugo_Symbol =['AAK1',float('nan'),'AAAS',float('nan')],
					  	          includeInPanel =[True, True, False,True],
					  	          ID=['foo','bar','baz','boo'],
					  	          SEQ_ASSAY_ID=['SAGE-TEST','SAGE-TEST','SAGE-TEST','SAGE-TEST'],
					  	          Feature_Type=['exon','intergenic','exon','exon'],
					  	          CENTER=['SAGE','SAGE','SAGE','SAGE']))

	expectedBedDf.sort_values("Chromosome",inplace=True)
	expectedBedDf.reset_index(drop=True, inplace=True)

	#symbols that can't be map should be null, includeInPanel column should be included if it exists
	bedDf = pd.DataFrame({0:['2', '9', '12', '19'],
				  	      1:[69688432, 1111, 53700240, 44080953],
				  	      2:[69689532, 1111, 53719548, 44084624],
				  	      3:['foo','bar','baz', 'boo'],
				  	      4:[True, True, False, True]})

	newBedDf = bedSPclass._process(bedDf, seq_assay_id, newPath, parentId, createPanel=False)
	newBedDf.sort_values("Chromosome",inplace=True)
	newBedDf.reset_index(drop=True, inplace=True)
	assert expectedBedDf.equals(newBedDf[expectedBedDf.columns])



def test_validation():
	def createMockTable(dataframe):
		table = mock.create_autospec(synapseclient.table.CsvFileTable)
		table.asDataFrame.return_value= dataframe
		return(table)

	def table_query_results(*args):
		return(table_query_results_map[args])

	symbols = pd.DataFrame(dict(hgnc_symbol=['AAK1','AAED1','AAAS','PINLYP','XRCC1'],
					  	        chromosome_name=['2', '9', '12','19','19'],
					  	        start_position=[69688532, 99401859, 53701240,44080952,44047192],
					  	        end_position=[69901481, 99417585, 53718648,44088116,44084625]))
	#This is the gene positions that all bed dataframe will be processed against
	table_query_results_map = {
	("SELECT * FROM syn11806563",) : createMockTable(symbols),
	}	


	syn = mock.create_autospec(synapseclient.Synapse) 
	syn.tableQuery.side_effect=table_query_results

	bedClass = bed(syn, "SAGE")
	bedSPclass = bedSP(syn,"SAGE")

	assert_raises(AssertionError, bedClass.validateFilename, ["foo"])
	assert_raises(AssertionError, bedClass.validateFilename, ["SAGE-test.txt"])
	assert bedClass.validateFilename(["SAGE-test.bed"]) == "bed"

	assert_raises(AssertionError, bedSPclass.validateFilename, ["foo"])
	assert_raises(AssertionError, bedSPclass.validateFilename, ["nonGENIE_SAGE-test.txt"])
	assert bedSPclass.validateFilename(["nonGENIE_SAGE-test.bed"]) == "bedSP"

	bedDf = pd.DataFrame(dict(a =['2', '9', '12'],
				  	          b =[69688533, 99401860, 53701241],
				  	          c =[69901480, 99417584, 53718647],
				  	          d =['AAK1','AAED1','AAAS'],
				  	          e =[True, True, True]))

	error, warning = bedClass._validate(bedDf)
	assert error == ""
	assert warning == ""	

	symbols = pd.DataFrame(dict(hgnc_symbol=['AAK1','AAED1','AAAS'],
					  	        chromosome_name=['2', '9', '12'],
					  	        start_position=[69688532, 99401859, 53701240],
					  	        end_position=[69901481, 99417585, 53718648]))

	#TEST 90% boundary
	bedDf = pd.DataFrame(dict(a =['2', '9', '12'],
				  	          b =[69688432, 99416585, 53700240],
				  	          c =[69689532, 99417685, 53719548],
				  	          d =['AAK1','AAED1','AAAS'],
				  	          e =[True, True, True]))
	error, warning = bedClass._validate(bedDf)
	assert error == ""
	assert warning == ""


	bedDf = pd.DataFrame(dict(b =[69688533, 99401860, 53701241],
				  	          c =[69901480, 99417584, 53718647],
				  	          d =['AAK1','AAED1','AAAS']))

	error, warning = bedClass._validate(bedDf)
	expectedErrors = ("Your BED file must at least have four columns in this order: Chromosome, Start_Position, End_Position, Hugo_Symbol.  Make sure there are no headers in your BED file.\n")
	assert error == expectedErrors
	assert warning == ""

	bedDf = pd.DataFrame(dict(a =['2', '9', '12'],
				  	          b =[69688533, 99401860, 53701241],
				  	          c =[69901480, 99417584, 53718647],
				  	          d =['+',float('nan'),'AAAS']))
	error, warning = bedClass._validate(bedDf)
	expectedErrors = ("You cannot submit any null symbols.\n"
					  "Fourth column must be the Hugo_Symbol column, not the strand column\n")
	assert error == expectedErrors
	assert warning == ""

	bedDf = pd.DataFrame(dict(a =['2', '9', '12'],
				  	          b =['69688533', 99401860, 53701241],
				  	          c =[69901480, '99417584', 53718647],
				  	          d =['AAK1','AAED1','AAAS']))

	error, warning = bedClass._validate(bedDf)
	expectedErrors = ("The Start_Position column must only be integers. Make sure there are no headers in your BED file.\n"
					  "The End_Position column must only be integers. Make sure there are no headers in your BED file.\n")
	assert error == expectedErrors
	assert warning == ""

	#Test 90% boundary failure boundary, with incorrect gene names
	bedDf = pd.DataFrame(dict(a =['2', '9', '12'],
				  	          b =[69901381, 4345, 11111],
				  	          c =[69911481, 99417590, 11113],
				  	          d =['foo','foo','AAAS']))

	error, warning = bedSPclass._validate(bedDf)
	expectedErrors = ("You have no correct gene symbols. Make sure your gene symbol column (4th column) is formatted like so: SYMBOL(;optionaltext).  Optional text can be semi-colon separated.\n")
	expectedWarnings = ("Any gene names that can't be remapped will be null.\n")
	assert error == expectedErrors
	assert warning == expectedWarnings

	#Test overlapping boundary with correct gene names
	bedDf = pd.DataFrame(dict(a =['2', '9'],
				  	          b =[1111, 4345],
				  	          c =[69880186, 99417590],
				  	          d =['AAK1','AAED1']))

	error, warning = bedSPclass._validate(bedDf)
	assert error == ""
	assert warning == ""

	#Test 2 gene symbols returned NULL
	bedDf = pd.DataFrame(dict(a =['19'],
				  	          b =[44080953],
				  	          c =[44084624],
				  	          d =['AAK1']))

	error, warning = bedSPclass._validate(bedDf)
	expectedErrors = ("You have no correct gene symbols. Make sure your gene symbol column (4th column) is formatted like so: SYMBOL(;optionaltext).  Optional text can be semi-colon separated.\n")
	expectedWarnings = ("Any gene names that can't be remapped will be null.\n")
	assert error == expectedErrors
	assert warning == expectedWarnings


#test_validation()
#test_processing()
