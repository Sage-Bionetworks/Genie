import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR,"../../processing"))

from mutationsInCis import mutationsInCis

def test_processing():
	pass

def test_validation():
	def createMockTable(dataframe):
		table = mock.create_autospec(synapseclient.table.CsvFileTable)
		table.asDataFrame.return_value= dataframe
		return(table)

	def table_query_results(*args):
		return(table_query_results_map[args])

	databaseMapping = pd.DataFrame(dict(Database=['mutationsInCis'],
					  	        		Id=['syn7765462']))
	mutCisDf = pd.DataFrame(dict(Flag=['mutationsInCis'],
				  	        	   Center=['SAGE'],
				  	        	   Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
				  	        	   Variant_Classification=["Nonsense_Mutation"],
				  	        	   Hugo_Symbol=["AKT1"],
				  	        	   HGVSp_Short=["p.1234"],
				  	        	   Chromosome=["1"],
				  	        	   Start_Position=[324234],
				  	        	   Reference_Allele=["AGC"],
				  	        	   Tumor_Seq_Allele2=["GCCCT"],
				  	        	   t_alt_count_num=[3],
				  	        	   t_depth=[234]))

	#This is the gene positions that all bed dataframe will be processed against
	table_query_results_map = {
	("SELECT * FROM syn10967259",) : createMockTable(databaseMapping),
	("select * from syn7765462 where Center = 'SAGE'",) : createMockTable(mutCisDf)
	}	

	syn = mock.create_autospec(synapseclient.Synapse) 
	syn.tableQuery.side_effect=table_query_results

	mutCis = mutationsInCis(syn, "SAGE")

	assert_raises(AssertionError, mutCis.validateFilename, ["foo"])
	assert mutCis.validateFilename(["mutationsInCis_filtered_samples.csv"]) == "mutationsInCis"


	mutCisDf = pd.DataFrame(dict(Flag=['mutationsInCis'],
				  	        	   Center=['SAGE'],
				  	        	   Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
				  	        	   Variant_Classification=["Nonsense_Mutation"],
				  	        	   Hugo_Symbol=["AKT1"],
				  	        	   HGVSp_Short=["p.1234"],
				  	        	   Chromosome=["1"],
				  	        	   Start_Position=[324234],
				  	        	   Reference_Allele=["AGC"],
				  	        	   Tumor_Seq_Allele2=["GCCCT"],
				  	        	   t_alt_count_num=[3],
				  	        	   t_depth=[234]))

	error, warning = mutCis._validate(mutCisDf)
	assert error == ""
	assert warning == ""
	mutCisDf = pd.DataFrame(dict(Flag=['mutationsInCis'],
			  	        	   Center=['SAGE'],
			  	        	   Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
			  	        	   Hugo_Symbol=["AKT1"],
			  	        	   HGVSp_Short=["p.1234"],
			  	        	   Chromosome=["1"],
			  	        	   Reference_Allele=["AGC"],
			  	        	   Tumor_Seq_Allele2=["GCCCT"],
			  	        	   t_alt_count_num=[3],
			  	        	   t_depth=[234]))

	error, warning = mutCis._validate(mutCisDf)
	expectedErrors = ("Mutations In Cis Filter File: Must at least have these headers: Variant_Classification,Start_Position.\n")
	assert error == expectedErrors

	mutCisDf = pd.DataFrame(dict(Flag=['mutationsInCis'],
			  	        	   Center=['SAGE'],
			  	        	   Tumor_Sample_Barcode=["GENIE-SAGE-ID1-1"],
			  	        	   Variant_Classification=["Nonsense_Mutation"],
			  	        	   Hugo_Symbol=["AKT1"],
			  	        	   HGVSp_Short=["foo"],
			  	        	   Chromosome=["1"],
			  	        	   Start_Position=[324234],
			  	        	   Reference_Allele=["AGC"],
			  	        	   Tumor_Seq_Allele2=["GCCCT"],
			  	        	   t_alt_count_num=[3],
			  	        	   t_depth=[234]))
	
	error, warning = mutCis._validate(mutCisDf)
	expectedErrors = ("Mutations In Cis Filter File: All variants must come from the original mutationInCis_filtered_samples.csv file in each institution's staging folder.\n")
	assert error == expectedErrors
