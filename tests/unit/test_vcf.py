import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises
import os
import sys
from genie import vcf

def test_processing():
	pass


def test_validation():


	syn = mock.create_autospec(synapseclient.Synapse) 

	vcfClass = vcf(syn, "SAGE")

	assert_raises(AssertionError, vcfClass.validateFilename, ["foo"])
	assert_raises(AssertionError, vcfClass.validateFilename, ["GENIE-SAGE-ID1.bed"])
	assert vcfClass.validateFilename(["GENIE-SAGE-ID1.vcf"]) == "vcf"

	vcfDf = pd.DataFrame({"#CHROM" :['2', '9', '12'],
				  	      "POS" :[69688533, 99401860, 53701241],
				  	      "ID" :['AAK1','AAED1','AAAS'],
				  	      "REF" :['AAK1','AAED1','AAAS'],
				  	      "ALT" :['AAK1','AAED1','AAAS'],
				  	      "QUAL":['AAK1','AAED1','AAAS'],
				  	      "FILTER":['AAK1','AAED1','AAAS'],
				  	      "INFO":['AAK1','AAED1','AAAS']})

	error, warning = vcfClass._validate(vcfDf)
	assert error == ""
	assert warning == ""	

	vcfDf = pd.DataFrame({"POS" :[69688533, 99401860, 53701241],
				  	      "ID" :['AAK1','AAED1','AAAS'],
				  	      "REF" :['AAK1','AAED1','AAAS'],
				  	      "ALT" :['AAK1','AAED1','AAAS'],
				  	      "QUAL":['AAK1','AAED1','AAAS'],
				  	      "FILTER":['AAK1','AAED1','AAAS'],
				  	      "INFO":['AAK1','AAED1','AAAS'],
				  	      "FOO":['AAK1','AAED1','AAAS'],
				  	      "DOO":['AAK1','AA ED1','AAAS']})

	error, warning = vcfClass._validate(vcfDf)
	expectedError = ("Your vcf file must have these headers: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.\n"
					 "Your vcf file must have FORMAT header if genotype columns exist.\n")
	expectedWarning = ("Your vcf file should not have any white spaces in any of the columns.\n")
	assert error == expectedError
	assert warning == expectedWarning	


	vcfDf = pd.DataFrame({"#CHROM" :['chr2', 'chrM', '12'],
			  	      "POS" :[69688533, 99401860, 53701241],
			  	      "ID" :['AAK1','AAED1','AAAS'],
			  	      "REF" :['AAK1','AAED1','AAAS'],
			  	      "ALT" :['AAK1','AAED1','AAAS'],
			  	      "QUAL":['AAK1','AAED1','AAAS'],
			  	      "FILTER":['AAK1','AAED1','AAAS'],
			  	      "INFO":['AAK1','AAED1','AAAS']})

	error, warning = vcfClass._validate(vcfDf)
	expectedError = ("Your vcf file must not have variants on chrM.\n")
	expectedWarning = ("Your vcf file should not have the chr prefix in front of chromosomes.\n")
	assert error == expectedError
	assert warning == expectedWarning