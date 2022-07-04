# import synapseclient
# import pandas as pd
# import mock
# from nose.tools import assert_raises
# import os
# import sys

# from genie.cbs import cbs

# def test_processing():

#   syn = mock.create_autospec(synapseclient.Synapse)

#   cbsClass = cbs(cbs, "SAGE")
#   pass

# def test_validation():

#   syn = mock.create_autospec(synapseclient.Synapse)

#   cbsClass = cbs(cbs, "SAGE")

#   assert_raises(AssertionError, cbsClass.validateFilename, ["foo"])
#   assert_raises(AssertionError, cbsClass.validateFilename, ["GENIE-test.txt"])
#   assert cbsClass.validateFilename(["GENIE-SAGE-test.cbs"]) == "cbs"

#   cbsDf = pd.DataFrame({"ID":['ID1','ID2','ID3','ID4','ID5'],
#                         "Chr":[1,2,3,4,5],
#                         "Start":[1,2,3,4,3],
#                         "End":[1,2,3,4,3],
#                         "Probes":[1,2,3,4,3],
#                         "Log2Ratio":[1,2,3,4,3]})

#   error, warning = cbsClass._validate(cbsDf)
#   assert error == ""
#   assert warning == ""

#   cbsDf = pd.DataFrame({"ID":['ID1','ID2','ID3','ID4','ID5'],
#                         "Chr":[1,2,float('nan'),4,5],
#                         "Start":[1,2,3,4,float('nan')],
#                         "End":[1,2,3,4,3],
#                         "Log2Ratio":[1,2,3,4,3]})
#   expectedErrors = ("Your cbs file must at least have these headers: Probes.\n"
#                     "cbs: No null or empty values allowed in column(s): Chr, Start.\n")
#   error, warning = cbsClass._validate(cbsDf)
#   assert error == expectedErrors
#   assert warning == ""

#   cbsDf = pd.DataFrame({"ID":['ID1','ID2','ID3','ID4','ID5'],
#                         "Chr":[1,2,3,4,5],
#                         "Start":[1,2,3,4.3,3],
#                         "End":[1,2,3.4,4,3],
#                         "Probes":[1,2,3,33.3,3],
#                         "Log2Ratio":[1,2,'f.d',4,3]})
#   error, warning = cbsClass._validate(cbsDf)
#   expectedErrors = ("cbs: Only integars allowed in these column(s): Start, End, Probes.\n"
#                     "cbs: Only numerical values allowed in Log2Ratio.\n")
#   assert error == expectedErrors
#   assert warning == ""
