import mock
import synapseclient
import pandas as pd

from genie import veoibd_clinical
import pytest
import datetime


def createMockTable(dataframe):
    table = mock.create_autospec(synapseclient.table.CsvFileTable)
    table.asDataFrame.return_value = dataframe
    return(table)


def table_query_results(*args):
    return(table_query_results_map[args])


with_nan = pd.DataFrame(dict(
    CODE=[1, 2, 3, 4, float('nan')],
    CBIO_LABEL=['Test', 'Why', 'foo', 'Me', 'Unknown'],
    DESCRIPTION=['non', 'asdf', 'asdf', 'asdff', 'asdfasdf']))

no_nan = pd.DataFrame(dict(
    CODE=[1, 2, 3, 4],
    CBIO_LABEL=['Test', 'Why', 'foo', 'Me'],
    DESCRIPTION=['non', 'asdf', 'asdf', 'asdff']))

sexdf = pd.DataFrame(dict(
    CODE=[1, 2, float('nan')],
    CBIO_LABEL=['Male', 'Female', 'Unknown'],
    DESCRIPTION=['Male', 'Female', 'Not coded']))

table_query_results_map = {
    ("SELECT * FROM syn7434222",): createMockTable(sexdf),
    ("SELECT * FROM syn7434236",): createMockTable(with_nan),
    ("SELECT * FROM syn7434242",): createMockTable(with_nan),
    ("SELECT * FROM syn7434273",): createMockTable(no_nan)}

syn = mock.create_autospec(synapseclient.Synapse)
syn.tableQuery.side_effect = table_query_results
clin_class = veoibd_clinical.ClinicalIndividual(syn, "SAGE")


def test_filetype():
    assert clin_class._fileType == "veoibd_clinical_individual"


@pytest.fixture(params=[
    (["clinical_individual.tsv"]),
    (["foo", "clinical_individual.csv"])
    ])
def filename_fileformat_map(request):
    return request.param


def test_incorrect_validatefilename(filename_fileformat_map):
    """Each parameter passed from filename_fileformat_map should fail with an assertion error.
    """
    filepath_list = filename_fileformat_map
    with pytest.raises(AssertionError):
        clin_class.validateFilename(filepath_list)


def test_correct_validatefilename():
    assert clin_class.validateFilename(["clinical_individual.csv"]) == "veoibd_clinical_individual"    

# def test_patient_fillvs__process():
#     '''
#     Test filling out of vital status values
#     This will be removed once vital status values are required
#     - capitalized column headers
#     - remapping of values
#     - Fill out CENTER column
#     - Append GENIE-CENTER-..
#     '''
#     expected_patientdf = pd.DataFrame(dict(
#         PATIENT_ID=["GENIE-SAGE-ID1", "GENIE-SAGE-ID2", "GENIE-SAGE-ID3",
#                     "GENIE-SAGE-ID4", "GENIE-SAGE-ID5"],
#         SEX=['Male', 'Female', 'Male', 'Female', 'Unknown'],
#         PRIMARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         SECONDARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         TERTIARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         ETHNICITY=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
#         CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
#         INT_DOD=['Not Collected', 'Not Collected', 'Not Collected',
#                  'Not Collected', 'Not Collected'],
#         INT_CONTACT=['Not Collected', 'Not Collected', 'Not Collected',
#                      'Not Collected', 'Not Collected'],
#         DEAD=['Not Collected', 'Not Collected', 'Not Collected',
#               'Not Collected', 'Not Collected'],
#         YEAR_DEATH=['Not Collected', 'Not Collected', 'Not Collected',
#                     'Not Collected', 'Not Collected'],
#         YEAR_CONTACT=['Not Collected', 'Not Collected', 'Not Collected',
#                       'Not Collected', 'Not Collected']))

#     patientdf = pd.DataFrame(dict(
#         PATIENT_Id=["ID1", "ID2", "ID3", "ID4", "ID5"],
#         sex=[1, 2, 1, 2, float('nan')],
#         PRIMARY_RACE=[1, 2, 3, 4, float('nan')],
#         Secondary_RACE=[1, 2, 3, 4, float('nan')],
#         TERTIARY_RACE=[1, 2, 3, 4, float('nan')],
#         ETHNICITY=[1, 2, 3, 4, float('nan')],
#         BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
#         CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"]))
#     patient_cols = [
#         "PATIENT_ID", "SEX", "PRIMARY_RACE", "SECONDARY_RACE",
#         "TERTIARY_RACE", "ETHNICITY", "BIRTH_YEAR", "CENTER",
#         'YEAR_CONTACT', 'YEAR_DEATH', 'INT_CONTACT', 'INT_DOD', 'DEAD']
#     clinical_template = pd.DataFrame(columns=patient_cols)
#     new_patientdf = clin_class._process(patientdf, clinical_template)

#     assert new_patientdf.columns.isin(expected_patientdf.columns).all()
#     assert expected_patientdf.equals(new_patientdf[expected_patientdf.columns])


# def test_patient_lesscoltemplate__process():
#     '''
#     Test scope is excluding values.
#     Only those value defined by the scope will be written out
#     '''
#     expected_patientdf = pd.DataFrame(dict(
#         PATIENT_ID=["GENIE-SAGE-ID1", "GENIE-SAGE-ID2", "GENIE-SAGE-ID3",
#                     "GENIE-SAGE-ID4", "GENIE-SAGE-ID5"],
#         SEX=['Male', 'Female', 'Male', 'Female', 'Unknown'],
#         PRIMARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         SECONDARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         TERTIARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         ETHNICITY=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
#         CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
#         INT_DOD=['Not Collected', 'Not Collected', 'Not Collected',
#                  'Not Collected', 'Not Collected'],
#         INT_CONTACT=['Not Collected', 'Not Collected', 'Not Collected',
#                      'Not Collected', 'Not Collected'],
#         DEAD=['Not Collected', 'Not Collected', 'Not Collected',
#               'Not Collected', 'Not Collected'],
#         YEAR_DEATH=['Not Collected', 'Not Collected', 'Not Collected',
#                     'Not Collected', 'Not Collected'],
#         YEAR_CONTACT=['Not Collected', 'Not Collected', 'Not Collected',
#                       'Not Collected', 'Not Collected']))
#     # TEST patient processing
#     patientdf = pd.DataFrame(dict(
#         PATIENT_Id=["ID1", "ID2", "ID3", "ID4", "ID5"],
#         sex=[1, 2, 1, 2, float('nan')],
#         PRIMARY_RACE=[1, 2, 3, 4, float('nan')],
#         Secondary_RACE=[1, 2, 3, 4, float('nan')],
#         TERTIARY_RACE=[1, 2, 3, 4, float('nan')],
#         ETHNICITY=[1, 2, 3, 4, float('nan')],
#         BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
#         CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"]))
#     patient_cols = [
#         "PATIENT_ID", "SEX", "PRIMARY_RACE", "SECONDARY_RACE",
#         "TERTIARY_RACE", "ETHNICITY", "BIRTH_YEAR", "CENTER",
#         'YEAR_CONTACT', 'YEAR_DEATH']
#     clinical_template = pd.DataFrame(columns=patient_cols)
#     new_patientdf = clin_class._process(patientdf, clinical_template)

#     assert new_patientdf.columns.isin(patient_cols).all()
#     assert expected_patientdf[patient_cols].equals(new_patientdf[patient_cols])


# def test_patient_vs__process():
#     '''
#     Test vital status columns being propogated with same data
#     '''
#     expected_patientdf = pd.DataFrame(dict(
#         PATIENT_ID=["GENIE-SAGE-ID1", "GENIE-SAGE-ID2", "GENIE-SAGE-ID3",
#                     "GENIE-SAGE-ID4", "GENIE-SAGE-ID5"],
#         SEX=['Male', 'Female', 'Male', 'Female', 'Unknown'],
#         PRIMARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         SECONDARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         TERTIARY_RACE=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         ETHNICITY=['Test', 'Why', 'foo', 'Me', 'Unknown'],
#         BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
#         CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
#         YEAR_DEATH=["Unknown", "Not Collected", "Not Applicable", 1990, 1990],
#         YEAR_CONTACT=["Unknown", "Not Collected", 1990, 1990, 1990],
#         INT_CONTACT=['>32485', '<6570', 'Unknown', 'Not Collected', 2000],
#         INT_DOD=['>32485', '<6570', 'Unknown', 'Not Collected',
#                  'Not Applicable'],
#         DEAD=[True, False, 'Unknown', 'Not Collected', True]))
#     # TEST patient processing
#     # Clinical file headers are capitalized prior to processing
#     patientdf = pd.DataFrame(dict(
#         PATIENT_Id=["ID1", "ID2", "ID3", "ID4", "ID5"],
#         sex=[1, 2, 1, 2, float('nan')],
#         PRIMARY_RACE=[1, 2, 3, 4, float('nan')],
#         Secondary_RACE=[1, 2, 3, 4, float('nan')],
#         TERTIARY_RACE=[1, 2, 3, 4, float('nan')],
#         ETHNICITY=[1, 2, 3, 4, float('nan')],
#         BIRTH_YEAR=[1990, 1990, 1990, 1990, 1990],
#         CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
#         YEAR_DEATH=["Unknown", "Not Collected", "Not Applicable",
#                     1990, 1990],
#         YEAR_CONTACT=["Unknown", "Not Collected", 1990, 1990, 1990],
#         INT_CONTACT=['>32485', '<6570', 'Unknown', 'Not Collected', 2000],
#         INT_DOD=['>32485', '<6570', 'Unknown', 'Not Collected',
#                  'Not Applicable'],
#         DEAD=[True, False, 'Unknown', 'Not Collected', True]))
#     patient_cols = [
#         "PATIENT_ID", "SEX", "PRIMARY_RACE", "SECONDARY_RACE",
#         "TERTIARY_RACE", "ETHNICITY", "BIRTH_YEAR", "CENTER",
#         'YEAR_CONTACT', 'YEAR_DEATH', 'INT_CONTACT', 'INT_DOD', 'DEAD']
#     clinical_template = pd.DataFrame(columns=patient_cols)
#     new_patientdf = clin_class._process(patientdf, clinical_template)
#     assert expected_patientdf.equals(new_patientdf[expected_patientdf.columns])


# def test_sample__process():
#     '''
#     Test sample processing
#     - column headers are capitalized
#     - SEQ_DATE is normalized (Mon-YYYY, Release)
#     - Allow UNKNOWN oncotree value
#     - Add on GENIE-CENTER-...
#     - Remapping of SAMPLE_TYPE/SAMPLE_TYPE_DETAILED value
#     - SEQ_YEAR from SEQ_DATE, nan if SEQ_DATE is Release

#     '''
#     expected_sampledf = pd.DataFrame(dict(
#         SAMPLE_ID=["GENIE-SAGE-ID1-1", "GENIE-SAGE-ID2-1", "GENIE-SAGE-ID3-1",
#                    "GENIE-SAGE-ID4-1", "GENIE-SAGE-ID5-1"],
#         PATIENT_ID=["GENIE-SAGE-ID1", "GENIE-SAGE-ID2", "GENIE-SAGE-ID3",
#                     "GENIE-SAGE-ID4", "GENIE-SAGE-ID5"],
#         AGE_AT_SEQ_REPORT=[100000, 100000, 100000, 100000, 100000],
#         ONCOTREE_CODE=['AMPCA', 'UNKNOWN', 'AMPCA', 'AMPCA', 'AMPCA'],
#         SAMPLE_TYPE=['Test', 'Why', 'foo', 'Me', 'Me'],
#         CENTER=["SAGE", "SAGE", "SAGE", "SAGE", "SAGE"],
#         SAMPLE_TYPE_DETAILED=['non', 'asdf', 'asdf', 'asdff', 'asdff'],
#         SEQ_ASSAY_ID=['SAGE-1', 'SAGE-1', 'SAGE-1', 'SAGE-1', 'SAGE-1'],
#         SEQ_DATE=['Jan-2012', 'Apr-2013', 'Jul-2014', 'Oct-2015', 'Release'],
#         SEQ_YEAR=[2012, 2013, 2014, 2015, pd.np.nan]))
#     sample_cols = [
#         "SAMPLE_ID", "PATIENT_ID", "AGE_AT_SEQ_REPORT", "ONCOTREE_CODE",
#         "SAMPLE_TYPE", "SEQ_ASSAY_ID", 'SEQ_DATE', 'SAMPLE_TYPE_DETAILED',
#         'SEQ_YEAR']

#     clinical_template = pd.DataFrame(columns=sample_cols)
#     # patient = False
#     sampledf = pd.DataFrame(dict(
#         SAMPLE_ID=["ID1-1", "ID2-1", "ID3-1", "ID4-1", "ID5-1"],
#         PATIENT_ID=["ID1", "ID2", "ID3", "ID4", "ID5"],
#         Age_AT_SEQ_REPORT=[100000, 100000, 100000, 100000, 100000],
#         ONCOTree_CODE=['AMPCA', ' UNKNOWN', 'AMPCA', 'AMPCA', 'AMPCA'],
#         SAMPLE_TYPE=[1, 2, 3, 4, 4],
#         SEQ_ASSAY_ID=['SAGE-1', 'SAGE-1', 'SAGE-1', 'SAGE-1', 'SAGE-1'],
#         SEQ_DATE=['Jan-2012', 'Apr-2013', 'JUL-2014', 'Oct-2015', 'release']))

#     new_sampledf = clin_class._process(sampledf, clinical_template)
#     assert new_sampledf.columns.isin(expected_sampledf.columns).all()
#     assert expected_sampledf.equals(new_sampledf[expected_sampledf.columns])


def test_perfect__validate_individual():
    '''
    Test perfect validation
    '''
    individualdf = pd.DataFrame(dict(individual_id=["ID1", "ID2", "ID3", "ID4", "ID5"],
                                     sex=['male', 'female', 'male', 'female', 'not collected'],
                                     age=[10, 12, 30, 1, 3],
                                     birth_country=['US', 'AU', 'ES', 'US', 'AR'],
                                     ethnicity=[1, 2, 3, 4, float('nan')],
                                     center=["FOO", "FOO", "FOO", "FOO", "FOO"],
                                     family_hx_ibd=[True, False, True, True, False],
                                     degree_one_with_ibd=[True, False, False, True, False],
                                     degree_two_with_ibd=[True, False, False, False, False],
                                     initial_dx=['foo', 'foo', 'foo', 'foo', 'foo'],
                                     gi_site=['', '', '', '', ''],
                                     eim=['', '', '', '', ''],
                                     dx_perianal=['', '', '', '', ''],
                                     dx_medication=['', '', '', '', ''],
                                     comments=['', '', '', '', '']))
                                    
    error, warning = clin_class._validate(individualdf)
    assert error == ""
    assert warning == ""

def test_duplicated__validate():
    '''
    Test for duplicated individual id.
    '''
    individualdf = pd.DataFrame(dict(individual_id=["ID1", "ID1", "ID3", "ID4", "ID5"],
                                     sex=['male', 'female', 'male', 'female', 'not collected'],
                                     age=[10, 12, 30, 1, 3],
                                     birth_country=['US', 'AU', 'ES', 'US', 'AR'],
                                     ethnicity=[1, 2, 3, 4, float('nan')],
                                     center=["FOO", "FOO", "FOO", "FOO", "FOO"],
                                     family_hx_ibd=[True, False, True, True, False],
                                     degree_one_with_ibd=[True, False, False, True, False],
                                     degree_two_with_ibd=[True, False, False, False, False],
                                     initial_dx=['foo', 'foo', 'foo', 'foo', 'foo'],
                                     gi_site=['', '', '', '', ''],
                                     eim=['', '', '', '', ''],
                                     dx_perianal=['', '', '', '', ''],
                                     dx_medication=['', '', '', '', ''],
                                     comments=['', '', '', '', '']))

    error, warning = clin_class._validate(individualdf)
    print(error)
    expectedErrors = ("Found duplicated ['individual_id']'s in the file.")

    assert error == expectedErrors
    assert warning == ""

# def test_missingcols__validate():
#     '''
#     Test for missing column errors
#     '''
#     clinicaldf = pd.DataFrame()
#     with mock.patch(
#             "genie.process_functions.get_oncotree_code_mappings",
#             return_value=onco_map_dict) as mock_get_onco_map:
#         error, warning = clin_class._validate(clinicaldf, json_oncotreeurl)
#         mock_get_onco_map.called_once_with(json_oncotreeurl)
#         expected_errors = (
#             "Sample: clinical file must have SAMPLE_ID column.\n"
#             "Patient: clinical file must have PATIENT_ID column.\n"
#             "Sample: clinical file must have AGE_AT_SEQ_REPORT column.\n"
#             "Sample: clinical file must have ONCOTREE_CODE column.\n"
#             "Sample: clinical file must have SAMPLE_TYPE column.\n"
#             "Sample: clinical file must have SEQ_ASSAY_ID column.\n"
#             "Sample: clinical file must SEQ_DATE column\n"
#             "Patient: clinical file must have BIRTH_YEAR column.\n"
#             "Patient: clinical file must have SEX column.\n")

#         expected_warnings = (
#             "Patient: Must have YEAR_DEATH column for 7...release uploads.\n"
#             "Patient: Must have YEAR_CONTACT column for 7...release uploads.\n"
#             "Patient: Must have INT_CONTACT column for 7...release uploads.\n"
#             "Patient: Must have INT_DOD column for 7...release uploads.\n"
#             "Patient: Must have DEAD column for 7...release uploads.\n"
#             "Patient: clinical file doesn't have PRIMARY_RACE column. "
#             "A blank column will be added\n"
#             "Patient: clinical file doesn't have SECONDARY_RACE column. "
#             "A blank column will be added\n"
#             "Patient: clinical file doesn't have TERTIARY_RACE column. "
#             "A blank column will be added\n"
#             "Patient: clinical file doesn't have ETHNICITY column. "
#             "A blank column will be added\n")
#         assert error == expected_errors
#         assert warning == expected_warnings


# def test_errors__validate():
#     '''
#     Test for validation errors
#     '''
#     sampleDf = pd.DataFrame(dict(
#         SAMPLE_ID=[float('nan'), "ID2-1", "ID3-1", "ID4-1", "ID5-1"],
#         PATIENT_ID=["ID6", "ID2", "ID3", float('nan'), "ID5"],
#         AGE_AT_SEQ_REPORT=[10, 100000, "doo", 100000, 100000],
#         ONCOTREE_CODE=['AMPCAD', 'TESTIS', 'AMPCA', 'AMPCA', 'UCEC'],
#         SAMPLE_TYPE=[1, 2, 3, 4, float('nan')],
#         SEQ_ASSAY_ID=[float('nan'), 'Sage-1', 'SAGE-1', 'S-SAGE-1', 'SAGE-1'],
#         SEQ_DATE=['Jane-2013', 'Jan-2013', 'Jan-2013', 'Jan-2013', 'Jan-2013'],
#         YEAR_DEATH=["Unknown", "Not Collected", "Not Applicable", 19930, 1990],
#         YEAR_CONTACT=["Unknown", "Not Collected", 1990, 1990, 19940],
#         INT_CONTACT=['>32485', '<6570', 'Unknown', 'Not Collected', "foobar"],
#         INT_DOD=['>32485', '<6570', 'Unknown', 'Not Collected', 'dense'],
#         DEAD=[1, False, 'Unknown', 'Not Collected', 'Not Applicable']))

#     patientDf = pd.DataFrame(dict(
#         PATIENT_ID=["ID6", "ID2", "ID3", float("nan"), "ID5"],
#         SEX=[1, 2, 1, 5, float('nan')],
#         PRIMARY_RACE=[1, 2, 3, 6, float('nan')],
#         SECONDARY_RACE=[1, 2, 3, 6, float('nan')],
#         TERTIARY_RACE=[1, 2, 3, 6, float('nan')],
#         ETHNICITY=[1, 2, 3, 6, float('nan')],
#         BIRTH_YEAR=[1990, 1990, datetime.datetime.utcnow().year + 1,
#                     1990, 1990],
#         CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"]))
#     clinicalDf = patientDf.merge(sampleDf, on="PATIENT_ID")
#     with mock.patch(
#             "genie.process_functions.get_oncotree_code_mappings",
#             return_value=onco_map_dict) as mock_get_onco_map:
#         error, warning = clin_class._validate(clinicalDf, json_oncotreeurl)
#         mock_get_onco_map.called_once_with(json_oncotreeurl)

#         expectedErrors = (
#             "Sample: PATIENT_ID's much be contained in the SAMPLE_ID's "
#             "(ex. SAGE-1 <-> SAGE-1-2)\n"
#             "Patient: All samples must have associated patient information "
#             "and no null patient ids allowed. "
#             "These samples are missing patient data: ID4-1\n"
#             "Sample: Please double check your AGE_AT_SEQ_REPORT.  "
#             "It must be an integer or 'Unknown'.\n"
#             "Sample: Please double check that all your ONCOTREE CODES exist "
#             "in the mapping. You have 1 samples that don't map. "
#             "These are the codes that don't map: AMPCAD\n"
#             "Sample: Please double check your SAMPLE_TYPE column. "
#             "No null values allowed.\n"
#             "Sample: Please double check your SEQ_ASSAY_ID columns, "
#             "there are empty rows.\n"
#             "Sample: Please make sure your SEQ_ASSAY_IDs start with your "
#             "center abbreviation: S-SAGE-1.\n"
#             "Sample: SEQ_DATE must be one of five values- "
#             "For Jan-March: use Jan-YEAR. "
#             "For Apr-June: use Apr-YEAR. "
#             "For July-Sep: use Jul-YEAR. "
#             "For Oct-Dec: use Oct-YEAR. (ie. Apr-2017) "
#             "For values that don't have SEQ_DATES that you want "
#             "released use 'release'.\n"
#             "Patient: Please double check your BIRTH_YEAR column, "
#             "it must be an integer in YYYY format > {year} or 'Unknown'.  "
#             "Support for blank values will be deprecated in 7...releases.\n"
#             "Patient: Please double check your YEAR_DEATH column, it must be "
#             "an integer in YYYY format, 'Unknown', 'Not Applicable' "
#             "or 'Not Collected'.\n"
#             "Patient: Please double check your YEAR_CONTACT column, it must "
#             "be an integer in YYYY format, 'Unknown' or 'Not Collected'.\n"
#             "Patient: Please double check your INT_CONTACT column, "
#             "it must be an integer, '>32485', '<6570', 'Unknown' or "
#             "'Not Collected'.\n"
#             "Patient: Please double check your INT_DOD column, it must be "
#             "an integer, '>32485', '<6570', 'Unknown', 'Not Collected' or "
#             "'Not Applicable'.\n"
#             "Patient: Please double check your DEAD column, it must be "
#             "True, False, 'Unknown' or 'Not Collected'.\n"
#             "Patient: Please double check your PRIMARY_RACE column.  "
#             "This column must be these values 1, 2, 3, 4, or blank.\n"
#             "Patient: Please double check your SECONDARY_RACE column.  "
#             "This column must be these values 1, 2, 3, 4, or blank.\n"
#             "Patient: Please double check your TERTIARY_RACE column.  "
#             "This column must be these values 1, 2, 3, 4, or blank.\n"
#             "Patient: Please double check your SEX column.  "
#             "This column must be these values 1, 2, or blank.\n"
#             "Patient: Please double check your ETHNICITY column.  "
#             "This column must be these values 1, 2, 3, 4, or blank.\n")
#         expectedWarnings = (
#             "Sample: All patients must have associated sample information. "
#             "These patients are missing sample data: ID6\n"
#             "Sample: Some SAMPLE_IDs have conflicting SEX and "
#             "ONCOTREE_CODES: ID2-1,ID5-1\n")
#         assert error == expectedErrors.format(
#             year=datetime.datetime.utcnow().year)
#         print(warning)
#         assert warning == expectedWarnings


# def test_duplicated__validate():
#     '''
#     Test for duplicated SAMPLE_ID error and and in the case that
#     both sample and patient
#     are uploaded, it could be a duplicated PATIENT_ID error
#     '''
#     patientDf = pd.DataFrame(dict(
#         PATIENT_ID=["ID1", "ID1", "ID3", "ID4", "ID5"],
#         SEX=[1, 2, 1, 2, float('nan')],
#         PRIMARY_RACE=[1, 2, 3, 4, float('nan')],
#         SECONDARY_RACE=[1, 2, 3, 4, float('nan')],
#         TERTIARY_RACE=[1, 2, 3, 4, float('nan')],
#         ETHNICITY=[1, 2, 3, 4, float('nan')],
#         BIRTH_YEAR=[float('nan'), 1990, 1990, 1990, 1990],
#         CENTER=["FOO", "FOO", "FOO", "FOO", "FOO"],
#         YEAR_DEATH=["Unknown", "Not Collected", "Not Applicable",
#                     1990, 1990],
#         YEAR_CONTACT=["Unknown", "Not Collected", 1990, 1990, 1990],
#         INT_CONTACT=['>32485', '<6570', 'Unknown', 'Not Collected', 2000],
#         INT_DOD=['>32485', '<6570', 'Unknown', 'Not Collected',
#                  'Not Applicable'],
#         DEAD=[True, False, 'Unknown', 'Not Collected', True]))

#     sampleDf = pd.DataFrame(dict(
#         SAMPLE_ID=["ID1-1", "ID3-1", "ID4-1", "ID5-1"],
#         PATIENT_ID=["ID1", "ID3", "ID4", "ID5"],
#         AGE_AT_SEQ_REPORT=[100000, 100000, 100000, 100000],
#         ONCOTREE_CODE=['AMPCA', 'UNKNOWN', 'AMPCA', 'AMPCA'],
#         SAMPLE_TYPE=[1, 3, 4, 4],
#         SEQ_ASSAY_ID=['SAGE-1-1', 'SAGE-1', 'SAGE-1', 'SAGE-1'],
#         SEQ_DATE=['Jan-2013', 'Jul-2013', 'Oct-2013', 'release']))

#     clinicalDf = patientDf.merge(sampleDf, on="PATIENT_ID")
#     with mock.patch(
#             "genie.process_functions.get_oncotree_code_mappings",
#             return_value=onco_map_dict) as mock_get_onco_map:
#         error, warning = clin_class._validate(clinicalDf, json_oncotreeurl)
#         mock_get_onco_map.called_once_with(json_oncotreeurl)
#         expectedErrors = (
#             "Sample: No duplicated SAMPLE_ID in the "
#             "sample file allowed.\n"
#             "If there are no duplicated SAMPLE_IDs, and both sample and "
#             "patient files are uploaded, then please check to make sure no "
#             "duplicated PATIENT_IDs exist in the patient file.\n")

#         assert error == expectedErrors
#         assert warning == ""
