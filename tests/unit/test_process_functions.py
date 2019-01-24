import synapseclient
import pandas as pd
import mock
from nose.tools import assert_raises, assert_equals
import genie

class test_update_database:
	
	def setup(self):
		self.databasedf = pd.DataFrame({'UNIQUE_KEY':['test1','test2','test3'],
							  "test":['test1','test2','test3'],
							  "foo":[1,2,3],
							  "baz":[float('nan'),float('nan'),float('nan')]})
		self.databasedf.index = ['1_3','2_3','3_5']
	
	def test_append_rows_to_database(self):

		new_datadf = pd.DataFrame({'UNIQUE_KEY':['test1','test2','test3','test4'],
							  "test":['test1','test2','test3','test4'],
							  "foo":[1,2,3,4],
							  "baz":[float('nan'),float('nan'),float('nan'),3.2]})
		expecteddf = pd.DataFrame({'test':['test4'],
										  'foo':[4],
										  'baz':[3.2]})
		append_rows = genie.process_functions._append_rows(new_datadf, self.databasedf, 'UNIQUE_KEY')
		append_rows.fillna('',inplace=True)
		expecteddf.fillna('',inplace=True)
		assert append_rows.equals(expecteddf[append_rows.columns])

	def test_append_no_rows(self):
		append_rows = genie.process_functions._append_rows(self.databasedf, self.databasedf, 'UNIQUE_KEY')
		assert append_rows.empty

	def test_update_rows_to_database(self):
		new_datadf = pd.DataFrame({'UNIQUE_KEY':['test1','test2','test3'],
							  "test":['test','test2','test3'],
							  "foo":[1,3,3],
							  "baz":[float('nan'),5,float('nan')]})

		expecteddf = pd.DataFrame({"test":['test','test2'],
							  "foo":[1,3],
							  "baz":['',5],
							  'ROW_ID':['1','2'],
							  'ROW_VERSION':['3','3']})
		update_rows = genie.process_functions._update_rows(new_datadf, self.databasedf, 'UNIQUE_KEY')
		assert update_rows.equals(expecteddf[update_rows.columns])

	def test_update_rows_to_database_maintaintype(self):
		new_datadf = pd.DataFrame({'UNIQUE_KEY':['test1','test2','test3'],
							  "test":['test1','test2','test3'],
							  "foo":[1,3,3],
							  "baz":[float('nan'),5,float('nan')]})
		#Test that the datatype passed into from new_datadf gets preserved
		expecteddf = pd.DataFrame({"test":['test2'],
							  "foo":[3],
							  "baz":[5],
							  'ROW_ID':['2'],
							  'ROW_VERSION':['3']})
		expecteddf = expecteddf.astype({'baz':object}) 
		update_rows = genie.process_functions._update_rows(new_datadf, self.databasedf, 'UNIQUE_KEY')
		assert update_rows.equals(expecteddf[update_rows.columns])

	def test_update_empty_database(self):
		new_datadf = pd.DataFrame({'UNIQUE_KEY':['test4'],
									  "test":['test'],
									  "foo":[1],
									  "baz":[float('nan')]})
		update_rows = genie.process_functions._update_rows(new_datadf, self.databasedf, 'UNIQUE_KEY')
		assert update_rows.empty

	def test_update_no_rows(self):
		update_rows = genie.process_functions._update_rows(self.databasedf, self.databasedf, 'UNIQUE_KEY')
		assert update_rows.empty

	def test_delete_rows_to_database(self):
		new_datadf = pd.DataFrame({'UNIQUE_KEY':['test1'],
							  "test":['test1'],
							  "foo":[1],
							  "baz":[float('nan')]})
		expecteddf = pd.DataFrame({0:['2','3'],
							  	   1:['3','5']})
		delete_rows = genie.process_functions._delete_rows(new_datadf, self.databasedf, 'UNIQUE_KEY')
		assert delete_rows.equals(expecteddf)

	def test_delete_no_rows(self):
		delete_rows = genie.process_functions._delete_rows(self.databasedf, self.databasedf, 'UNIQUE_KEY')
		assert delete_rows.empty

