from __future__ import absolute_import
from genie import FileTypeFormat, process_functions
import logging
import os
import pandas as pd
import re
logger = logging.getLogger(__name__)


class patientCounts(FileTypeFormat):

    _fileType = "patientCounts"

    _process_kwargs = ["newPath", "oncotreeLink", "databaseSynId"]
    
    _validation_kwargs = ["oncotreeLink"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "patient_counts.txt"

    def _process(self, patientCountsDf, oncotreeLink):
        patientCountsDf['CENTER'] = self.center
        #This needs to be passed into this function
        #oncotreeMapping = pd.read_csv(oncotreeLink,sep="\t")
        # primaries = [re.sub(".+[(](.+)[)]","\\1",primary) if not pd.isnull(primary) else primary for primary in oncotreeMapping.level_1]
        # secondaries = [re.sub(".+[(](.+)[)]","\\1",secondary) if not pd.isnull(secondary) else secondary for secondary in oncotreeMapping.level_2]
        # terts = [re.sub(".+[(](.+)[)]","\\1",tert) if not pd.isnull(tert) else tert for tert in oncotreeMapping.level_3]
        # quarts = [re.sub(".+[(](.+)[)]","\\1",quart) if not pd.isnull(quart) else quart for quart in oncotreeMapping.level_4]
        # quins = [re.sub(".+[(](.+)[)]","\\1",quin) if not pd.isnull(quin) else quin for quin in oncotreeMapping.level_5]
        # patientCounts['PRIMARY_CODE'] = process_functions.getPrimary(patientCounts.ONCOTREE_CODE, pd.Series(primaries), pd.Series(secondaries), pd.Series(terts), pd.Series(quarts), pd.Series(quins))
        
        oncotreeMapping = pd.read_csv(oncotreeLink,sep="\t")
        if oncotreeMapping.empty:
            oncotree_mapping_dict = process_functions.get_oncotree_code_mappings(oncotreeLink)
            patientCountsDf['PRIMARY_CODE'] = [oncotree_mapping_dict[i.upper()]['ONCOTREE_PRIMARY_NODE'] for i in patientCountsDf.ONCOTREE_CODE]
        else:
            levels = [col for col in oncotreeMapping.columns if "level_" in col]
            oncotreeDict = {}
            for level in levels:
                oncotreeDict[level] = pd.Series([re.sub(".+[(](.+)[)]","\\1",code) if not pd.isnull(code) else '' for code in oncotreeMapping[level]])

            primary = oncotreeDict.pop('level_1')
            #Need to optimize this at some point
            patientCountsDf['PRIMARY_CODE'] = patientCountsDf.ONCOTREE_CODE.apply(lambda code: process_functions.getPrimary(code, oncotreeDict, primary))
        return(patientCountsDf)

    def process_steps(self, patientCountsDf, newPath, oncotreeLink, databaseSynId):
        patientCountsDf = self._process(patientCountsDf, oncotreeLink)
        process_functions.updateData(self.syn, databaseSynId, patientCountsDf, self.center)
        patientCountsDf.to_csv(newPath, sep="\t",index=False)
        return(newPath)

    def _validate(self, patCountsDf, oncotreeLink):
        total_error = ""
        warning = ""
        oncotree_mapping = process_functions.get_oncotree_codes(oncotreeLink)
        if oncotree_mapping.empty:
            oncotree_mapping_dict = process_functions.get_oncotree_code_mappings(oncotreeLink)
            oncotree_mapping['ONCOTREE_CODE'] = oncotree_mapping_dict.keys()
        haveColumn = process_functions.checkColExist(patCountsDf, "ONCOTREE_CODE")
        if haveColumn:
            if sum(patCountsDf['ONCOTREE_CODE'].duplicated()) > 0:
                total_error += "Patient Counts: Must not have any duplicated ONCOTREE CODES.\n"
            if not all(patCountsDf['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])):
                unmapped_oncotrees = patCountsDf['ONCOTREE_CODE'][~patCountsDf['ONCOTREE_CODE'].isin(oncotree_mapping['ONCOTREE_CODE'])]
                total_error += "Patient Counts: Please double check that all your ONCOTREE CODES exist in the mapping. You have %d codes that don't map. These are the codes that don't map: %s\n" % (len(unmapped_oncotrees),",".join(set(unmapped_oncotrees)))
        else:
            total_error += "Patient Counts: File must have ONCOTREE_CODE column.\n"
        
        haveColumn = process_functions.checkColExist(patCountsDf, "NUM_PATIENTS_PD1_PDL1")
        if haveColumn:
            if not all([isinstance(i, int) for i in patCountsDf['NUM_PATIENTS_PD1_PDL1']]):
                total_error += "Patient Counts: Must not have any null values, and must be all integers.\n"
        else:
            total_error += "Patient Counts: File must have NUM_PATIENTS_PD1_PDL1 column.\n"

        return(total_error, warning)
