from __future__ import absolute_import
from genie import maf, process_functions
import os
import logging
import pandas as pd
logger = logging.getLogger(__name__)

class mafSP(maf):

    _fileType = "mafSP"

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "nonGENIE_data_mutations_extended_%s.txt" % self.center

    def formatMAF(self, mafDf):
        #keepMafColumns = ['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2','Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status','Mutation_Status','Sequencing_Phase','Sequence_Source','Validation_Method','Score','BAM_File','Sequencer','HGVSp_Short','t_ref_count','t_alt_count','n_ref_count','n_alt_count','Protein_position','Codons','SWISSPROT','RefSeq','t_depth','n_depth','FILTER']
        #mafDf = mafDf[keepMafColumns]
        mafDf['Center'] = self.center
        mafDf['Tumor_Sample_Barcode'] = [process_functions.checkGenieId(i,self.center) for i in mafDf['Tumor_Sample_Barcode']]
        mafDf['Sequence_Source'] = pd.np.nan
        mafDf['Sequencer'] = pd.np.nan
        mafDf['Validation_Status'][mafDf['Validation_Status'].isin(["Unknown","unknown"])] = ''
        return(mafDf)

    def storeProcessedMaf(self, filePath, mafSynId, centerMafSynId, isNarrow=False):
        logger.info('STORING %s' % filePath)
        database = self.syn.get(mafSynId)
        mafDataFrame = pd.read_csv(filePath,sep="\t")
        process_functions.updateData(self.syn, mafSynId, mafDataFrame, self.center, filterByColumn="Center", toDelete=True)
        return(filePath)