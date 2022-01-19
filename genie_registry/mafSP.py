import os
import logging

import pandas as pd
from genie import process_functions

from .maf import maf

logger = logging.getLogger(__name__)


class mafSP(maf):
    """
    MAF SP file format validation / processing
    """

    _fileType = "mafSP"

    def _validateFilename(self, filePath):
        """
        Validates filename.  Should be
        nonGENIE_data_mutations_extended_CENTER.txt
        """
        assert os.path.basename(
            filePath[0]
        ) == "nonGENIE_data_mutations_extended_{}.txt".format(self.center)

    def formatMAF(self, mafDf):
        """
        The sponsored project maf file doesn't have less columns
        """
        mafDf["Center"] = self.center
        mafDf["Tumor_Sample_Barcode"] = [
            process_functions.checkGenieId(i, self.center)
            for i in mafDf["Tumor_Sample_Barcode"]
        ]
        mafDf["Sequence_Source"] = float("nan")
        mafDf["Sequencer"] = float("nan")
        mafDf["Validation_Status"][
            mafDf["Validation_Status"].isin(["Unknown", "unknown"])
        ] = ""
        return mafDf

    def storeProcessedMaf(self, filePath, mafSynId, centerMafSynId, isNarrow=False):
        """
        Stores SP maf to database
        """
        logger.info("STORING %s" % filePath)
        mafDataFrame = pd.read_csv(filePath, sep="\t")
        process_functions.updateData(
            self.syn,
            mafSynId,
            mafDataFrame,
            self.center,
            filterByColumn="Center",
            toDelete=True,
        )
        return filePath
