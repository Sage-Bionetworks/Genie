import logging
import os

import pandas as pd

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions

logger = logging.getLogger(__name__)


class seg(FileTypeFormat):

    _fileType = "seg"

    _process_kwargs = ["newPath", "databaseSynId"]

    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "genie_data_cna_hg19_%s.%s" % (
            self.center,
            self._fileType,
        )

    def _process(self, seg):
        seg.columns = [col.upper() for col in seg.columns]
        seg = seg.drop_duplicates()
        seg = seg.rename(
            columns={
                "LOC.START": "LOCSTART",
                "LOC.END": "LOCEND",
                "SEG.MEAN": "SEGMEAN",
                "NUM.MARK": "NUMMARK",
            }
        )
        seg["CHROM"] = [str(chrom).replace("chr", "") for chrom in seg["CHROM"]]
        seg["CENTER"] = self.center
        seg["LOCSTART"] = seg["LOCSTART"].astype(int)
        seg["LOCEND"] = seg["LOCEND"].astype(int)
        seg["NUMMARK"] = seg["NUMMARK"].astype(int)
        return seg

    def process_steps(self, seg, newPath, databaseSynId):
        seg = self._process(seg)
        process_functions.updateData(
            self.syn, databaseSynId, seg, self.center, toDelete=True
        )
        seg.to_csv(newPath, sep="\t", index=False)
        return newPath

    def _validate(self, segDF):
        total_error = ""
        warning = ""
        segDF.columns = [col.upper() for col in segDF.columns]

        REQUIRED_HEADERS = pd.Series(
            ["ID", "CHROM", "LOC.START", "LOC.END", "NUM.MARK", "SEG.MEAN"]
        )

        if not all(REQUIRED_HEADERS.isin(segDF.columns)):
            total_error += "Your seg file is missing these headers: %s.\n" % ", ".join(
                REQUIRED_HEADERS[~REQUIRED_HEADERS.isin(segDF.columns)]
            )
        else:
            intCols = ["LOC.START", "LOC.END", "NUM.MARK"]
            nonInts = [col for col in intCols if segDF[col].dtype != int]
            if len(nonInts) > 0:
                total_error += (
                    "Seg: Only integars allowed in these column(s): %s.\n"
                    % ", ".join(sorted(nonInts))
                )
            if not segDF["SEG.MEAN"].dtype in [float, int]:
                total_error += "Seg: Only numerical values allowed in SEG.MEAN.\n"

        checkNA = segDF.isna().apply(sum)
        nullCols = [ind for ind in checkNA.index if checkNA[ind] > 0]
        if len(nullCols) > 0:
            total_error += (
                "Seg: No null or empty values allowed in column(s): %s.\n"
                % ", ".join(sorted(nullCols))
            )

        if process_functions.checkColExist(segDF, "ID"):
            error = process_functions.validate_genie_identifier(
                identifiers=segDF["ID"], center=self.center, filename="Seg", col="ID"
            )
            total_error += error

        return (total_error, warning)
