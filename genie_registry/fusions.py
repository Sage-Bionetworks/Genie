from functools import partial
import os
import logging

import pandas as pd

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions

logger = logging.getLogger(__name__)


def validateSymbol(x, bedDf, returnMappedDf=True):
    valid = False
    gene = x["HUGO_SYMBOL"]
    if sum(bedDf["Hugo_Symbol"] == gene) > 0:
        valid = True
    elif sum(bedDf["ID"] == gene) > 0:
        mismatch = bedDf[bedDf["ID"] == gene]
        mismatch.drop_duplicates(inplace=True)
        logger.info(
            "%s will be remapped to %s" % (gene, mismatch["Hugo_Symbol"].values[0])
        )
        x["HUGO_SYMBOL"] = mismatch["Hugo_Symbol"].values[0]
    # else:
    #    logger.warning("%s cannot be remapped and will not be released. The symbol must exist in your seq assay ids (bed files) and must be mappable to a gene." % gene)
    #    x['HUGO_SYMBOL'] = float('nan')
    # x['FUSION'] = x['FUSION'].replace("%s-" % gene,"%s-" % x['HUGO_SYMBOL'])
    # x['COMMENTS'] = str(x['COMMENTS']).replace("-%s" % gene,"-%s" % str(x['COMMENTS']))
    if returnMappedDf:
        return x
    else:
        return valid


# Remap fusion's FUSION column
def remapFusion(gene_dict, DF, col):
    nonmapped = []
    total = []
    for each in DF[col]:
        for key in gene_dict:
            value = gene_dict[key]
            if value == False:
                nonmapped.append(key)
            else:
                each = each.replace("%s-" % key, "%s-" % value)
                each = each.replace("-%s fusion" % key, "-%s fusion" % value)
        total.append(each)
    DF[col] = total
    return (DF, nonmapped)


class fusions(FileTypeFormat):

    _fileType = "fusions"

    _process_kwargs = ["newPath", "databaseSynId", "databaseToSynIdMappingDf"]

    _validation_kwargs = ["nosymbol_check", "project_id"]

    # VALIDATE FILENAME
    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "data_fusions_%s.txt" % self.center

    def _process(self, fusion, databaseToSynIdMappingDf):
        fusion.columns = [col.upper() for col in fusion.columns]
        fusion["CENTER"] = self.center
        newsamples = [
            process_functions.checkGenieId(i, self.center)
            for i in fusion["TUMOR_SAMPLE_BARCODE"]
        ]
        fusion["TUMOR_SAMPLE_BARCODE"] = newsamples

        # This is temporary, because comments column will be removed
        # if fusion.get("COMMENTS") is None:
        #    fusion['COMMENTS'] = ""
        # #Will remove comments column
        # fusion['COMMENTS'] = ""
        fusion["ENTREZ_GENE_ID"] = fusion["ENTREZ_GENE_ID"].fillna(0)
        fusion = fusion.drop_duplicates()
        fusion["ID"] = fusion["HUGO_SYMBOL"].copy()
        bedSynId = process_functions.getDatabaseSynId(
            self.syn, "bed", databaseToSynIdMappingDf=databaseToSynIdMappingDf
        )
        bed = self.syn.tableQuery(
            "select Hugo_Symbol, ID from %s where CENTER = '%s'"
            % (bedSynId, self.center)
        )
        bedDf = bed.asDataFrame()
        fusion = fusion.apply(lambda x: validateSymbol(x, bedDf), axis=1)
        # Create nonmapped gene dict
        temp = fusion[fusion["HUGO_SYMBOL"] != fusion["ID"]]
        foo = temp[~temp.HUGO_SYMBOL.isnull()]
        temp = foo[["HUGO_SYMBOL", "ID"]]
        temp.drop_duplicates(inplace=True)
        temp.index = temp.ID
        del temp["ID"]
        # fusion = fusion[~fusion['HUGO_SYMBOL'].isnull()]
        fusion["FUSION"] = fusion["FUSION"].fillna("")
        fusion, nonmapped = remapFusion(temp.to_dict()["HUGO_SYMBOL"], fusion, "FUSION")
        # Fill in blank hugo symbol columns with original symbol
        null_symbols_idx = fusion["HUGO_SYMBOL"].isnull()
        fusion["HUGO_SYMBOL"][null_symbols_idx] = fusion["ID"][null_symbols_idx]
        # fusion, nonmapped = remapFusion(temp.to_dict()['HUGO_SYMBOL'], fusion, "COMMENTS")
        fusion["ENTREZ_GENE_ID"] = [int(float(i)) for i in fusion["ENTREZ_GENE_ID"]]
        return fusion

    # PROCESSING
    def process_steps(self, fusion, databaseSynId, newPath, databaseToSynIdMappingDf):
        fusion = self._process(fusion, databaseToSynIdMappingDf)
        process_functions.updateData(
            self.syn, databaseSynId, fusion, self.center, toDelete=True
        )
        fusion.to_csv(newPath, sep="\t", index=False)
        return newPath

    def _validate(self, fusionDF, nosymbol_check, project_id):
        total_error = ""
        warning = ""

        # Frame: "in-frame" or "frameshift".
        # Fusion_Status (OPTIONAL): An assessment of the mutation type (i.e., "SOMATIC", "GERMLINE", "UNKNOWN", or empty)

        fusionDF.columns = [col.upper() for col in fusionDF.columns]

        REQUIRED_HEADERS = pd.Series(
            [
                "HUGO_SYMBOL",
                "ENTREZ_GENE_ID",
                "CENTER",
                "TUMOR_SAMPLE_BARCODE",
                "FUSION",
                "DNA_SUPPORT",
                "RNA_SUPPORT",
                "METHOD",
                "FRAME",
            ]
        )
        if fusionDF.get("COMMENTS") is None:
            fusionDF["COMMENTS"] = float("nan")
        if not all(REQUIRED_HEADERS.isin(fusionDF.columns)):
            total_error += (
                "Your fusion file must at least have these headers: %s.\n"
                % ",".join(REQUIRED_HEADERS[~REQUIRED_HEADERS.isin(fusionDF.columns)])
            )
        if (
            process_functions.checkColExist(fusionDF, "HUGO_SYMBOL")
            and not nosymbol_check
        ):
            # logger.info("VALIDATING %s GENE SYMBOLS" % os.path.basename(filePath))
            # invalidated_genes = fusionDF["HUGO_SYMBOL"].drop_duplicates().apply(validateSymbol)
            databaseToSynIdMappingDf = process_functions.get_synid_database_mappingdf(
                self.syn, project_id
            )
            bedSynId = process_functions.getDatabaseSynId(
                self.syn, "bed", databaseToSynIdMappingDf=databaseToSynIdMappingDf
            )
            bed = self.syn.tableQuery(
                "select Hugo_Symbol, ID from %s where CENTER = '%s'"
                % (bedSynId, self.center)
            )
            bedDf = bed.asDataFrame()
            # invalidated_genes = self.pool.map(process_functions.validateSymbol, fusionDF["HUGO_SYMBOL"].drop_duplicates())
            if fusionDF["HUGO_SYMBOL"].isnull().any():
                total_error += (
                    "Your fusion file should not have any NA/blank Hugo Symbols.\n"
                )
            # fusionDF = fusionDF.drop_duplicates("HUGO_SYMBOL").apply(lambda x: validateSymbol(x, bedDf), axis=1)

        # if process_functions.checkColExist(fusionDF, "DNA_SUPPORT"):
        #     if not fusionDF.DNA_SUPPORT.isin(["yes","no","unknown"]).all():
        #         total_error += "Your fusion file's DNA_SUPPORT column must be 'yes', 'no', or 'unknown'"

        # if process_functions.checkColExist(fusionDF, "RNA_SUPPORT"):
        #     if not fusionDF.RNA_SUPPORT.isin(["yes","no","unknown"]).all():
        #         total_error += "Your fusion file's RNA_SUPPORT column must be 'yes', 'no', or 'unknown'"

        # if process_functions.checkColExist(fusionDF, "FRAME"):
        #     if not fusionDF.FRAME.isin(["in-frame","frameshift"]).all():
        #         total_error += "Your fusion file's FRAME column must be 'in-frame', or 'frameshift'"

        return (total_error, warning)
