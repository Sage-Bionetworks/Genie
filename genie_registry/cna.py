import logging
import os

import pandas as pd
import synapseclient

from genie.example_filetype_format import FileTypeFormat
from genie import process_functions

logger = logging.getLogger(__name__)


def validateSymbol(gene, bedDf, returnMappedDf=True):
    """
    Validate gene symbol

    Args:
        gene: Gene name
        bedDf: Bed pandas dataframe
        returnMappedDf: Return a mapped gene. Defaults to True

    Returns:
        gene name or boolean for whether a gene is valid
    """
    valid = False
    if sum(bedDf["Hugo_Symbol"] == gene) > 0:
        valid = True
    elif sum(bedDf["ID"] == gene) > 0:
        mismatch = bedDf[bedDf["ID"] == gene]
        mismatch.drop_duplicates(inplace=True)
        logger.info(
            "{} will be remapped to {}".format(gene, mismatch["Hugo_Symbol"].values[0])
        )
        gene = mismatch["Hugo_Symbol"].values[0]
    else:
        logger.warning(
            "{} cannot be remapped and will not be released. The symbol "
            "must exist in your seq assay ids (bed files) and must be "
            "mappable to a gene.".format(gene)
        )
        gene = float("nan")
    if returnMappedDf:
        return gene
    else:
        return valid


def makeCNARow(row, symbols):
    """
    Make CNA Row (Deprecated function)

    CNA values are no longer stored in the database

    Args:
        row: one row in the CNA file
        symbols:  list of Gene symbols
    """
    totalrow = "{symbols}\n{values}".format(
        symbols=",".join(symbols), values=",".join(row.astype(str))
    )
    totalrow = totalrow.replace(".0", "")
    return totalrow


def mergeCNAvalues(x):
    """Merge CNA values, make sure if there are two rows that are the
    same gene, the values are merged"""
    # Change into its own series, because sometimes doing an apply
    # will cause there to be a missing index value which will
    # cause dropna() to fail.
    values = pd.Series(x.values)
    values.dropna(inplace=True)
    uniqueValues = set(values.unique())
    if len(uniqueValues) == 1:
        returnVal = x.tolist()[0]
    elif len(uniqueValues) <= 2:
        uniqueValues.discard(0)
        if len(uniqueValues) == 1:
            returnVal = list(uniqueValues)[0]
        else:
            returnVal = float("nan")
    else:
        returnVal = float("nan")
    return returnVal


def checkIfOneZero(x):
    assert len(set(x.tolist())) == 1, "Can only be one unique value"


class cna(FileTypeFormat):

    _fileType = "cna"

    _process_kwargs = ["newPath"]

    _validation_kwargs = ["nosymbol_check"]

    # VALIDATE FILENAME
    def _validateFilename(self, filePath):
        assert os.path.basename(filePath[0]) == "data_CNA_{}.txt".format(self.center)

    def _process(self, cnaDf):
        cnaDf.rename(columns={cnaDf.columns[0]: cnaDf.columns[0].upper()}, inplace=True)
        cnaDf.rename(columns={"HUGO_SYMBOL": "Hugo_Symbol"}, inplace=True)

        index = [
            i for i, col in enumerate(cnaDf.columns) if col.upper() == "ENTREZ_GENE_ID"
        ]
        if len(index) > 0:
            del cnaDf[cnaDf.columns[index][0]]

        bedSynId = self.genie_config["bed"]
        bed = self.syn.tableQuery(
            f"select Hugo_Symbol, ID from {bedSynId} where CENTER = '{self.center}'"
        )
        bedDf = bed.asDataFrame()
        cnaDf["Hugo_Symbol"] = cnaDf["Hugo_Symbol"].apply(
            lambda x: validateSymbol(x, bedDf)
        )
        order = cnaDf.columns
        cnaDf = cnaDf[~cnaDf["Hugo_Symbol"].isnull()]
        # cnaDf = cnaDf.applymap(str)
        duplicatedGenes = pd.DataFrame()
        duplicated_symbols = cnaDf["Hugo_Symbol"][
            cnaDf["Hugo_Symbol"].duplicated()
        ].unique()
        for i in duplicated_symbols:
            dups = cnaDf[cnaDf["Hugo_Symbol"] == i]
            newVal = dups[dups.columns[dups.columns != "Hugo_Symbol"]].apply(
                mergeCNAvalues
            )
            temp = pd.DataFrame(newVal).transpose()
            temp["Hugo_Symbol"] = i
            duplicatedGenes = pd.concat([duplicatedGenes, temp], sort=False)
        cnaDf.drop_duplicates("Hugo_Symbol", keep=False, inplace=True)
        cnaDf = pd.concat([cnaDf, duplicatedGenes], sort=False)
        cnaDf = cnaDf[order]
        return cnaDf

    def process_steps(self, cnaDf, newPath):
        newCNA = self._process(cnaDf)

        centerMafSynId = self.genie_config["centerMaf"]
        if not newCNA.empty:
            cnaText = process_functions.removePandasDfFloat(newCNA)
            # Replace blank with NA's
            cnaText = (
                cnaText.replace("\t\t", "\tNA\t")
                .replace("\t\t", "\tNA\t")
                .replace("\t\n", "\tNA\n")
            )
            with open(newPath, "w") as cnaFile:
                cnaFile.write(cnaText)
            self.syn.store(synapseclient.File(newPath, parent=centerMafSynId))
        return newPath

    def _validate(self, cnvDF, nosymbol_check):
        total_error = ""
        warning = ""
        cnvDF.columns = [col.upper() for col in cnvDF.columns]

        if cnvDF.columns[0] != "HUGO_SYMBOL":
            total_error += "Your cnv file's first column must be Hugo_Symbol\n"
        haveColumn = process_functions.checkColExist(cnvDF, "HUGO_SYMBOL")
        if haveColumn:
            keepSymbols = cnvDF["HUGO_SYMBOL"]
            cnvDF.drop("HUGO_SYMBOL", axis=1, inplace=True)

        # if sum(cnvDF.apply(lambda x: sum(x.isnull()))) > 0:
        #   total_error += "Your cnv file must not have any empty values\n"

        if process_functions.checkColExist(cnvDF, "ENTREZ_GENE_ID"):
            del cnvDF["ENTREZ_GENE_ID"]
        error = process_functions.validate_genie_identifier(
            identifiers=cnvDF.columns, center=self.center, filename="cnv", col="samples"
        )
        total_error += error
        # cnvDF = cnvDF.fillna('')
        allowed_values = [
            "-2.0",
            "-2",
            "-1.5",
            "-1.0",
            "-1",
            "0.0",
            "0",
            "0.5",
            "1.0",
            "1",
            "1.5",
            "2",
            "2.0",
            "nan",
        ]
        if not all(cnvDF.applymap(lambda x: str(x) in allowed_values).all()):
            total_error += (
                "All values must be NA/blank, -2, -1.5, -1, -0.5, "
                "0, 0.5, 1, 1.5, or 2.\n"
            )
        else:
            cnvDF["HUGO_SYMBOL"] = keepSymbols
            if haveColumn and not nosymbol_check:
                bedSynId = self.genie_config["bed"]
                bed = self.syn.tableQuery(
                    f"select Hugo_Symbol, ID from {bedSynId} "
                    f"where CENTER = '{self.center}'"
                )
                bedDf = bed.asDataFrame()
                cnvDF["remapped"] = cnvDF["HUGO_SYMBOL"].apply(
                    lambda x: validateSymbol(x, bedDf)
                )
                cnvDF = cnvDF[~cnvDF["remapped"].isnull()]

                # Do not allow any duplicated genes after symbols
                # have been remapped
                if sum(cnvDF["remapped"].duplicated()) > 0:
                    duplicated = cnvDF["remapped"].duplicated(keep=False)
                    total_error += (
                        "Your CNA file has duplicated Hugo_Symbols "
                        "(After remapping of genes): {} -> {}.\n".format(
                            ",".join(cnvDF["HUGO_SYMBOL"][duplicated]),
                            ",".join(cnvDF["remapped"][duplicated]),
                        )
                    )
        return (total_error, warning)
