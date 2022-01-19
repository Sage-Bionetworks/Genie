import logging
import os

import pandas as pd

logger = logging.getLogger(__name__)


class FileTypeFormat(object):

    _process_kwargs = ["newPath", "databaseSynId"]

    _fileType = "fileType"

    _validation_kwargs = []

    def __init__(self, syn, center, poolSize=1):
        self.syn = syn
        self.center = center

        # self.pool = multiprocessing.Pool(poolSize)

    def _get_dataframe(self, filePathList):
        """
        This function by defaults assumes the filePathList is length of 1
        and is a tsv file.  Could change depending on file type.

        Args:
            filePathList:  A list of file paths (Max is 2 for the two
                           clinical files)

        Returns:
            df: Pandas dataframe of file
        """
        filePath = filePathList[0]
        df = pd.read_csv(filePath, sep="\t", comment="#")
        return df

    def read_file(self, filePathList):
        """
        Each file is to be read in for validation and processing.
        This is not to be changed in any functions.

        Args:
            filePathList:  A list of file paths (Max is 2 for the two
                           clinical files)

        Returns:
            df: Pandas dataframe of file
        """
        df = self._get_dataframe(filePathList)
        return df

    def _validateFilename(self, filePath):
        """
        Function that changes per file type for validating its filename
        Expects an assertion error.

        Args:
            filePath: Path to file
        """
        # assert True
        raise NotImplementedError

    def validateFilename(self, filePath):
        """
        Validation of file name.  The filename is what maps the file
        to its validation and processing.

        Args:
            filePath: Path to file

        Returns:
            str: file type defined by self._fileType
        """
        self._validateFilename(filePath)
        return self._fileType

    def process_steps(self, df, **kwargs):
        """
        This function is modified for every single file.
        It reformats the file and stores the file into database and Synapse.
        """
        pass

    def preprocess(self, newpath):
        """
        This is for any preprocessing that has to occur to the entity name
        to add to kwargs for processing.  entity name is included in
        the new path

        Args:
            newpath: Path to file
        """
        return dict()

    def process(self, filePath, **kwargs):
        """
        This is the main processing function.

        Args:
            filePath: Path to file
            kwargs: The kwargs are determined by self._process_kwargs

        Returns:
            str: file path of processed file
        """
        preprocess_args = self.preprocess(kwargs.get("newPath"))
        kwargs.update(preprocess_args)
        mykwargs = {}
        for required_parameter in self._process_kwargs:
            assert required_parameter in kwargs.keys(), (
                "%s not in parameter list" % required_parameter
            )
            mykwargs[required_parameter] = kwargs[required_parameter]
        logger.info("PROCESSING %s" % filePath)
        # This is done because the clinical files are being merged into a list
        if self._fileType == "clinical":
            path_or_df = self.read_file(filePath)
        # If file type is vcf or maf file, processing requires a filepath
        elif self._fileType not in ["vcf", "maf", "mafSP", "md"]:
            path_or_df = self.read_file([filePath])
        else:
            path_or_df = filePath
        path = self.process_steps(path_or_df, **mykwargs)
        return path

    def _validate(self, df, **kwargs):
        """
        This is the base validation function.
        By default, no validation occurs.

        Args:
            df: A dataframe of the file
            kwargs: The kwargs are determined by self._validation_kwargs

        Returns:
            tuple: The errors and warnings as a file from validation.
                   Defaults to blank strings
        """
        errors = ""
        warnings = ""
        logger.info("NO VALIDATION for %s files" % self._fileType)
        return (errors, warnings)

    def validate(self, filePathList, **kwargs):
        """
        This is the main validation function.
        Every file type calls self._validate, which is different.

        Args:
            filePathList: A list of file paths.
            kwargs: The kwargs are determined by self._validation_kwargs

        Returns:
            tuple: The errors and warnings as a file from validation.
        """
        mykwargs = {}
        for required_parameter in self._validation_kwargs:
            assert required_parameter in kwargs.keys(), (
                "%s not in parameter list" % required_parameter
            )
            mykwargs[required_parameter] = kwargs[required_parameter]

        errors = ""

        try:
            df = self.read_file(filePathList)
        except Exception as e:
            errors = "The file(s) ({filePathList}) cannot be read. Original error: {exception}".format(
                filePathList=filePathList, exception=str(e)
            )
            warnings = ""

        if not errors:
            logger.info("VALIDATING %s" % os.path.basename(",".join(filePathList)))
            errors, warnings = self._validate(df, **mykwargs)

        valid = errors == ""

        return valid, errors, warnings
