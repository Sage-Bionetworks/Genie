from __future__ import absolute_import
from genie import example_filetype_format, sampleRetraction
import logging
import os
logger = logging.getLogger(__name__)

class patientRetraction(sampleRetraction.sampleRetraction):

    _fileType = "patientRetraction"