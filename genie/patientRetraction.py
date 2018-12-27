from __future__ import absolute_import
from genie import sampleRetraction
import logging
import os
logger = logging.getLogger(__name__)

class patientRetraction(sampleRetraction):

	_fileType = "patientRetraction"