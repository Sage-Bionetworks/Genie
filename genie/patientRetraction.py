from __future__ import absolute_import
from genie import example_filetype_format
from genie import sampleRetraction

import logging
import os
#import pandas as pd
#import synapseclient
#import datetime
#logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class patientRetraction(sampleRetraction.sampleRetraction):

	_fileType = "patientRetraction"