import logging
import os
import pandas as pd
import process_functions
import synapseclient
import datetime
import sampleRetraction
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class patientRetraction(sampleRetraction.sampleRetraction):

	_fileType = "patientRetraction"