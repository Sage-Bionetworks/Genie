# Import logging last to not take in synapseclient logging
import logging

from . import process_functions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("genie")
