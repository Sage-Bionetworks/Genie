import logging

from .sampleRetraction import sampleRetraction

logger = logging.getLogger(__name__)


class patientRetraction(sampleRetraction):

    _fileType = "patientRetraction"
