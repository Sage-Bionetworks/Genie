# Import logging last to not take in synapseclient logging
import logging

from . import vcf
from . import clinical
from . import cna
from . import maf
from . import mafSP
from . import clinicalSP
from . import mutationsInCis

from .__version__ import __version__

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
