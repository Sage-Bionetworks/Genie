"""Initialize GENIE registry"""
# Import logging last to not take in synapseclient logging
import logging

from . import bed
from . import vcf
from . import bedSP
from . import workflow
from . import clinical
from . import seg
from . import cbs
from . import maf
from . import mafSP
from . import clinicalSP
from . import cna
from . import fusions
from . import sampleRetraction
from . import patientRetraction
from . import mutationsInCis
from . import assay

from genie.__version__ import __version__

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logging.getLogger("keyring").setLevel(logging.WARNING)
