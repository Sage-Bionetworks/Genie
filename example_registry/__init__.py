"""Initialize GENIE registry"""
# Import logging last to not take in synapseclient logging
import logging

from . import bed
from . import vcf
from . import workflow
from . import clinical
from . import seg
from . import cbs
from . import maf
from . import cna
from . import fusions
from . import sample_retraction
from . import patient_retraction
from . import mutations_in_cis
from . import assay
from . import validate

from .__version__ import __version__

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
