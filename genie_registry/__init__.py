"""Initialize GENIE registry"""
# Import logging last to not take in synapseclient logging
import logging

from . import (
    assay,
    bed,
    cbs,
    clinical,
    cna,
    fusions,
    maf,
    mutationsInCis,
    patientRetraction,
    sampleRetraction,
    seg,
    structural_variant,
    vcf,
    workflow,
)

# Set genie registry to be the same as main genie codebase
from genie import __version__

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logging.getLogger("keyring").setLevel(logging.WARNING)

__all__ = [
    "__version__",
    "assay",
    "bed",
    "cbs",
    "clinical",
    "cna",
    "fusions",
    "maf",
    "mutationsInCis",
    "patientRetraction",
    "sampleRetraction",
    "seg",
    "structural_variant",
    "vcf",
    "workflow",
]
