from __future__ import absolute_import

from . import bed
from . import bedSP
from . import workflow
from . import clinical
from . import seg
from . import cbs
from . import maf
from . import mafSP
from . import clinicalSP
from . import vcf
from . import cna
from . import fusions
from . import sampleRetraction
from . import patientRetraction
from . import patientCounts
from . import mutationsInCis
from . import vitalStatus
from . import process_functions
from . import example_filetype_format
from . import create_case_lists

PROCESS_FILES = {'bed':bed.bed,
				 'bedSP':bedSP.bedSP,
				 'maf':maf.maf,
				 'mafSP':mafSP.mafSP,
				 'clinical':clinical.clinical,
				 'clinicalSP':clinicalSP.clinicalSP,
				 'vcf':vcf.vcf,
				 'cbs':cbs.cbs,
				 'cna':cna.cna,
				 'fusions':fusions.fusions,
				 'md':workflow.workflow,
				 'seg':seg.seg,
				 'patientRetraction':patientRetraction.patientRetraction,
				 'sampleRetraction':sampleRetraction.sampleRetraction,
				 'patientCounts':patientCounts.patientCounts,
				 'mutationsInCis':mutationsInCis.mutationsInCis,
				 'vitalStatus':vitalStatus.vitalStatus}
#Must import validate after the fact
from . import validate
