######################################################################
# Main app information.
__author__ = "Vidyasagar"
__copyright__ = "Copyright 2024, Vidyasagar"
__credits__ = ["Vidyasagar"]
__maintainer__ = "Vidyasagar"
__version__ = "0.0.2"
__licence__ = "MIT"

##############################################################################
# Local imports.
from .logging import setup_logging

from .exec.common import get_input, find_matching_files, file_suggester
from .exec.hometab  import run_FastQC, run_MultiQC, run_bwa_mem, run_download, run_bwa_index
from .exec.bcftools import run_bcf_call, run_bcf_norm, run_bcf_stats, run_bcf_filter, run_bcf_mpileup
from .exec.samtools  import run_samtools_view, run_samtools_sort, run_samtools_index

from .widgets.Home import HomeWidgets
from .widgets.Samtools import SamtoolsWidgets
from .widgets.Bcftools import BcftoolsWidgets
from .widgets.Pipeline import PipelineWidgets

from .Help.Help import HelpMarkdown

##############################################################################
# Export the imports.
__all__ = [
    "setup_logging",

    "get_input", "find_matching_files", "file_suggester",
    "run_bcf_call", "run_bcf_norm", "run_bcf_stats", "run_bcf_filter", "run_bcf_mpileup", 
    "run_bwa_index", "run_bwa_mem", "run_download", "run_MultiQC", "run_FastQC",
    "run_samtools_view", "run_samtools_sort", "run_samtools_index",

    "HomeWidgets",
    "SamtoolsWidgets",
    "BcftoolsWidgets",
    "PipelineWidgets",

    "HelpMarkdown",
]

### __init__.py ends here
