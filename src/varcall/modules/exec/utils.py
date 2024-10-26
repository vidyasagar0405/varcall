import glob
from datetime import datetime
import os.path

from pathlib import Path

from textual.suggester import Suggester
from textual.widgets import Input

# from varcall.modules.widgets import YesOrNo
# from varcall.main import Varcall



def get_random_file_name_current(ext:str) -> str:

    now = datetime.now()
    date_time = now.strftime("%H%M_%d%m%Y")
    path = os.getcwd()
    path = f"{path}/{date_time}.{ext}"

    return path

def default_file_name(self, tool:str, input:str , has_cwd_set: bool) -> str:

    input_cwd_path_tuple = get_basename_and_ext(input)
    f"{input_cwd_path_tuple[0]}.sorted{input_cwd_path_tuple[1]}"

    if has_cwd_set:
        cwd_path = self.workingDir

    cwd_path = os.getcwd()
    now = datetime.now()
    date_time = now.strftime("%H%M_%d%m%Y")

    default_path = {
        "download": f"{cwd_path}/{date_time}.download",

        "fastqc": f"{cwd_path}results/fastqc/",
        "multiqc": f"{cwd_path}results/multiqc{date_time}/",

        "bwa_index": f"{cwd_path}data/reference/",
        "bwa_mem": f"{cwd_path}results/sam/aligned{date_time}.sam",

        "samtools_sort": f"{cwd_path}results/sam/aligned{date_time}.sorted.sam",
        "samtools_view": f"{cwd_path}results/bam/aligned{date_time}.bam",
        "samtools_index": f"{cwd_path}results/sam/",

        "bcf_call": f"{cwd_path}results/bcf/{date_time}.call.bcf",
        "bcf_filter": f"{cwd_path}results/multiqc/",
        "bcf_norm": f"{cwd_path}results/multiqc/",
        "bcf_stats": f"{cwd_path}results/multiqc/",
    }
    return default_path[tool] if True else f"Error: '{self.name}' command not found."



def get_output_name(self, tool: str, input) -> str:
    if not self.workingDir:
        return default_file_name(self, tool, input, True)
    else:
        return default_file_name(self, tool, input, False)



class FileSuggester(Suggester):
    async def get_suggestion(self, value: str) -> str | None:
        """
        Provides file path suggestions based on the input value.
        """
        # Handle absolute paths
        if value.startswith("/"):
            return " "
        # Handle invalid pattern '**' anywhere in the input
        if "**" in value:
            return " "
        path = next(Path().glob(f"{value}*"), None)
        return str(path) if path else None


file_suggester = FileSuggester(use_cache=False, case_sensitive=True)


def get_basename_and_ext(path) -> tuple:
    """returns a tuple with basename and extension (has only two elements) e.g: ("basename", ".ext")
    returns ("", "") if input is false"""
    filename = os.path.splitext(os.path.basename(path))
    return filename

def get_file_extension(path) -> str:
    filename = os.path.splitext(os.path.basename(path))
    return filename[-1]


def find_matching_files(path_pattern: str) -> str:
    """Find all files matching the given wildcard path pattern."""
    matching_files = glob.glob(path_pattern)
    return " ".join(matching_files)


def get_input(self, input_widget_id) -> str:
    try:
        raw_input = str(self.query_one(f"#{input_widget_id}", Input).value).strip()
        return find_matching_files(raw_input)
    except Exception:
        return ""



