import glob
from datetime import datetime
import os.path

import logging
import subprocess
import threading
from pathlib import Path

from varcall.modules.logging import setup_logging
from textual.suggester import Suggester
from textual.widgets import Input

setup_logging()


def get_random_file_name_current(ext: str) -> str:

    now = datetime.now()
    date_time = now.strftime("%H%M_%d%m%Y")
    path = os.getcwd()
    path = f"{path}/{date_time}.{ext}"

    return path


def default_file_name(self, tool: str, input: str, has_cwd_set: bool) -> str:

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


class Command:
    def __init__(
        self,
        workingDir: str,
        name: str,
        first_path: str,
        second_path: str = "",
        third_input: str = "",
        fourth_input: str = "",
        fifth_input=None,
    ) -> None:
        self.name = name
        self.workingDir = workingDir
        self.first_path = first_path
        self.second_path = second_path
        self.third_input = third_input
        self.fourth_input = fourth_input
        self.fifth_input = fifth_input

    def get_command(self) -> str:
        # Generate commands dynamically based on instance variables
        commands = {
            "download": f"curl -L {self.first_path} -o {self.second_path or get_output_name(self, "download", self.first_path)}",
            "fastqc": f"fastqc {self.first_path} -o {self.second_path or get_output_name(self, "fastqc", self.first_path)}",
            "multiqc": f"multiqc {self.first_path} -o {self.second_path or get_output_name(self, "multiqc", self.first_path)}",
            "bwa_index": f"bwa index {self.first_path or get_output_name(self, "bwa_index", self.first_path)}",
            "bwa_mem": f"bwa mem -t {self.second_path or 4} {self.first_path} {self.third_input} {self.fourth_input} -o {self.fifth_input or get_output_name(self, "bwa_mem", self.first_path)}",
            "samtools_sort": f"samtools sort {self.first_path} -o {self.second_path or get_output_name(self, "samtools_sort", self.first_path)}",
            "samtools_view": f"samtools view -b -h {self.first_path} -o {self.second_path or get_output_name(self, "samtools_view", self.first_path)} region={self.third_input}",
            "samtools_index": f"samtools index {self.first_path or get_output_name(self, "samtools_index", self.first_path)}",
            "bcf_call": f"bcftools call -mv -Oz -o {self.second_path or get_output_name(self, "bcf_call", self.first_path)} {self.first_path}",
            "bcf_filter": f"bcftools filter -s LOWQUAL -e '{self.third_input}' {self.first_path} -o {self.second_path or get_output_name(self, "bcf_filter", self.first_path)}",
            "bcf_norm": f"bcftools norm -f {self.first_path} {self.second_path or get_output_name(self, "bcf_norm", self.first_path)}",
            "bcf_stats": f"bcftools stats {self.first_path} > {self.second_path or get_output_name(self, "bcf_stats", self.first_path)}",
        }
        return commands.get(self.name, f"Error: '{self.name}' command not found.")


def run_command(
    self,
    name: str,
    first_path: str,
    second_path: str = "",
    third_input: str = "",
    fourth_input: str = "",
    fifth_input: str = "",
):

    _name = name.replace(" ", "_")

    first_path = get_input(self, f"{_name}_first_path")
    second_path = get_input(self, f"{_name}_second_path")
    third_input = get_input(self, f"{_name}_third_input")
    fourth_input = get_input(self, f"{_name}_fourth_inpu")
    fifth_input = get_input(self, f"{_name}_fifth_input")

    if not first_path and (name != "fastqc" or "multiqc" or "bwa_index"):
        self.notify("Please provide a valid path", severity="warning", title=name)
        return

    self.notify(f"{name} {str(first_path )}...", title=name)
    logging.info(f"{name} {str(first_path )}...")
    self.query_one(f"#{_name}_horizontal").add_class("running")
    cmd = Command(
        name,
        first_path,
        second_path,
        third_input,
        fourth_input,
        fifth_input,
    ).get_command()
    threading.Thread(
        target=_run_command,
        args=(
            self,
            cmd,
            name,
        ),
    ).start()


def _run_command(self, cmd: str, name: str) -> None:
    """ """
    _name = name.replace(" ", "_")
    try:
        self.notify(cmd, title=name)
        logging.info("Running command: " + cmd)
        subprocess.run(
            cmd.split(),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        logging.info(f"{name} completed")
        self.notify(f"{name} completed", title=name)
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during {name}: {e.stderr}")
        self.notify(
            f"An error occurred during {name}: {e.stderr}",
            severity="error",
            timeout=10.0,
            title=name,
        )

    finally:
        self.query_one(f"#{_name}_horizontal").remove_class("running")
