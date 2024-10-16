from textual.widgets import Input
from textual.suggester import Suggester
from pathlib import Path
import glob
# from ...src.main import Varcall

import logging
import subprocess
import threading

from modules.logging import setup_logging

setup_logging()


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


def find_matching_files(path_pattern: str) -> str:
    """Find all files matching the given wildcard path pattern."""
    matching_files = glob.glob(path_pattern)
    return " ".join(matching_files)


def get_input(self, input_widget_id) -> str:
    try:
        raw_input = str(self.query_one(f"#{input_widget_id}", Input).value).strip()
        return find_matching_files(raw_input)
    except:
        return ""


class Command:
    def __init__(self, name: str, input_path: str, output_path: str = "", third_input: str = "", fourth_input: str = "", fifth_input = None) -> None:
        self.name = name
        self.input_path = input_path
        self.output_path = output_path
        self.third_input = third_input
        self.fourth_input = fourth_input
        self.fifth_input = fifth_input

    def get_command(self) -> str:
    # Generate commands dynamically based on instance variables
        commands = {
            
            "download": f"curl -L {self.input_path} -o {self.output_path}",
            "fastqc": f"fastqc {self.input_path} -o {self.output_path}",
            "multiqc": f"multiqc {self.input_path} -o {self.output_path}",
            "bwa_index": f"bwa index {self.input_path}",

            # no_of_threads → third_input
            # ref_path → fourth_input
            # read_path_2 → fifth_input
            # view_region → third_input
            # filter → third_input

            "bwa_mem": f"bwa mem -t {self.third_input} {self.input_path} {self.fourth_input} {self.fifth_input} -o {self.output_path}",
            "multiqc_flagstat": f"multiqc {self.third_input}/results/sam/aligned.sam.flagstat -o {self.third_input}/results/multiqc/sam",
            "samtools_sort": f"samtools sort {self.input_path} -o {self.output_path}",
            "samtools_view": f"samtools view -b -h {self.input_path} -o {self.output_path} region={self.third_input}",
            "samtools_index": f"samtools index {self.input_path}",
            "bcf_call": f"bcftools call -mv -Oz -o {self.output_path} {self.input_path}",
            "bcf_filter": f"bcftools filter -s LOWQUAL -e '{self.third_input}' {self.input_path} -o {self.output_path}",
            "bcf_norm": f"bcftools norm -f {self.input_path} {self.output_path}",
            "bcf_stats": f"bcftools stats {self.input_path} > {self.output_path}"
        }
        return commands.get(self.name, f"Error: '{self.name}' command not found.")




def run_command(self, name: str, input_path: str, output_path: str = "", third_input: str = "", fourth_input: str = "", fifth_input = None):
    """ """
    input_path = get_input(self, "bcf_call_input_input")
    output_path = get_input(self, "bcf_call_output_input")
    if not input_path:
        self.notify(
            "Please provide a valid path", severity="warning", title="Samtools view"
        )
        return
    if not output_path:
        self.notify(
            "Please provide a valid path", severity="warning", title="Samtools view"
        )
        return

    self.notify(f"bcf call {str(input_path)}...", title="bcf call")
    logging.info(f"bcf call {str(input_path)}...")
    self.query_one("#bcf_call_horizontal").add_class("running")
    bcf_call_cmd = f"bcftools call -mv -Oz -o {output_path} {input_path}"
    threading.Thread(
        target=_run_command,
        args=(
            self,
            cmd,
            name,
        ),
    ).start()


def _run_command(self, cmd: str, name:str) -> None:
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
