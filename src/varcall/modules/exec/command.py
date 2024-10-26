import logging
import subprocess
import threading
from varcall.modules.logging import setup_logging
from varcall.modules.exec.utils import get_input, get_output_name


setup_logging()

class Command:
    def __init__(self, workingDir: str, name: str, first_path: str, second_path: str = "", third_input: str = "", fourth_input: str = "", fifth_input = None) -> None:
        # super().__init__()
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
            "bcf_stats": f"bcftools stats {self.first_path} > {self.second_path or get_output_name(self, "bcf_stats", self.first_path)}"
        }
        return commands.get(self.name, f"Error: '{self.name}' command not found.")




def run_command(self, Command: Command):

    _name = Command.name.replace(" ", "_")

    first_path = get_input(self, f"{_name}_first_path")
    second_path = get_input(self, f"{_name}_second_path")
    third_input = get_input(self, f"{_name}_third_input")
    fourth_input  = get_input(self, f"{_name}_fourth_inpu")
    fifth_input = get_input(self, f"{_name}_fifth_input")

    if not first_path and (Command.name != ("fastqc" or "multiqc" or "bwa_index")):
        self.notify( "Please provide a valid path", severity="warning", title=Command.name)
        return

    self.notify(f"{Command.name} {str(first_path )}...", title=Command.name)
    logging.info(f"{Command.name} {str(first_path )}...")
    self.query_one(f"#{_name}_horizontal").add_class("running")
    cmd = Command(
        Command.name,
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
            Command.name,
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

