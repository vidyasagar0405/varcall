import threading
import subprocess
import logging
from pathlib import Path
from modules.logging import setup_logging
from textual.widgets import Input
from textual.suggester import Suggester

setup_logging()

class FileSuggester(Suggester):
    async def get_suggestion(self, value: str) -> str | None:
        # Handle absolute paths
        if value.startswith('/'):
            return " "
        # Handle invalid pattern '**' anywhere in the input
        if '**' in value:
            return " "
        path = next(Path().glob(f"{value}*"), None)
        return str(path) if path else None

def run_download(self):
        raw_url = self.query_one("#Input_url", Input)
        download_url = raw_url.value.strip()
        if not download_url:
            self.notify("Please provide a valid url", severity="warning", title="Download")
            return
        elif not self.full_output_path:
            self.notify("Please provide a valid url", severity="warning", title="Download")
            return
        else:
            self.notify(f"Downloading {str(download_url)}...", title="Download")
            logging.info(f"Downloading {str(download_url)}...")
            self.query_one("#Download_options").add_class("running")
            threading.Thread(target=_run_download, args=(self, download_url,)).start()

def _run_download(self, download_url: str) -> None:
    try:
        download_cmd = f"curl -L {str(download_url)} -o {self.full_output_path}"
        self.notify(download_cmd, title="Download")
        logging.info(f"Running command: " + download_cmd)
        result = subprocess.run(
            [download_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True
        )
        logging.info(f"Downloaded {str(download_url)}")
        self.notify(f"Downloaded {str(download_url)}", title="Download")
        self.query_one("#Download_options").remove_class("running")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while Downloading: {e.stderr}")
        self.notify(f"An error occurred while Downloading: {e.stderr}", severity="error", title="Download")
        self.query_one("#Download_options").remove_class("running")

def run_FastQC(self):
        input_path = self.query_one("#FastQC_Input", Input)
        input_path = input_path.value.strip()
        if not input_path:
            input_path = f"{self.workingDir}/data/reads/*"
            self.notify(f"FastQC {str(input_path)}...", title="FastQC")
            logging.info(f"FastQC {str(input_path)}...")
            self.query_one("#FastQC_Horizontal").add_class("running")
            threading.Thread(target=_run_FastQC, args=(self, input_path,)).start()
        else:
            self.notify(f"FastQC {str(input_path)}...", title="FastQC")
            logging.info(f"FastQC {str(input_path)}...")
            self.query_one("#FastQC_Horizontal").add_class("running")
            threading.Thread(target=_run_FastQC, args=(self, input_path,)).start()

def _run_FastQC(self, input_path: str) -> None:
    try:
        FastQC_cmd = f"fastqc {str(input_path)} -o {self.workingDir}/results/fastqc/"
        self.notify(FastQC_cmd, title="FastQC")
        logging.info(f"Running command: " + FastQC_cmd)
        result = subprocess.run(
            [FastQC_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True
        )
        logging.info(f"FastQC completed for {str(input_path)}")
        self.notify(f"FastQC completed for {str(input_path)}", title="FastQC")
        self.query_one("#FastQC_Horizonal").remove_class("running")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during FastQC: {e.stderr}")
        self.notify(f"An error occurred during FastQc: {e.stderr}", severity="error", title="FastQC")
        self.query_one("#FastQC_Horizonal").remove_class("running")

def run_MultiQC(self):
        input_path = self.query_one("#MultiQC_Input", Input)
        input_path = input_path.value.strip()
        if not input_path:
            input_path = f"{self.workingDir}/results/fastqc/"

        self.notify(f"MultiQC {str(input_path)}...", title="MultiQC")
        logging.info(f"MultiQC {str(input_path)}...")
        self.query_one("#MultiQC_Horizontal").add_class("running")
        threading.Thread(target=_run_MultiQC, args=(self, input_path,)).start()

def _run_MultiQC(self, input_path: str) -> None:
    try:
        MultiQC_cmd = f"multiqc {str(input_path)} -o {self.workingDir}/results/multiqc/"
        self.notify(MultiQC_cmd, title="MultiQC")
        logging.info(f"Running command: " + MultiQC_cmd)
        result = subprocess.run(
            [MultiQC_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True
        )
        logging.info(f"MultiQC completed for {str(input_path)}")
        self.notify(f"MultiQC completed for {str(input_path)}", title="MultiQC")
        self.query_one("#MultiQC_Horizonal").remove_class("running")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during MultiQC: {e.stderr}")
        self.notify(f"An error occurred during MultiQC: {e.stderr}", severity="error", title="MultiQC")
        self.query_one("#MultiQC_Horizonal").remove_class("running")

def run_bwa_index(self):
        input_path = self.query_one("#bwa_ref_Input", Input)
        input_path = input_path.value.strip()
        if not input_path:
            input_path = f"{self.workingDir}/data/reference/*"
        self.notify(f"bwa indexing {str(input_path)}...", title="bwa index")
        logging.info(f"bwa indexing {str(input_path)}...")
        self.query_one("#bwa_Horizontal").add_class("running")
        threading.Thread(target=_run_bwa_index, args=(self, input_path,)).start()

def _run_bwa_index(self, input_path: str) -> None:
    try:
        bwa_index_cmd = f"bwa index {str(input_path)}"
        self.notify(bwa_index_cmd, title="bwa Index")
        logging.info(f"Running command: " + bwa_index_cmd)
        result = subprocess.run(
            [bwa_index_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True
        )
        logging.info(f"bwa Index completed for {str(input_path)}")
        self.notify(f"bwa Index completed for {str(input_path)}", title="bwa Index")
        self.query_one("#bwa_Horizonal").remove_class("running")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during bwa Index: {e.stderr}")
        self.notify(f"An error occurred during bwa Index: {e.stderr}", severity="error", title="bwa Index")
        self.query_one("#bwa_Horizonal").remove_class("running")

def run_bwa_mem(self):
        read_path = self.query_one("#bwa_reads_Input", Input)
        ref_path = self.query_one("#bwa_ref_Input", Input)
        read_path = read_path.value.strip()
        ref_path = ref_path.value.strip()
        if read_path != 2:
            self.notify("Enter two sample reads file path", title="bwa mem")
        if not ref_path:
            self.notify("Enter reference genome file path", title="bwa mem")
        self.notify(f"aligning reads using bwa mem {str(read_path)}...", title="bwa mem")
        logging.info(f"aligning reads using bwa mem {str(read_path)}...")
        self.query_one("#bwa_Horizontal").add_class("running")
        threading.Thread(target=_run_bwa_index, args=(self, ref_path, read_path,)).start()

def _run_bwa_index(self, ref_path: str, read_path: str) -> None:
    try:
        bwa_mem_cmd = f"bwa mem {str(ref_path)} {str(read_path)} -o trial/results/sam/aligned.sam"
        self.notify(bwa_mem_cmd, title="bwa mem")
        logging.info(f"Running command: " + bwa_mem_cmd)
        result = subprocess.run(
            [bwa_mem_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True
        )
        logging.info(f"Completed aligning reads using bwa mem {str(read_path)}")
        self.notify(f"Completed aligning reads using bwa mem {str(read_path)}", title="bwa mem")
        self.query_one("#bwa_Horizonal").remove_class("running")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error aligning reads using bwa mem: {e.stderr}")
        self.notify(f"Error aligning reads using bwa mem: {e.stderr}", severity="error", title="bwa mem")
        self.query_one("#bwa_Horizonal").remove_class("running")
