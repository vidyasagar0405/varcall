import threading
import subprocess
import logging
from modules.logging import setup_logging
from modules.exec_func.samtools_funcs import flagstat_bam
from modules.exec_func.common_funcs import get_input

setup_logging()


def run_download(self):
    """
    Initiates the download process for a given URL.

    This function retrieves the URL from the input field, validates it,
    and starts a new thread to handle the download process.

    Args:
        self: The instance of the class calling this function.
    """
    download_url = get_input(self, "Input_url")
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
        threading.Thread(
            target=_run_download,
            args=(
                self,
                download_url,
            ),
        ).start()


def _run_download(self, download_url: str) -> None:
    """
    Executes the download command using subprocess.

    This function runs the curl command to download the file from the given URL
    and saves it to the specified output path.

    Args:
        self: The instance of the class calling this function.
        download_url (str): The URL from which to download the file.
    """
    try:
        download_cmd = f"curl -L {str(download_url)} -o {self.full_output_path}"
        self.notify(download_cmd, title="Download")
        logging.info("Running command: " + download_cmd)
        subprocess.run(
            [download_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"Downloaded {str(download_url)}")
        self.notify(f"Downloaded {str(download_url)}", title="Download")
        self.query_one("#Download_options").remove_class("running")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while Downloading: {e.stderr}")
        self.notify(
            f"An error occurred while Downloading: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="Download",
        )
    finally:
        self.query_one("#Download_options").remove_class("running")


def run_FastQC(self):
    """
    Initiates the FastQC process for quality control of sequencing data.

    This function retrieves the input file path from the input field, validates it,
    and starts a new thread to handle the FastQC process.

    Args:
        self: The instance of the class calling this function.
    """
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "FastQC_Input")
    if not input_path:
        input_path = f"{self.workingDir}/data/reads/*"
    self.notify(f"FastQC {str(input_path)}...", title="FastQC")
    logging.info(f"FastQC {str(input_path)}...")
    self.query_one("#FastQC_Horizontal").add_class("running")
    threading.Thread(
        target=_run_FastQC,
        args=(
            self,
            input_path,
        ),
    ).start()


def _run_FastQC(self, input_path: str) -> None:
    """
    Executes the FastQC command using subprocess.

    This function runs the FastQC command to perform quality control on the input file
    and saves the results to the specified output directory.

    Args:
        self: The instance of the class calling this function.
        input_path (str): The path to the input file for FastQC.
    """
    try:
        if not self.workingDir:
            self.workingDir = __file__
        FastQC_cmd = f"fastqc {str(input_path)} -o {self.workingDir}/results/fastqc/"
        self.notify(FastQC_cmd, title="FastQC")
        logging.info("Running command: " + FastQC_cmd)
        subprocess.run(
            [FastQC_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"FastQC completed for {str(input_path)}")
        self.notify(f"FastQC completed for {str(input_path)}", title="FastQC")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during FastQC: {e.stderr}")
        self.notify(
            f"An error occurred during FastQc: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="FastQC",
        )

    finally:
        self.query_one("#FastQC_Horizontal").remove_class("running")


def run_MultiQC(self):
    """
    Initiates the MultiQC process for aggregating results from multiple tools.

    This function retrieves the input file path from the input field, validates it,
    and starts a new thread to handle the MultiQC process.

    Args:
        self: The instance of the class calling this function.
    """
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "MultiQC_Input")
    if not input_path:
        input_path = f"{self.workingDir}/results/fastqc/"

    self.notify(f"MultiQC {str(input_path)}...", title="MultiQC")
    logging.info(f"MultiQC {str(input_path)}...")
    self.query_one("#MultiQC_Horizontal").add_class("running")
    threading.Thread(
        target=_run_MultiQC,
        args=(
            self,
            input_path,
        ),
    ).start()


def _run_MultiQC(self, input_path: str) -> None:
    """
    Executes the MultiQC command using subprocess.

    This function runs the MultiQC command to aggregate results from multiple tools
    and saves the results to the specified output directory.

    Args:
        self: The instance of the class calling this function.
        input_path (str): The path to the input directory for MultiQC.
    """
    try:
        if not self.workingDir:
            self.workingDir = __file__
        MultiQC_cmd = f"multiqc {str(input_path)} -o {self.workingDir}/results/multiqc/"
        self.notify(MultiQC_cmd, title="MultiQC")
        logging.info("Running command: " + MultiQC_cmd)
        subprocess.run(
            [MultiQC_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"MultiQC completed for {str(input_path)}")
        self.notify(f"MultiQC completed for {str(input_path)}", title="MultiQC")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during MultiQC: {e.stderr}")
        self.notify(
            f"An error occurred during MultiQC: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="MultiQC",
        )

    finally:
        self.query_one("#MultiQC_Horizontal").remove_class("running")


def run_bwa_index(self):
    """
    Initiates the BWA indexing process for the reference genome.

    This function retrieves the input file path from the input field, validates it,
    and starts a new thread to handle the BWA indexing process.

    Args:
        self: The instance of the class calling this function.
    """
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "bwa_ref_Input")
    if not input_path:
        input_path = f"{self.workingDir}/data/reference/*"
    self.notify(f"bwa indexing {str(input_path)}...", title="bwa index")
    logging.info(f"bwa indexing {str(input_path)}...")
    self.query_one("#bwa_Horizontal").add_class("running")
    threading.Thread(
        target=_run_bwa_index,
        args=(
            self,
            input_path,
        ),
    ).start()


def _run_bwa_index(self, input_path: str) -> None:
    """
    Executes the BWA index command using subprocess.

    This function runs the BWA index command to index the reference genome
    and saves the results to the specified output directory.

    Args:
        self: The instance of the class calling this function.
        input_path (str): The path to the reference genome file for BWA indexing.
    """
    try:
        if not self.workingDir:
            self.workingDir = __file__
        bwa_index_cmd = f"bwa index {str(input_path)}"
        self.notify(bwa_index_cmd, title="bwa Index")
        logging.info("Running command: " + bwa_index_cmd)
        subprocess.run(
            [bwa_index_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"bwa Index completed for {str(input_path)}")
        self.notify(f"bwa Index completed for {str(input_path)}", title="bwa Index")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during bwa Index: {e.stderr}")
        self.notify(
            f"An error occurred during bwa Index: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="bwa Index",
        )

    finally:
        self.query_one("#bwa_Horizontal").remove_class("running")


def run_bwa_mem(self):
    """
    Initiates the BWA MEM process for aligning reads to the reference genome.

    This function retrieves the input file paths and number of threads from the input fields,
    validates them, and starts a new thread to handle the BWA MEM process.

    Args:
        self: The instance of the class calling this function.
    """
    if not self.workingDir:
        self.workingDir = __file__
    ref_path = get_input(self, "bwa_ref_Input")
    read_path_1 = get_input(self, "bwa_reads_Input_1")
    read_path_2 = get_input(self, "bwa_reads_Input_2")
    no_of_threads = get_input(self, "bwa_threads_Input")
    if not (read_path_1 and read_path_2):
        self.notify("Enter reads file path", title="bwa mem")
    if not ref_path:
        self.notify("Enter reference genome file path", title="bwa mem")
    if not no_of_threads:
        no_of_threads = 4
    self.notify(
        f"aligning reads {str(read_path_1)} {str(read_path_2)} to {str(ref_path)} using bwa mem",
        title="bwa mem",
    )
    logging.info(
        f"aligning reads {str(read_path_1)} {str(read_path_2)} to {str(ref_path)} using bwa mem"
    )
    self.query_one("#bwa_Horizontal").add_class("running")
    threading.Thread(
        target=_run_bwa_mem,
        args=(
            self,
            ref_path,
            read_path_1,
            read_path_2,
            no_of_threads,
        ),
    ).start()


def _run_bwa_mem(
    self, ref_path: str, read_path_1: str, read_path_2: str, no_of_threads: str
) -> None:
    """
    Executes the BWA MEM command using subprocess.

    This function runs the BWA MEM command to align the reads to the reference genome
    and saves the results to the specified output directory.

    Args:
        self: The instance of the class calling this function.
        ref_path (str): The path to the reference genome file for BWA MEM.
        read_path (str): The path to the reads file for BWA MEM.
        no_of_threads (str): The number of threads to use for BWA MEM.
    """
    try:
        if not self.workingDir:
            self.workingDir = __file__
        bwa_mem_cmd = f"bwa mem -t {str(no_of_threads)} {str(ref_path)} {str(read_path_1)} {str(read_path_2)} -o {self.workingDir}/results/sam/aligned.sam"
        self.notify(bwa_mem_cmd, title="bwa mem")
        logging.info("Running command: " + bwa_mem_cmd)
        subprocess.run(
            [bwa_mem_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(
            f"Completed aligning reads {str(read_path_1)} {str(read_path_2)} to {str(ref_path)} using bwa mem"
        )
        self.notify(
            f"Completed aligning reads {str(read_path_1)} {str(read_path_2)} to {str(ref_path)} using bwa mem",
            title="bwa mem",
        )

        flagstat_bam(
            f"{self.workingDir}/results/sam/aligned.sam",
            f"{self.workingDir}/results/sam/aligned.sam.flagstat",
        )

        MultiQC_cmd = f"multiqc {self.workingDir}/results/sam/aligned.sam.flagstat -o {self.workingDir}/results/multiqc/sam"
        subprocess.run(
            [MultiQC_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Error aligning using bwa mem: {e.stderr}")
        self.notify(
            f"Error aligning using bwa mem: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="bwa mem",
        )

    finally:
        self.query_one("#bwa_Horizontal").remove_class("running")


# def run_Trimmomatic(self):
#     input_path = self.query_one("#bwa_ref_Input", Input)
#     input_path = input_path.value.strip()
#     if not input_path:
#         input_path = f"{self.workingDir}/data/reference/*"
#     self.notify(f"bwa indexing {str(input_path)}...", title="bwa index")
#     logging.info(f"bwa indexing {str(input_path)}...")
#     self.query_one("#bwa_Horizontal").add_class("running")
#     threading.Thread(
#         target=_run_Trimmomatic,
#         args=(
#             self,
#             input_path,
#         ),
#     ).start()
#
#
# def _run_Trimmomatic(self, input_path: str) -> None:
#     try:
#         Trimmomatic_cmd = (f"trimmomatic SE {input_fastq} {output_fastq} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
#         bwa_index_cmd = f"bwa index {str(input_path)}"
#         self.notify(bwa_index_cmd, title="bwa Index")
#         logging.info("Running command: " + bwa_index_cmd)
#         subprocess.run(
#             [Trimmomatic_cmd],
#             check=True,
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE,
#             text=True,
#             shell=True,
#         )
#         logging.info(f"bwa Index completed for {str(input_path)}")
#         self.notify(f"bwa Index completed for {str(input_path)}", title="bwa Index")
#         self.query_one("#bwa_Horizontal").remove_class("running")
#     except subprocess.CalledProcessError as e:
#         logging.error(f"An error occurred during bwa Index: {e.stderr}")
#         self.notify(
#             f"An error occurred during bwa Index: {e.stderr}",
#             severity="error", timeout=10.0,
#             title="bwa Index",
#         )
