import subprocess
import threading
import glob
import logging
from pathlib import Path
from textual import on
from textual.app import App, ComposeResult
from textual.containers import HorizontalScroll, ScrollableContainer, Horizontal
from textual.widgets import Button, Footer, Header, TabbedContent, TabPane, Input, LoadingIndicator, Label
from modules.widgets import InputBlock, InputProjectName
from modules.logging_util import setup_logging
from modules.Help import HelpMarkdown
from textual_fspicker import SelectDirectory

setup_logging()

class Varcall(App):
    BINDINGS = [
        ("d", "toggle_dark_mode", "Dark/Light mode"),
        ("f", "run_fastqc", "Run FastQC")
    ]
    CSS_PATH = "Varcall.css"

    def compose(self) -> ComposeResult:
        with ScrollableContainer(id="ScrollableContainer"):
            yield Header(show_clock=True)
            yield Footer(show_command_palette=True)
            with TabbedContent():
                with TabPane("Home", id="HomeTab"):
                    yield InputProjectName()
                    yield InputBlock("Select Reference Genome file", "ReferenceGenome")
                    yield InputBlock("Select Reads (fastq files)", "Reads")
                    with HorizontalScroll(id="HomeTab_FastQC"):
                        yield Label("FastQC description", id="FastQC_description")
                        yield Button("FastQC", id="FastQCButton")
                        yield Button("View Results", id="FastQCResultsButton")
                        yield LoadingIndicator(id="FastQCLoadingIndicator")
                    with HorizontalScroll(id="HomeTab_MultiQC"):
                        yield Label("MultiQC description", id="MultiQC_description")
                        yield Button("MultiQC", id="MultiQCButton")
                        yield Button("View Results", id="MultiQCResultsButton")
                        yield LoadingIndicator(id="MultiQCLoadingIndicator")
                    with HorizontalScroll(id="HomeTab_Bwa"):
                        yield Label("Bwa description", id="Bwa_description")
                        yield Button("Index Ref", id="BwaIndexButton")
                        yield Button("Alignment", id="AlignmentButton")
                        yield LoadingIndicator(id="BwaLoadingIndicator")
                with TabPane("Samtools", id="SamtoolsTab"):
                    with ScrollableContainer(id="SamtoolsTabScrollableContainer"):
                        yield HelpMarkdown()
    
                with TabPane("Bcftools", id="BcftoolsTab"):
                    with ScrollableContainer(id="BcftoolsTabScrollableContainer"):
                        yield HelpMarkdown()
    
                with TabPane("Help", id="HelpTab"):
                    with ScrollableContainer(id="HelpTabScrollableContainer"):
                        yield HelpMarkdown()
    
    def  show_selected(self, to_show: Path | None) -> None:
        self.query_one("#Label_selected_Reads_directory").update("Cancelled" if to_show is None else str("Selected directory: " + str(to_show)))

    @on(Button.Pressed, "#Select_Directory_Button")
    def select_directory(self) -> None:
        self.push_screen(
            SelectDirectory("."),
            callback=self.show_selected,
            )
        self.notify("Not fully implemented yet", severity="warning")
        logging.warning("Select directory function not fully implemented yet.")

    @on(Button.Pressed, "#download_Button")
    def run_download_url(self) -> None:
        input_file_name_widget = self.query_one("#Input_Reads_download_url", Input)
        download_url = input_file_name_widget.value.strip()
        # download_url = self.expand_wildcard(download_url)

        if not download_url:
            self.notify("Please provide a valid file path", severity="warning")
            return
        else:
            self.notify(f"Downloading {' '.join(download_url)}...")
            logging.info(f"Downloading {' '.join(download_url)}...")
            threading.Thread(target=self._run_download_url, args=(download_url,)).start()

    def _run_download_url(self, download_url: str) -> None:
        try:
            logging.info(f"Running command: " + "cdownload_url " + "-L " +''.join(download_url))
            result = subprocess.run(
                ["cdownload_url", "-L", ''.join(download_url), "-o", ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logging.info(f"Downloaded {' '.join(download_url)}")
            self.notify(f"Downloaded {' '.join(download_url)}")
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while Downloading: {e.stderr}")
            self.notify(f"An error occurred while Downloading: {e.stderr}", severity="error")


    @on(Button.Pressed, "#FastQCButton")
    def run_fastqc(self) -> None:
        input_file_name_widget = self.query_one("#Input_Reads_file_name", Input)
        path_to_reads = input_file_name_widget.value.strip()
        path_to_reads = self.expand_wildcard(path_to_reads)

        if not path_to_reads:
            self.notify("Please provide a valid file path", severity="warning")
            return
        else:
            self.query_one("#HomeTab_FastQC").add_class("running")
            self.notify(f"Starting FastQC analysis for {' '.join(path_to_reads)}...")
            logging.info(f"Starting FastQC analysis for {' '.join(path_to_reads)}...")
            threading.Thread(target=self._run_fastqc, args=(path_to_reads,)).start()

    def expand_wildcard(self, path: str):
        if "*" in path:
            matched_files = glob.glob(path)
            if not matched_files:
                self.notify("No files found matching the pattern", severity="warning")
                logging.warning("No files found for the pattern: " + path)
                return []
            return [str(file) for file in matched_files]
        return [path]

    def _run_fastqc(self, path_to_reads: str) -> None:
        try:
            logging.info(f"Running command: fastqc {' '.join(path_to_reads)}")
            result = subprocess.run(
                ["fastqc"] + path_to_reads,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logging.info(f"FastQC analysis completed for {' '.join(path_to_reads)}")
            self.notify(f"FastQC analysis completed for {' '.join(path_to_reads)}")
            self.query_one("#HomeTab_FastQC").remove_class("running")
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while running FastQC: {e.stderr}")
            self.notify(f"An error occurred while running FastQC: {e.stderr}", severity="error")
            self.query_one("#HomeTab_FastQC").remove_class("running")

    @on(Button.Pressed, "#MultiQCButton")
    def run_multiqc(self) -> None:
        input_file_name_widget = self.query_one("#Input_Reads_file_name", Input)
        path_to_reads = input_file_name_widget.value.strip()
        path_to_reads = self.expand_wildcard(path_to_reads)

        if not path_to_reads:
            self.notify("Please provide a valid file path", severity="warning")
            return
        self.query_one("#HomeTab_MultiQC").add_class("running")
        self.notify(f"Starting MultiQC analysis for {' '.join(path_to_reads)}...")
        logging.info(f"Starting MultiQC analysis for {' '.join(path_to_reads)}...")
        threading.Thread(target=self._run_multiqc, args=(path_to_reads,)).start()

    def _run_multiqc(self, path_to_reads: str) -> None:
        try:
            logging.info(f"Running command: multiqc {' '.join(path_to_reads)}")
            result = subprocess.run(
                ["multiqc"] + path_to_reads,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logging.info(f"MultiQC analysis completed for {' '.join(path_to_reads)}")
            self.notify(f"MultiQC analysis completed for {' '.join(path_to_reads)}")
            self.query_one("#HomeTab_MultiQC").remove_class("running")
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while running MultiQC: {e.stderr}")
            self.notify(f"An error occurred while running MultiQC: {e.stderr}", severity="error")
            self.query_one("#HomeTab_MultiQC").remove_class("running")

    @on(Button.Pressed, "#BwaIndexButton")
    def run_bwa_index(self) -> None:
        input_file_name_widget = self.query_one("#Input_ReferenceGenome_file_name", Input)
        reference_file = input_file_name_widget.value.strip()

        if not reference_file:
            self.notify("Please provide a valid reference file path", severity="warning")
            return

        self.notify(f"Starting BWA index for {reference_file}...")
        logging.info(f"Starting BWA index for {reference_file}...")
        threading.Thread(target=self._run_bwa_index, args=(reference_file,)).start()

    def _run_bwa_index(self, reference_file: str) -> None:
        try:
            logging.info(f"Running command: bwa index {reference_file}")
            result = subprocess.run(
                ["bwa", "index", reference_file],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logging.info(f"BWA index completed for {reference_file}")
            self.notify(f"BWA index completed for {reference_file}")
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while running BWA index: {e.stderr}")
            self.notify(f"An error occurred while running BWA index: {e.stderr}", severity="error")

    @on(Button.Pressed, "#AlignmentButton")
    def run_alignment(self) -> None:
        input_file_name_widget = self.query_one("#Input_Reads_file_name", Input)
        reference_file_widget = self.query_one("#Input_ReferenceGenome_file_name", Input)
        path_to_reads = input_file_name_widget.value.strip().split()  # Split into separate files

        # Ensure there are exactly two read files
        if len(path_to_reads) != 2 or not reference_file_widget.value.strip():
            self.notify("Please provide valid file paths for both reads and a reference", severity="warning")
            return

        reference_file = reference_file_widget.value.strip()

        self.query_one("#HomeTab_Bwa").add_class("running")
        self.notify(f"Starting alignment for {path_to_reads[0]} and {path_to_reads[1]} against {reference_file}...")
        logging.info(f"Starting alignment for {path_to_reads[0]} and {path_to_reads[1]} against {reference_file}...")
        threading.Thread(target=self._run_alignment, args=(path_to_reads[0], path_to_reads[1], reference_file)).start()

    def _run_alignment(self, read_file_1: str, read_file_2: str, reference_file: str) -> None:
        # Check if the reference file exists
        if not Path(reference_file).exists():
            logging.error(f"Reference file does not exist: {reference_file}")
            self.notify(f"Reference file does not exist: {reference_file}", severity="error")
            return

        # Check if the read files exist
        for read_file in [read_file_1, read_file_2]:
            if not Path(read_file).exists():
                logging.error(f"Read file does not exist: {read_file}")
                self.notify(f"Read file does not exist: {read_file}", severity="error")
                return

        try:
            logging.info(f"Running command: bwa mem {reference_file} {read_file_1} {read_file_2}")
            result = subprocess.run(
                ["bwa", "mem", reference_file, read_file_1, read_file_2, "-o", "aligned.sam"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            self.query_one("#HomeTab_Bwa").remove_class("running")
            logging.info(f"Alignment completed for {read_file_1} and {read_file_2} against {reference_file}")
            self.notify(f"Alignment completed for {read_file_1} and {read_file_2} against {reference_file}")
        except subprocess.CalledProcessError as e:
            self.query_one("#HomeTab_Bwa").remove_class("running")
            logging.error(f"An error occurred while running alignment: {e.stderr}")
            self.notify(f"An error occurred while running alignment: {e.stderr}", severity="error")

    @on(Button.Pressed, "#FastQCResultsButton")
    def open_fastqc_results(self):
        try:
            result = subprocess.run(
                    ["firefox samples/reads/*fastqc.html"],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=True
                )
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred: {e.stderr}")
            self.notify(f"An error occurred: {e.stderr}", severity="error")

    @on(Button.Pressed, "#MultiQCResultsButton")
    def open_multiqc_results(self):
        try:
            result = subprocess.run(
                    ["./open_qc.sh"],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=True
                )
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred: {e.stderr}")
            self.notify(f"An error occurred: {e.stderr}", severity="error")

    def action_toggle_dark_mode(self):
        self.dark = not self.dark

    def action_run_fastqc(self):
        self.run_fastqc()

if __name__ == "__main__":
    logging.info("Application started.")
    Varcall().run()
    logging.info("Application closed.")
