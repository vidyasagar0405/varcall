import logging
import subprocess
from os import listdir
from webbrowser import open_new_tab

from textual import on
from textual.app import App, ComposeResult
from textual.containers import ScrollableContainer
from textual.reactive import reactive
from textual.widgets import (
    Button,
    Footer,
    Header,
    Input,
    Label,
    Select,
    TabbedContent,
    TabPane,
)

# Importing custom modules and functions
from modules.Bcftools_widgets import BcftoolsWidgets
from modules.exec_func.hometab_funcs import *  # noqa: F403
from modules.Help import HelpMarkdown
from modules.Home_widgets import HomeWidgets
from modules.logging import setup_logging
from modules.Pipeline_widgets import PipelineWidgets
from modules.Samtools_widgets import SamtoolsWidgets
from modules.YesOrNo import YesOrNo

# Setup logging configuration
setup_logging()


class Varcall(App[None]):
    # Key bindings for the application
    BINDINGS = [
        ("f1", "show_help", "Help"),
        ("q", "exit_app", "Exit App"),
        ("D", "toggle_dark_mode", "Dark/Light mode"),
        ("d", "run_download", "Download"),
        ("f", "run_fastqc", "Run FastQC"),
        ("m", "run_multiqc", "Run MultiQC"),
        ("i", "run_index", "Index reference genome"),
        ("a", "run_alignment", "Align reads"),
        ("d", "run_download", "Download"),
    ]
    CSS_PATH = "main.css"

    def on_mount(self) -> None:
        self.action_show_help()

    # Reactive variables to track state
    workingDir = reactive("")
    outputfile_name = reactive("")
    outputdir_name = reactive("")
    full_output_path = reactive("")

    def compose(self) -> ComposeResult:
        """
        Composes the UI layout of the application.

        This method sets up the main components of the UI, including the header,
        footer, and tabbed content for different functionalities.
        """
        with ScrollableContainer(id="ScrollableContainer"):
            yield Header(show_clock=True)
            yield Footer(show_command_palette=True)
            with TabbedContent():
                with TabPane("Home", id="HomeTab"):
                    yield HomeWidgets()
                with TabPane("Samtools", id="SamtoolsTab"):
                    yield SamtoolsWidgets()
                with TabPane("Bcftools", id="BcftoolsTab"):
                    yield BcftoolsWidgets()
                with TabPane("Pipeline", id="Pipeline"):
                    yield PipelineWidgets()
                with TabPane("Help", id="HelpTab"):
                    yield HelpMarkdown()

    # Event handler for project name input submission
    @on(Input.Submitted, "#Input_Project_name")
    def update_working_dir(self, event: Input.Submitted) -> None:
        """
        Updates the working directory based on the project name input.

        This method is triggered when the project name input is submitted.
        It updates the working directory, creates necessary subdirectories,
        and notifies the user.

        Args:
            event (Input.Submitted): The event object containing the input value.
        """
        self.workingDir = str(event.input.value)
        self.update_full_output_path()
        self.notify(
            f"New Project Started and working directory updated to {self.workingDir}"
        )
        # Create necessary directories for the project
        makedir = f"mkdir -p {self.workingDir} {self.workingDir}/results {self.workingDir}/data {self.workingDir}/data/reads {self.workingDir}/data/reference {self.workingDir}/results/fastqc {self.workingDir}/results/multiqc {self.workingDir}/results/multiqc/sam {self.workingDir}/results/sam {self.workingDir}/results/bam {self.workingDir}/results/vcf {self.workingDir}/results/bcf"
        self.notify(makedir)
        try:
            subprocess.run(
                makedir.split(),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred: {e.stderr}")
            self.notify(f"An error occurred: {e.stderr}", severity="error")

    # Event handler for output file name input change
    @on(Input.Changed, "#Input_outputfile_name")
    def update_outputfile_name(self, event: Input.Changed) -> None:
        """
        Updates the output file name based on the input change.

        This method is triggered when the output file name input is changed.

        Args:
            event (Input.Changed): The event object containing the input value.
        """
        self.outputfile_name = event.input.value
        self.update_full_output_path()

    # Event handler for output directory selection change
    @on(Select.Changed, "#Select_outputdir")
    def update_outputdir_name(self, event: Select.Changed) -> None:
        """
        Updates the output directory name based on the selection change.

        This method is triggered when the output directory selection is changed.
        It updates the output directory name and the full output path.

        Args:
            event (Select.Changed): The event object containing the selected value.
        """
        self.outputdir_name = event.select.value
        self.update_full_output_path()

    # Event handler for download button press
    @on(Button.Pressed, "#Download_button")
    def call_run_download(self) -> None:
        """
        Initiates the download process when the download button is pressed.

        This method is triggered when the download button is pressed.
        It calls the run_download function to start the download process.
        """
        run_download(self)  # noqa: F405

    # Event handler for FastQC button press
    @on(Button.Pressed, "#FastQC_Button")
    def call_run_FastQC(self) -> None:
        """
        Initiates the FastQC process when the FastQC button is pressed.

        This method is triggered when the FastQC button is pressed.
        It calls the run_FastQC function to start the FastQC process.
        """
        run_FastQC(self)  # noqa: F405

    # Event handler for MultiQC button press
    @on(Button.Pressed, "#MultiQC_Button")
    def call_run_MultiQC(self) -> None:
        """
        Initiates the MultiQC process when the MultiQC button is pressed.

        This method is triggered when the MultiQC button is pressed.
        It calls the run_MultiQC function to start the MultiQC process.
        """
        run_MultiQC(self)  # noqa: F405

    # Event handler for BWA index button press
    @on(Button.Pressed, "#bwa_index_Button")
    def call_run_bwa_index(self) -> None:
        """
        Initiates the BWA indexing process when the BWA index button is pressed.

        This method is triggered when the BWA index button is pressed.
        It calls the run_bwa_index function to start the BWA indexing process.
        """
        run_bwa_index(self)  # noqa: F405

    # Event handler for BWA align button press
    @on(Button.Pressed, "#bwa_align_Button")
    def call_run_bwa_mem(self) -> None:
        """
        Initiates the BWA MEM process when the BWA align button is pressed.

        This method is triggered when the BWA align button is pressed.
        It calls the run_bwa_mem function to start the BWA MEM process.
        """
        run_bwa_mem(self)  # noqa: F405

    # Event handler to view FastQC results
    @on(Button.Pressed, "#view_FastQC_res")
    def view_FastQC_res(self) -> None:
        """
        Opens the FastQC results in a web browser.

        This method is triggered when the view FastQC results button is pressed.
        It opens the FastQC results HTML files present
        in the workingDir/results/fastqc directory in a new browser tab.
        """
        dir_path = f"{self.workingDir}/results/fastqc/"
        for file in listdir(dir_path):
            if file.endswith(".html"):
                open_new_tab(f"{dir_path}/{file}")

    # Event handler to view MultiQC results
    @on(Button.Pressed, "#view_MultiQC_res")
    def view_MultiQC_res(self) -> None:
        """
        Opens the MultiQC results in a web browser.

        This method is triggered when the view MultiQC results button is pressed.
        It opens the MultiQC results HTML files present
        in the workingDir/results/multiqc directory in a new browser tab.
        """
        dir_path = f"{self.workingDir}/results/multiqc/"
        for file in listdir(dir_path):
            if file.endswith(".html"):
                open_new_tab(f"{dir_path}/{file}")

    @on(Button.Pressed, "#view_bwa_res")
    def view_bwa_res(self) -> None:
        """
        This method is triggered when the view bwa results button is pressed.
        Executes flagstat_bam() and multiqc on the result file and then,
        Opens the it in a web browser.
        """

        dir_path = f"{self.workingDir}/results/multiqc/sam/"
        for file in listdir(dir_path):
            if file.endswith(".html"):
                open_new_tab(f"{dir_path}/{file}")

    # Handels the key press for Toggle dark mode
    def action_toggle_dark_mode(self) -> None:
        self.dark = not self.dark

    # Handels the key press for Run download action
    def action_run_download(self) -> None:
        self.call_run_download()

    # Handels the key press for Run FastQC action
    def action_run_fastqc(self) -> None:
        self.call_run_FastQC()

    # Handels the key press for Run MultiQC action
    def action_run_multiqc(self) -> None:
        self.call_run_MultiQC()

    # Handels the key press for Run alignment action
    def action_run_alignment(self) -> None:
        self.call_run_bwa_mem()

    # Handels the key press for Run index action
    def action_run_index(self) -> None:
        self.call_run_bwa_index()

    # Handels the key press for Show help action
    def action_show_help(self) -> None:
        self.query_one(TabbedContent).active = "SamtoolsTab"

    # Handels the key press for Exit application action
    def action_exit_app(self) -> None:
        self.push_screen(
            YesOrNo("Do you want to exit application?"), self.maybe_exit_app
        )

    # Confirm exit application
    def maybe_exit_app(self, bool) -> None:
        if bool:
            self.exit()

    # Update the full output path based on current state
    def update_full_output_path(self) -> None:
        self.full_output_path = (
            f"{self.workingDir}/data/{self.outputdir_name}/{self.outputfile_name}"
        )

    # Watcher for full output path changes
    def watch_full_output_path(self, full_output_path: str) -> None:
        self.query_one("#output_path_label", Label).update(
            f"Saved in: {full_output_path}"
        )


# Main entry point of the application
if __name__ == "__main__":
    logging.info("Application started.")
    Varcall().run()
    logging.info("Application closed.")
