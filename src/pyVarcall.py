import logging
import subprocess
from os import listdir
from webbrowser import open_new_tab

# Importing custom modules and functions
from modules.Bcftools_widgets import BcftoolsWidgets
from modules.Samtools_widgets import SamtoolsWidgets
from modules.exec_func.hometab_funcs import (
    run_download,
    run_FastQC,
    run_MultiQC,
    run_bwa_index,
    run_bwa_mem,
)
from modules.Help import HelpMarkdown
from modules.Home_widgets import HomeWidgets
from modules.logging import setup_logging
from modules.YesOrNo import YesOrNo
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

# Setup logging configuration
setup_logging()

class Varcall(App[None]):
    # Key bindings for the application
    BINDINGS = [
        ("f1", "show_help", "Help"),
        ("ctrl+c", "exit_app", "Exit App"),
        ("D", "toggle_dark_mode", "Dark/Light mode"),
        ("d", "run_download", "Download"),
        ("f", "run_fastqc", "Run FastQC"),
        ("m", "run_multiqc", "Run MultiQC"),
        ("i", "run_index", "Index reference genome"),
        ("a", "run_alignment", "Align reads"),
        ("d", "run_download", "Download"),
    ]
    CSS_PATH = "pyVarcall.css"

    # Reactive variables to track state
    workingDir = reactive("")
    outputfile_name = reactive("")
    outputdir_name = reactive("")
    full_output_path = reactive("")

    # NOTE: disabled during development for comvenience
    # def on_mount(self) -> None:
    #     self.action_show_help()

    # def on_mount(self) -> None:
    #     self.push_screen(HomeWidgets())

    # Compose the UI layout
    def compose(self) -> ComposeResult:
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
                with TabPane("Help", id="HelpTab"):
                    yield HelpMarkdown()

    # Event handler for project name input submission
    @on(Input.Submitted, "#Input_Project_name")
    def update_working_dir(self, event: Input.Submitted) -> None:
        self.workingDir = str(event.input.value)
        self.update_full_output_path()
        self.notify(
            f"New Project Started and working directory updated to {self.workingDir}"
        )
        # Create necessary directories for the project
        makedir = f"mkdir -p {self.workingDir} {self.workingDir}/results {self.workingDir}/data {self.workingDir}/data/reads {self.workingDir}/data/reference {self.workingDir}/results/fastqc {self.workingDir}/results/multiqc {self.workingDir}/results/sam {self.workingDir}/results/bam {self.workingDir}/results/vcf {self.workingDir}/results/bcf"
        self.notify(makedir)
        try:
            subprocess.run(
                [makedir],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                shell=True,
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred: {e.stderr}")
            self.notify(f"An error occurred: {e.stderr}", severity="error")

    # Event handler for output file name input change
    @on(Input.Changed, "#Input_outputfile_name")
    def update_outputfile_name(self, event: Input.Changed) -> None:
        self.outputfile_name = event.input.value
        self.update_full_output_path()

    # Event handler for output directory selection change
    @on(Select.Changed, "#Select_outputdir")
    def update_outputdir_name(self, event: Select.Changed) -> None:
        self.outputdir_name = event.select.value
        self.update_full_output_path()

    # Event handler for download button press
    @on(Button.Pressed, "#Download_button")
    def call_run_download(self) -> None:
        run_download(self)

    # Event handler for FastQC button press
    @on(Button.Pressed, "#FastQC_Button")
    def call_run_FastQC(self) -> None:
        run_FastQC(self)

    # Event handler for MultiQC button press
    @on(Button.Pressed, "#MultiQC_Button")
    def call_run_MultiQC(self) -> None:
        run_MultiQC(self)

    # Event handler for BWA index button press
    @on(Button.Pressed, "#bwa_index_Button")
    def call_run_bwa_index(self) -> None:
        run_bwa_index(self)

    # Event handler for BWA align button press
    @on(Button.Pressed, "#bwa_align_Button")
    def call_run_bwa_mem(self) -> None:
        run_bwa_mem(self)

    # Event handler to view FastQC results
    @on(Button.Pressed, "#view_FastQC_res")
    def view_FastQC_res(self) -> None:
        dir_path = f"{self.workingDir}/results/fastqc/"
        for file in listdir(dir_path):
            if file.endswith(".html"):
                open_new_tab(f"{dir_path}/{file}")

    # Event handler to view MultiQC results
    @on(Button.Pressed, "#view_MultiQC_res")
    def view_MultiQC_res(self) -> None:
        dir_path = f"{self.workingDir}/results/multiqc/"
        for file in listdir(dir_path):
            if file.endswith(".html"):
                open_new_tab(f"{dir_path}/{file}")

    # Toggle dark mode
    def action_toggle_dark_mode(self) -> None:
        self.dark = not self.dark

    # Run download action
    def action_run_download(self) -> None:
        self.call_run_download()

    # Run FastQC action
    def action_run_fastqc(self) -> None:
        self.call_run_FastQC()

    # Run MultiQC action
    def action_run_multiqc(self) -> None:
        self.call_run_MultiQC()

    # Run alignment action
    def action_run_alignment(self) -> None:
        self.call_run_bwa_mem()

    # Run index action
    def action_run_index(self) -> None:
        self.call_run_bwa_index()

    # Show help action
    def action_show_help(self) -> None:
        self.query_one(TabbedContent).active = "HelpTab"

    # Exit application action
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
