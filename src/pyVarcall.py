import subprocess
# import threading
# import glob
import logging
# from pathlib import Path
from webbrowser import open_new_tab
from os import listdir
from textual import on
from textual.app import App, ComposeResult
from textual.containers import ScrollableContainer, Horizontal, Container
from textual.widgets import Button, Footer, Header, Select, TabbedContent, TabPane, Input, Label
from textual.screen import ModalScreen
from modules.logging import setup_logging
from modules.Help import HelpMarkdown
from modules.Home_widgets import HomeWidgets
from modules.exec_func import run_download, run_FastQC, run_MultiQC
from textual.reactive import reactive

setup_logging()

class YesOrNo(ModalScreen):
    def __init__(self, Title: str="Input") -> None:
        self.Title = Title
        super().__init__()

    def compose(self) -> ComposeResult:
        with Container(id="yesorno"):
            yield Label(self.Title)
            with Horizontal(id="horizontal_yesorno"):
                yield Button.success("Yes", id="yes")
                yield Button.error("No", id="no")


    @on(Button.Pressed)
    def close_modalscreen(self, event: Button.Pressed) -> None:
        self.dismiss(event.button.id == "yes")

class Varcall(App[None]):

    BINDINGS = [
        ("f1", "show_help", "Help"),
        ("ctrl+c", "exit_app" , "Exit App"),
        ("D", "toggle_dark_mode", "Dark/Light mode"),
        ("d", "run_download", "Download"),
        ("f", "run_fastqc", "Run FastQC"),
        ("m", "run_multiqc", "Run MultiQC"),
        ("i", "run_index", "Index reference genome"),
        ("a", "run_alignment", "Align reads"),
        ("d", "run_download", "Download"),
    ]
    CSS_PATH = "pyVarcall.css"


    workingDir = reactive("")
    outputfile_name = reactive("")
    outputdir_name = reactive("")
    full_output_path = reactive("")

    # NOTE: disabled during development for comvenience
    # def on_mount(self) -> None:
    #     self.action_show_help()


    def compose(self) -> ComposeResult:
        with ScrollableContainer(id="ScrollableContainer"):
            yield Header(show_clock=True)
            yield Footer(show_command_palette=True)
            with TabbedContent():
                with TabPane("Home", id="HomeTab"):
                    yield HomeWidgets()
                with TabPane("Help", id="HelpTab"):
                    yield HelpMarkdown()

    @on(Input.Submitted, "#Input_Project_name")
    def update_working_dir(self, event: Input.Submitted) -> None:
        self.workingDir = str(event.input.value)
        self.update_full_output_path()
        self.notify(f"New Project Started and working directory updated to {self.workingDir}")
        makedir = f"mkdir -p {self.workingDir} {self.workingDir}/results {self.workingDir}/data {self.workingDir}/data/reads {self.workingDir}/data/reference {self.workingDir}/results/fastqc {self.workingDir}/results/multiqc {self.workingDir}/results/sam {self.workingDir}/results/bam {self.workingDir}/results/vcf {self.workingDir}/results/bcf"
        self.notify(makedir)
        try:
            result = subprocess.run(
                    [makedir],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=True
                )
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred: {e.stderr}")
            self.notify(f"An error occurred: {e.stderr}", severity="error")

    @on(Input.Changed, "#Input_outputfile_name")
    def update_outputfile_name (self, event: Input.Changed) -> None:
        self.outputfile_name = event.input.value
        self.update_full_output_path()

    @on(Select.Changed, "#Select_outputdir")
    def update_outputdir_name (self, event: Select.Changed) -> None:
        self.outputdir_name = event.select.value
        self.update_full_output_path()

    @on(Button.Pressed, "#Download_button")
    def call_run_download(self) -> None:
        run_download(self)

    @on(Button.Pressed, "#FastQC_Button")
    def call_run_FastQC(self) -> None:
        run_FastQC(self)

    @on(Button.Pressed, "#MultiQC_Button")
    def call_run_MultiQC(self) -> None:
        run_MultiQC(self)


    @on(Button.Pressed, "#view_FastQC_res")
    def view_FastQC_res(self) -> None:
        dir_path = f"{self.workingDir}/results/fastqc/"
        for file in listdir(dir_path):
            if file.endswith(".html"):
                open_new_tab(f"{dir_path}/{file}")

    @on(Button.Pressed, "#view_MultiQC_res")
    def view_MultiQC_res(self) -> None:
        dir_path = f"{self.workingDir}/results/multiqc/"
        for file in listdir(dir_path):
            if file.endswith(".html"):
                open_new_tab(f"{dir_path}/{file}")


    def action_toggle_dark_mode(self) -> None:
        self.dark = not self.dark

    def action_run_download(self) -> None:
        self.call_run_download()

    def action_run_fastqc(self) -> None:
        self.call_run_FastQC()

    def action_run_multiqc(self) -> None:
        self.call_run_MultiQC()

    def action_run_alignment(self) -> None:
        self.call_run_FastQC()

    def action_run_index(self) -> None:
        self.call_run_FastQC()

    def action_show_help(self) -> None:
        self.query_one(TabbedContent).active = "HelpTab"

    def action_exit_app(self) -> None:
        self.push_screen(YesOrNo("Do you want to exit application?"), self.maybe_exit_app)

    def maybe_exit_app (self, bool) -> None:
        if bool:
            self.exit()

    def update_full_output_path(self) -> None:
        self.full_output_path = f"{self.workingDir}/data/{self.outputdir_name}/{self.outputfile_name}"

    def watch_full_output_path(self, full_output_path: str) -> None:
        self.query_one("#output_path_label", Label).update(f"Saved in: {full_output_path}")

if __name__ == "__main__":
    logging.info("Application started.")
    Varcall().run()
    logging.info("Application closed.")
