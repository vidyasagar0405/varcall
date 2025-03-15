from pathlib import Path

from textual.app import App, ComposeResult
from textual.containers import ScrollableContainer
from textual.widgets import (
    Footer,
    Header,
    Markdown,
    MarkdownViewer,
    TabbedContent,
    TabPane,
)

from varcall.components.input_block import ProcessWidgets
from varcall.process.process_dict import BCFTOOLS_PROCESS, GATK_PROCESS, HOME_PROCESSES, PIPLELINES, SAMTOOLS_PROCESSES
from varcall.help import help_path

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
    ]

    CSS_PATH = "./app.css"

    def compose(self) -> ComposeResult:
        with ScrollableContainer(id="ScrollableContainer"):
            yield Header()
            yield Footer(show_command_palette=True)
            with TabbedContent():
                with TabPane("Home", id="HomeTab"):
                    yield ProcessWidgets(HOME_PROCESSES)
                with TabPane("Samtools", id="SamtoolsTab"):
                    yield ProcessWidgets(SAMTOOLS_PROCESSES)
                with TabPane("Bcftools", id="BcftoolsTab"):
                    yield ProcessWidgets(BCFTOOLS_PROCESS)
                with TabPane("GATK", id="GATKTab"):
                    yield ProcessWidgets(GATK_PROCESS)
                with TabPane("Pipeline", id="Pipeline"):
                    yield ProcessWidgets(PIPLELINES)
                with TabPane("Help", id="HelpTab"):
                    yield MarkdownViewer(Path(help_path).read_text())



# Main entry point of the application
def main():
    Varcall().run()


if __name__ == "__main__":
    main()
