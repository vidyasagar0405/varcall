from textual.app import App, ComposeResult
from textual.containers import ScrollableContainer
from textual.widgets import (
    Button,
    Footer,
    Header,
    TabbedContent,
    TabPane,
)

from varcall.components.input_block import ProcessWidgets

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

    CSS = """
#--content-tab-HomeTab{color: green;}
#--content-tab-HelpTab{color: #1691CE;}
    """

    def compose(self) -> ComposeResult:
        with ScrollableContainer(id="ScrollableContainer"):
            yield Header()
            yield Footer(show_command_palette=True)
            with TabbedContent():
                with TabPane("Home", id="HomeTab"):
                    yield ProcessWidgets()
                with TabPane("Samtools", id="SamtoolsTab"):
                    yield Button("Samtools")
                with TabPane("Bcftools", id="BcftoolsTab"):
                    yield Button("Bcftools")
                with TabPane("Pipeline", id="Pipeline"):
                    yield Button("Pipeline")
                with TabPane("Help", id="HelpTab"):
                    yield Button("Help")



# Main entry point of the application
def main():
    Varcall().run()


if __name__ == "__main__":
    main()
