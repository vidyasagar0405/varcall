from pathlib import Path
from textual.app import ComposeResult
from textual.containers import ScrollableContainer
from textual.suggester import Suggester
from textual.widgets import Button, Input, Label, LoadingIndicator, Static

from varcall.process.process_class import ProcessConfig


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


class ProcessWidgets(Static):
    DEFAULT_CSS = """
    ProcessWidgets {
        width: 100%;
        height: auto;
        padding: 1;
    }

    ScrollableContainer {
        width: 100%;
        height: auto;
        margin: 0 1;
    }

    .process-widget {
        layout: grid;
        grid-size: 1;
        grid-rows: 4;
        background: $boost;
        height: auto;
        min-height: 12;
        margin: 1 0;
        padding: 1;
        border: tall $primary;
    }

    .process-title {
        text-style: bold;
        background: $panel;
        color: $text;
        width: 100%;
        height: 3;
        content-align: center middle;
        border: tall $primary-darken-2;
    }

    .process-description {
        margin: 1;
        padding: 0 2;
        color: $text-muted;
        width: 100%;
    }

    Input {
        margin: 0 1;
        width: 98%;
        border: tall $primary-darken-2;
    }

    Input:focus {
        border: tall $accent;
    }

    Button {
        margin: 1;
        min-width: 16;
        background: $primary;
    }

    Button:hover {
        background: $primary-lighten-2;
    }

    .action_buttons {
        dock: left;
    }

    .view_results {
        dock: right;
        background: $success;
    }

    .view_results:hover {
        background: $success-lighten-2;
    }

    LoadingIndicator {
        width: 16;
        height: 3;
        dock: right;
        display: none;
        margin: 1;
        color: $warning;
    }

    .running LoadingIndicator {
        display: block;
    }

    .process-widget:focus-within {
        border: tall $accent;
    }

    .process-widget.running {
        border: tall $warning;
    }
    """



    def __init__(self, processes: dict[str, ProcessConfig]):
        self.processes = processes
        super().__init__()

    def compose(self) -> ComposeResult:
        for process_name, config in self.processes.items():
            with ScrollableContainer (classes="process-widget", id=f"{process_name}_widget"):
                yield Label(config.name, classes="process-title")

                # Input fields
                for field in config.input_fields:
                    yield Input(
                        placeholder=field.replace("_", " ").title(),
                        id=f"{process_name}_{field}_input",
                        suggester=file_suggester,
                    )

                # Description
                if config.description:
                    yield Label(config.description, classes="process-description")

                # Process button
                yield Button(
                    config.name,
                    id=f"{process_name}_button",
                    classes="action_buttons",
                )

                # View results button (only for processes that generate viewable results)
                if process_name in ["fastqc", "multiqc", "bcftools_stats"]:
                    yield Button(
                        "View Results",
                        id=f"view_{process_name}_results",
                        classes="view_results"
                    )

                # Loading indicator
                yield LoadingIndicator(id=f"{process_name}_loading")
