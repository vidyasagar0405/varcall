from textual.app import ComposeResult
from textual.containers import ScrollableContainer, Horizontal
from textual.widgets import Button, Input, Label, LoadingIndicator, Static

from varcall.process.process_class import ProcessConfig
from varcall.components.file_suggester import file_suggester


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
        border: tall $primary;
        background: $boost;
    }

    .process-widget {
        layout: vertical;
        width: 100%;
        height: auto;
        min-height: 12;
        margin: 1 0;
        padding: 1;
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

    .input-field {
        margin: 1 0;
        width: 98%;
    }

    .process-description {
        margin: 1;
        padding: 0 2;
        color: $text-muted;
        width: 98%;
    }

    .button-container {
        width: 100%;
        content-align: right middle;
        height: auto;
        margin: 1 0;
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
        margin-right: 1;
    }

    .view_results {
        background: $success;
    }

    .view_results:hover {
        background: $success-lighten-2;
    }

    LoadingIndicator {
        width: 16;
        height: 3;
        display: none;
        margin: 1;
        color: $warning;
    }

    .running LoadingIndicator {
        display: block;
    }

    ScrollableContainer:focus-within {
        border: tall $accent;
    }

    ScrollableContainer.running {
        border: tall $warning;
    }
    """

    def __init__(self, processes: dict[str, ProcessConfig]):
        self.processes = processes
        super().__init__()

    def compose(self) -> ComposeResult:
        for process_name, config in self.processes.items():
            with ScrollableContainer(classes="process-widget", id=f"{process_name}_widget"):
                # Title at top
                yield Label(config.name, classes="process-title")

                # Input fields in middle
                for field in config.input_fields:
                    yield Input(
                        placeholder=field.replace("_", " ").title(),
                        id=f"{process_name}_{field}_input",
                        suggester=file_suggester,
                        classes="input-field"
                    )

                # Description
                if config.description:
                    yield Label(config.description, classes="process-description")

                # Button container for horizontal alignment
                with Horizontal(classes="button-container"):
                    yield LoadingIndicator(id=f"{process_name}_loading")

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
