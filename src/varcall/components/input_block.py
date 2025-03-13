from pathlib import Path
from textual.app import ComposeResult
from textual.containers import Grid
from textual.suggester import Suggester
from textual.widgets import Button, Input, Label, LoadingIndicator, Static

from varcall.process.process_dict import PROCESSES


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
    def compose(self) -> ComposeResult:
        for process_name, process_config in PROCESSES.items():
            with Grid(id=f"{process_name}_widget"):
                yield Label(process_config.name, id=f"{process_name}_title")
                for field in process_config.input_fields:
                    yield Input(
                        placeholder=field.replace("_", " ").title(),
                        id=f"{process_name}_{field}_input",
                        suggester=file_suggester,
                    )

                yield Label(
                    process_config.description, id=f"{process_name}_description"
                )
                yield Button(
                    process_config.name,
                    id=f"{process_name}_button",
                    classes="action_buttons",
                )
                yield Button("View Results", id="view_FastQC_res")
                yield LoadingIndicator(id=f"{process_name}_loading")
