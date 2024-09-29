from textual.app import ComposeResult
from textual.widgets import Button, Static


class PipelineWidgets(Static):
    def compose(self) -> ComposeResult:
        yield Button("Pipeline")
