from textual.app import ComposeResult
from textual.widgets import Button, Static


class BcftoolsWidgets(Static):
    def compose(self) -> ComposeResult:
        yield Button()
        yield Button()
        yield Button()
