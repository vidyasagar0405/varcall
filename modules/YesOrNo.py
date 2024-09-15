from textual.screen import ModalScreen
from textual.widgets import Button, Label
from textual.app import ComposeResult
from textual.containers import Container, Horizontal
from textual import on


class YesOrNo(ModalScreen):
    def __init__(self, Title: str = "Input") -> None:
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
