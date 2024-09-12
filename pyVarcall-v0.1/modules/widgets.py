from pathlib import Path
from textual.containers import ScrollableContainer, Horizontal
from textual.widgets import Button, Static, Input, Label
from textual.suggester import Suggester

class FileSuggester(Suggester):
    async def get_suggestion(self, value: str) -> str | None:
        # Handle absolute paths
        if value.startswith('/'):
            return " "
        # Handle invalid pattern '**' anywhere in the input
        if '*' in value:
            return " "
        path = next(Path().glob(f"{value}*"), None)
        return str(path) if path else None

class InputProjectName(Static):

    def compose(self):
        with Horizontal(id="InputProjectName"):
            yield Label("Enter project name (Enter absolute path)", id="Project_name")
            yield Input(id="Input_Project_name")
            yield Label("This will be your working directory", id="Project_location")

class InputBlock(Static):
    def __init__(self, BlockTitle: str="Input", Name: str="Input") -> None:
        self.BlockTitle = BlockTitle
        self.Name = Name
        super().__init__()

    def compose(self):
        yield Label(f"{self.BlockTitle}", id="InputBlockTitle")
        with Horizontal(id="InputurlContainer"):
            yield Input(placeholder="Enter url (uses curl to fetch files)", id=f"Input_{self.Name}_url", classes="Input_url")
            yield Label("-o", id="-o")
            yield Input(placeholder="output file name", id=f"Input_{self.Name}_output_file", suggester=FileSuggester(use_cache=False), classes="Input_output_file")
            yield Button("Download", id="download_Button")
        with Horizontal(id="Inputfile"):
            yield Input(placeholder="Enter file name", id=f"Input_{self.Name}_file_name", suggester=FileSuggester(use_cache=False), classes="Input_file_name")
            yield Button("select dir", id="Select_Directory_Button")
        yield Label("Selected directory:", id=f"Label_selected_{self.Name}_directory")
