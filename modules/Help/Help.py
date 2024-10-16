from pathlib import Path
from textual.widgets import Markdown
from textual.app import ComposeResult


# Function to read the Markdown file
def read_markdown_file(file_path: str) -> str:
    with open(file_path, "r") as file:
        return file.read()


# Specify the path to the Help.md file
help_markdown_path = Path("./modules/Help/Help.md")

# Read the content of the Markdown file
help_markdown_content = read_markdown_file(help_markdown_path)


class HelpMarkdown(Markdown):
    def compose(self) -> ComposeResult:
        yield Markdown(help_markdown_content)
