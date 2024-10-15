from textual.widgets import Input
from textual.suggester import Suggester
from pathlib import Path
import glob


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


def find_matching_files(path_pattern: str) -> str:
    """Find all files matching the given wildcard path pattern."""
    matching_files = glob.glob(path_pattern)
    return (" ".join(matching_files))


def get_input(self, input_widget_id) -> str:
    raw_input = str(self.query_one(f"#{input_widget_id}", Input).value).strip()
    return find_matching_files(raw_input)
