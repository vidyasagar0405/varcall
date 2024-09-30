from textual.widgets import Input
from textual.suggester import Suggester
from pathlib import Path


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


def get_input(self, input_widget_id):
    return str(self.query_one(f"#{input_widget_id}", Input).value).strip()
