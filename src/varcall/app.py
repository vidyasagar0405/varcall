from pathlib import Path

from textual.app import App, ComposeResult
from textual.containers import ScrollableContainer
from textual.widgets import ( Button, Footer, Header, MarkdownViewer,
                              TabbedContent, TabPane,
                            )

from varcall.components.input_block import ProcessWidgets
from varcall.components.yes_or_no import YesOrNo
from varcall.process.process_class import Process
from varcall.config.yaml_parcing import MASTER_CONFIG
from varcall.help import help_path


class Varcall(App[None]):
    # Key bindings for the application
    BINDINGS = [
        ("f1", "show_help", "Help"),
        ("Q", "exit_app", "Exit App"),
        ("T", "change_theme", "Change Theme"),
        ("P", "print_theme_data", "Print Theme data"),
    ]

    CSS_PATH = "./app.css"

    def compose(self) -> ComposeResult:
        with ScrollableContainer(id="ScrollableContainer"):
            yield Header()
            with TabbedContent():
                for tab, processes_dict in MASTER_CONFIG.items():
                    tab = tab.title()
                    with TabPane(tab, id=f"{tab}Tab"):
                        yield ProcessWidgets(processes_dict)
                with TabPane("Help", id="HelpTab"):
                    yield MarkdownViewer(Path(help_path).read_text())
            yield Footer(show_command_palette=True)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        button_id = event.button.id
        if button_id and button_id.endswith("_button"):
            process_name = button_id.replace("_button", "")
            for tab in MASTER_CONFIG.values():
                if process_name in tab:
                    Process(self, tab[process_name]).run()

    # Confirm exit application
    def maybe_exit_app(self, bool) -> None:
        if bool:
            self.exit()

    def action_exit_app(self) -> None:
        self.push_screen(
            YesOrNo("Do you want to exit application?"), self.maybe_exit_app
        )

    # Handels the key press for Toggle dark mode
    def action_change_theme(self) -> None:
        self.search_themes()

    # Handels the key press for Show help action
    def action_show_help(self) -> None:
        self.query_one(TabbedContent).active = "HelpTab"

    def action_print_theme_data(self):
        with open("themes.txt", "a") as f:
            f.write(f"\n\n{self.app.theme}\n\n")
            for key, value in self.app.theme_variables.items():
                f.write(f"{key}: {value}\n")


# Main entry point of the application
def main():
    Varcall().run()


if __name__ == "__main__":
    main()
