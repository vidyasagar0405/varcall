from textual.app import ComposeResult
from textual.widgets import Button, Input, Label, Static
from textual.containers import Container, Horizontal
from modules.exec_func.samtools_funcs import *  # noqa: F403
from modules.exec_func.common_funcs import file_suggester
from modules.logging import setup_logging

setup_logging()


class SamtoolsWidgets(Static):

    def compose(self) -> ComposeResult:

        with Container(id="sam_view_container"):
            yield Label("Samtools view")
            with Horizontal(id="sam_view_horizontal"):
                yield Input( placeholder="input sam/bam file name", id="sam_view_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bam file name", id="sam_view_output_input", suggester=file_suggester,)
                yield Input( placeholder="region (optional)", id="sam_view_region_input", suggester=file_suggester,)
                yield Button( "Convert sam to bam", id="sam_view_button",)
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Convert .sam files to .bam files. you can also select only to convert a certain region", classes="description")

        with Container(id="sam_sort_container"):
            yield Label("Samtools sort")
            with Horizontal(id="sam_sort_horizontal"):
                yield Input( placeholder="input sam/bam file name", id="sam_sort_input_input", suggester=file_suggester,)
                yield Input( placeholder="output sam/bam file name", id="sam_sort_output_input", suggester=file_suggester,)
                yield Button("sort bam file", id="sam_sort_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC].sam and .bam files need to be sorted before futher processing", classes="description")

        with Container(id="sam_index_container"):
            yield Label("Samtools index")
            with Horizontal(id="sam_index_horizontal"):
                yield Input( placeholder="input bam file name", id="sam_index_input_input", suggester=file_suggester,)
                yield Button("index bam file", id="sam_index_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]index the .bam file for futher analysis", classes="description")
