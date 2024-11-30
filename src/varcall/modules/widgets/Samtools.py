from textual.app import ComposeResult
from textual.widgets import Button, Input, Label, Static
from textual.containers import Grid
from varcall.modules.exec.samtools  import *  # noqa: F403
from varcall.modules.exec.common  import file_suggester
from varcall.modules.logging import setup_logging

setup_logging()


class SamtoolsWidgets(Static):

    def compose(self) -> ComposeResult:

        with Grid(id="sam_view_container"):
            yield Label("Samtools view")

            yield Input( placeholder="input sam/bam file name", id="sam_view_input_input", suggester=file_suggester,)
            yield Input( placeholder="output bam file name", id="sam_view_output_input", suggester=file_suggester,)
            yield Input( placeholder="region (optional)", id="sam_view_region_input", suggester=file_suggester,)
            yield Button( "Convert sam to bam", id="sam_view_button",)
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Convert .sam files to .bam files. you can also select only to convert a certain region", classes="description")

        with Grid(id="sam_sort_container"):
            yield Label("Samtools sort")

            yield Input( placeholder="input sam/bam file name", id="sam_sort_input_input", suggester=file_suggester,)
            yield Input( placeholder="output sam/bam file name", id="sam_sort_output_input", suggester=file_suggester,)
            yield Button("sort bam file", id="sam_sort_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC].sam and .bam files need to be sorted before futher processing", classes="description")

        with Grid(id="sam_index_container"):
            yield Label("Samtools index")

            yield Input( placeholder="input bam file name", id="sam_index_input_input", suggester=file_suggester,)
            yield Button("index bam file", id="sam_index_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]index the .bam file for futher analysis", classes="description")
