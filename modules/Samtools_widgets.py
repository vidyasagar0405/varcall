import threading
from pysam.utils import SamtoolsError
import logging
from textual.app import ComposeResult, on
from textual.widgets import Button, Input, Label, Static
from textual.containers import Container, Horizontal
from modules.exec_func.samtools_funcs import *  # noqa: F403
from modules.exec_func.common_funcs import get_input
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
    #
    # @on(Button.Pressed, "#sam_view_button")
    # def sam_view(self) -> None:
    #     if not self.workingDir:
    #         self.workingDir = __file__
    #     input = get_input(self, "sam_view_input_input")
    #     ouput = get_input(self, "sam_view_output_input")
    #     sam_range = get_input(self, "sam_view_range_input")
    #     if not sam_range:
    #         sam_range = None
    #     self.notify("samtools view starting", title="Samtools view")
    #     try:
    #         threading.Thread( target=view_bam, args=( input, ouput, sam_range,),).start()  # noqa: F405
    #         self.notify("samtools view completed", title="Samtools view")
    #         logging.info
    #     except SamtoolsError as err:
    #         self.notify(f"An error occured: {err}", title="Samtools view")
    #         logging.error(f"An error occured: {err}")
    #
    # @on(Button.Pressed, "#sam_sort_button")
    # def sam_sort(self) -> None:
    #     if not self.workingDir:
    #         self.workingDir = __file__
    #     input = get_input(self, "sam_sort_input_input")
    #     output = get_input(self, "sam_sort_output_input")
    #     self.notify("samtools sort starting", title="Samtools sort")
    #     try:
    #         threading.Thread( target=sort_bam, args=( input, output,),).start()  # noqa: F405
    #         self.notify("samtools sort completed", title="Samtools sort")
    #     except SamtoolsError as err:
    #         self.notify(f"An error occured: {err}", title="Samtools sort")
    #         logging.error(f"An error occured: {err}")
    #
    # @on(Button.Pressed, "#sam_index_button")
    # def sam_index(self) -> None:
    #     if not self.workingDir:
    #         self.workingDir = __file__
    #     input = get_input(self, "sam_index_input_input")
    #     self.notify("samtools index starting", title="Samtools index")
    #     try:
    #         threading.Thread(target=index_bam, args=(input,)).start()  # noqa: F405
    #         self.notify("samtools index completed", title="Samtools index")
    #     except SamtoolsError as err:
    #         self.notify(f"An error occured: {err}", title="Samtools index")
    #         logging.error(f"An error occured: {err}")
