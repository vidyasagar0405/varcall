import threading
from pysam.utils import SamtoolsError
import logging
from textual.app import ComposeResult, on
from textual.widgets import Button, Input, Label, Static
from textual.containers import Container, Horizontal
from modules.exec_func.samtools_funcs import *  # noqa: F403
from modules.exec_func.hometab_funcs import FileSuggester
from modules.logging import setup_logging

setup_logging()

class SamtoolsWidgets(Static):
    
    file_suggester = FileSuggester(use_cache=False, case_sensitive=True)

    def compose(self) -> ComposeResult:
        with Container(id="sam_view_container"):
            yield Label("Samtools view")
            with Horizontal(id="sam_view_horizontal"):
                yield Input(placeholder="input sam/bam file name", id="sam_view_input_input", suggester=self.file_suggester) 
                yield Input(placeholder="output bam file name", id="sam_view_output_input", suggester=self.file_suggester)
                yield Input(placeholder="range (optional)", id="sam_view_range_input", suggester=self.file_suggester)
            yield Button("Convert sam to bam", id="sam_view_button",)

        with Container(id="sam_sort_container"):
            yield Label("Samtools sort")
            with Horizontal(id="sam_sort_horizontal"):
                yield Input(placeholder="input sam/bam file name", id="sam_sort_input_input", suggester=self.file_suggester)
                yield Input(placeholder="output sam/bam file name", id="sam_sort_output_input", suggester=self.file_suggester)
            yield Button("sort bam file", id="sam_sort_button")

        with Container(id="sam_index_container"):
            yield Label("Samtools index")
            with Horizontal(id="sam_index_horizontal"):
                yield Input(placeholder="input bam file name", id="sam_index_input_input", suggester=self.file_suggester)
            yield Button("index bam file", id="sam_index_button")

    @on(Button.Pressed, "#sam_view_button")
    def sam_view(self) -> None:
        input = str(self.query_one("#sam_view_input_input", Input).value).strip()
        ouput = str(self.query_one("#sam_view_output_input", Input).value).strip()
        sam_range = str(self.query_one("#sam_view_range_input", Input).value).strip()
        if not sam_range:
            sam_range = None
        self.notify("samtools view starting", title="Samtools view")
        try:
            threading.Thread(target=view_bam, args=(input, ouput, sam_range,)).start()  # noqa: F405
            self.notify("samtools view completed", title="Samtools view")
        except SamtoolsError as err:
            self.notify(f"An error occured: {err}", title="Samtools view")
            logging.error(f"An error occured: {err}")


    @on(Button.Pressed, "#sam_sort_button")
    def sam_sort(self) -> None:
        input = str(self.query_one("#sam_sort_input_input", Input).value).strip()
        ouput = str(self.query_one("#sam_sort_output_input", Input).value).strip()
        self.notify("samtools sort starting", title="Samtools sort")
        try:
            threading.Thread(target=sort_bam, args=(input, ouput,)).start()  # noqa: F405
            self.notify("samtools sort completed", title="Samtools sort")
        except SamtoolsError as err:
            self.notify(f"An error occured: {err}", title="Samtools sort")
            logging.error(f"An error occured: {err}")


    @on(Button.Pressed, "#sam_index_button")
    def sam_index(self) -> None:
        input = str(self.query_one("#sam_index_input_input", Input).value).strip()
        self.notify("samtools index starting", title="Samtools index")
        try:
            threading.Thread(target=index_bam, args=(input,)).start()  # noqa: F405
            self.notify("samtools index completed", title="Samtools index")
        except SamtoolsError as err:
            self.notify(f"An error occured: {err}", title="Samtools index")
            logging.error(f"An error occured: {err}")

# TODO:

"""
# Samtools Commands:
1. samtools view       # Convert, filter, or view alignment files (BAM/CRAM/SAM).
2. samtools sort       # Sort alignment files by coordinates or read names.
3. samtools index      # Index a BAM or CRAM file to allow fast access to specific regions.
4. samtools mpileup    # Generate pileup data for variant calling (deprecated in favor of bcftools mpileup).
5. samtools flagstat   # Report alignment statistics.
6. samtools stats      # Generate detailed statistics on the alignments.
7. samtools depth      # Compute the depth of coverage.
"""
