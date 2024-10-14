from textual.app import ComposeResult
from textual.containers import Container, Horizontal
from textual.widgets import Button, Static, Label, Input

from modules.exec_func.common_funcs import file_suggester, get_input


class BcftoolsWidgets(Static):
    def compose(self) -> ComposeResult:

        with Container(id="bam_mpileup_container"):
            yield Label("bamtools mpileup")
            with Horizontal(id="bam_mpileup_horizontal"):
                yield Input( placeholder="reference", id="bam_mpileup_reference_input", suggester=file_suggester,)
                yield Input( placeholder="input bam file name", id="bam_mpileup_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bam file name", id="bam_mpileup_output_input", suggester=file_suggester,)
                yield Button( "mpileup", id="bam_mpileup_button",)
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Convert .bam files to .bam files. you can also select only to convert a certain region", classes="description")

        with Container(id="bam_call_container"):
            yield Label("bamtools call")
            with Horizontal(id="bam_call_horizontal"):
                yield Input( placeholder="input bam/bam file name", id="bam_call_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bam/bam file name", id="bam_call_output_input", suggester=file_suggester,)
                yield Button("call bam file", id="bam_call_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC].bam and .bam files need to be called before futher processing", classes="description")

        with Container(id="bam_filter_container"):
            yield Label("bamtools filter")
            with Horizontal(id="bam_filter_horizontal"):
                yield Input( placeholder="filter", id="bam_filter_filter_input", suggester=file_suggester,)
                yield Input( placeholder="input bam file name", id="bam_filter_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bam file name", id="bam_filter_output_input", suggester=file_suggester,)
                yield Button( "filter", id="bam_filter_button",)
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Convert .bam files to .bam files. you can also select only to convert a certain region", classes="description")

        with Container(id="bam_norm_container"):
            yield Label("bamtools norm")
            with Horizontal(id="bam_norm_horizontal"):
                yield Input( placeholder="input bam/bam file name", id="bam_norm_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bam/bam file name", id="bam_norm_output_input", suggester=file_suggester,)
                yield Button("norm bam file", id="bam_call_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC].bam and .bam files need to be normed before futher processing", classes="description")

        with Container(id="bam_stats_container"):
            yield Label("bamtools stats")
            with Horizontal(id="bam_stats_horizontal"):
                yield Input( placeholder="input bam/bam file name", id="bam_stats_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bam/bam file name", id="bam_stats_output_input", suggester=file_suggester,)
                yield Button("stats bam file", id="bam_stats_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC].bam and .bam files need to be statsed before futher processing", classes="description")

# TODO:

"""
# Bcftools Commands:
1. bcftools mpileup    # Generate VCF/BCF files from alignment files.
2. bcftools call       # Call variants (SNPs/indels).
3. bcftools filter     # Apply filters to VCF/BCF files.
4. bcftools norm       # Normalize variants (e.g., left-align indels).
5. bcftools merge      # Merge multiple VCF/BCF files.
6. bcftools annotate   # Annotate VCF files with additional information.
7. bcftools stats      # Generate statistics for VCF/BCF files.
"""
