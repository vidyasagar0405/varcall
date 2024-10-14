from textual.app import ComposeResult
from textual.containers import Container, Horizontal
from textual.widgets import Button, Static, Label, Input

from modules.exec_func.common_funcs import file_suggester


class BcftoolsWidgets(Static):
    def compose(self) -> ComposeResult:

        with Container(id="bcf_mpileup_container"):
            yield Label("bcftools mpileup")
            with Horizontal(id="bcf_mpileup_horizontal"):
                yield Input( placeholder="reference", id="bcf_mpileup_reference_input", suggester=file_suggester,)
                yield Input( placeholder="input bcf file name", id="bcf_mpileup_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bcf file name", id="bcf_mpileup_output_input", suggester=file_suggester,)
                yield Button( "mpileup", id="bcf_mpileup_button",)
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Generate VCF or BCF containing genotype likelihoods for bcf or CRAM files.", classes="description")

        with Container(id="bcf_call_container"):
            yield Label("bcftools call")
            with Horizontal(id="bcf_call_horizontal"):
                yield Input( placeholder="input bcf/vcf file name", id="bcf_call_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bcf/vcf file name", id="bcf_call_output_input", suggester=file_suggester,)
                yield Button("call bcf file", id="bcf_call_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]compiles the genotype likelihoods and calls only the genetic variants", classes="description")

        with Container(id="bcf_filter_container"):
            yield Label("bcftools filter")
            with Horizontal(id="bcf_filter_horizontal"):
                yield Input( placeholder="filter", id="bcf_filter_filter_input", suggester=file_suggester,)
                yield Input( placeholder="input bcf file name", id="bcf_filter_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bcf file name", id="bcf_filter_output_input", suggester=file_suggester,)
                yield Button( "filter", id="bcf_filter_button",)
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Apply filters to the called bcf/vcf data", classes="description")

        with Container(id="bcf_norm_container"):
            yield Label("bcftools norm")
            with Horizontal(id="bcf_norm_horizontal"):
                yield Input( placeholder="input bcf/bcf file name", id="bcf_norm_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bcf/bcf file name", id="bcf_norm_output_input", suggester=file_suggester,)
                yield Button("norm bcf file", id="bcf_call_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Normalise the data", classes="description")

        with Container(id="bcf_stats_container"):
            yield Label("bcftools stats")
            with Horizontal(id="bcf_stats_horizontal"):
                yield Input( placeholder="input bcf/bcf file name", id="bcf_stats_input_input", suggester=file_suggester,)
                yield Input( placeholder="output bcf/bcf file name", id="bcf_stats_output_input", suggester=file_suggester,)
                yield Button("stats bcf file", id="bcf_stats_button")
            yield Label("[bold #01B0DC]Description: [/bold #01B0DC]Get the stats of bvf/vcf files", classes="description")

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
