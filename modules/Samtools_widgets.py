from textual.app import ComposeResult
from textual.widgets import Button, Static


class SamtoolsWidgets(Static):
    def compose(self) -> ComposeResult:
        yield Button()
        yield Button()


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
