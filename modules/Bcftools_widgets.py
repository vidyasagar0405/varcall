from textual.app import ComposeResult
from textual.widgets import Button, Static


class BcftoolsWidgets(Static):
    def compose(self) -> ComposeResult:
        yield Button()
        yield Button()
        yield Button()


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
