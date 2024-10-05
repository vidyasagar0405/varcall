import pysam

from modules import logging
from modules.logging import setup_logging
from pathlib import Path

setup_logging()


# 1. samtools view equivalent in pysam
def view_bam(input_file, output_file, region=None):
    """
    The view_bam function reads a BAM file, optionally filters the
    reads based on a specified genomic region,
    and writes the filtered (or unfiltered) reads to an output BAM file.

    Also converts a SAM file input to BAM file
    """
    with pysam.AlignmentFile(input_file, "rb") as infile, pysam.AlignmentFile(
        output_file, "wb", header=infile.header
    ) as outfile:
        for read in infile.fetch(region=region) if region else infile:
            outfile.write(read)


# 2. samtools sort equivalent in pysam
def sort_bam(input_file, output_file):
    """
    The sort_bam function sorts the input BAM file by coordinate order and
    writes the sorted data to the output BAM file.
    This is equivalent to `samtools sort` which sorts alignments to ensure
    they are in order by reference positions.
    """
    pysam.sort("-o", output_file, input_file)
    logging.info("async sort completed")


# 3. samtools index equivalent in pysam
def index_bam(input_file):
    """
    The index_bam function creates an index file for the input BAM file,
    which allows for quick access to data from specific regions.
    This is equivalent to `samtools index` and generates a .bai file
    which is essential for many downstream analyses.
    """
    pysam.index(input_file)


# 4. samtools flagstat equivalent in pysam
def flagstat_bam(input_file, out_file):
    """
    The flagstat_bam function generates a summary of alignment statistics
    for the input BAM file, such as total number of reads, mapped reads,
    and properly paired reads.
    This is equivalent to `samtools flagstat` and is useful for assessing
    the quality of alignments.
    """
    stats = pysam.flagstat(input_file)
    out_file = Path(out_file)
    out_file.touch(exist_ok=True)
    with open(out_file, "w") as outfile:
        outfile.write(str(stats))


# 5. samtools stats equivalent in pysam
def stats_bam(input_file, output_file):
    """
    The stats_bam function calculates various statistics from the input BAM file
    and writes the results to the output file. The statistics include
    information on read lengths, GC content, and quality scores, among others.
    This is equivalent to `samtools stats`, providing detailed insights
    into the sequencing data.
    """
    with open(output_file, "w") as out:
        stats = str(pysam.stats(input_file))
        out.write(stats)


# if __name__ == "__main__":
#     sort_bam("aligned.sam", "aligned.sorted.sam")
#     index_bam("aligned.sorted.sam")
