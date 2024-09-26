import pysam


# 1. samtools view equivalent in pysam
def view_bam(input_file, output_file, region=None):
    with pysam.AlignmentFile(input_file, "rb") as infile, pysam.AlignmentFile(
        output_file, "wb", header=infile.header
    ) as outfile:
        for read in infile.fetch(region=region) if region else infile:
            outfile.write(read)


# 2. samtools sort equivalent in pysam
def sort_bam(input_file, output_file):
    pysam.sort("-o", output_file, input_file)


# 3. samtools index equivalent in pysam
def index_bam(input_file):
    pysam.index(input_file)


# 4. samtools flagstat equivalent in pysam
def flagstat_bam(input_file):
    stats = pysam.flagstat(input_file)
    print(stats)


# 5. samtools stats equivalent in pysam
def stats_bam(input_file, output_file):
    with open(output_file, "w") as out:
        stats = str(pysam.stats(input_file))
        out.write(stats)
