# **Help**

Log file is located at `/tmp/varcall.log`

>
> **NOTE**
> It is recommended to run `varcall` in a conda environment
>

To install conda refer [miniconda](https://docs.anaconda.com/miniconda/)

`varcall` requires variant calling tools to run, you can download them by:
```console
conda create -n varcall -y
conda install -c bioconda fastqc multiqc bwa samtools bcftools
```

## File structure

- Project Name (*working dir*)
    - Data
         - Reference Genome
         - Reads
            - Sample_1_1.fastq
            - Sample_1_2.fastq
            - Sample_2_1.fastq
            - Sample_2_2.fastq
        - BED
    - Results
        - FastQC
            - Raw
            - Trimmed
        - MultiQC
            - Alignment
            - Raw
            - Trimmed
        - Samfile
        - Bamfile
        - VCF
        - BCf

