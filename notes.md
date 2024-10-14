# Variant calling paper ref

## The evaluation of Bcftools mpileup and GATK HaplotypeCaller for variant calling in non‑human species 

[Paper DOI link](https://doi.org/10.1038/s41598-022-15563-2)

**Steps:**

- aligned using bwa-mem2
    - `bwa-mem2 mem ref.genome.fa read1.fq.gz read2.fq.gz|samtools view -bS |samtools sort -o Bam/a1.bam`

- he resulting bam files were indexed using samtools

- `bcftools mpileup -f REFERENCE LIST_OF_BAM | bcftools call -mv -Oz -o VCFFILE`

- Filtering was performed by discarding SNVs in which the variant calling score at QUAL field is lower than 20
