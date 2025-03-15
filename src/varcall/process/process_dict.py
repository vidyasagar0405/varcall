from varcall.process.process_class import ProcessConfig

# Dictionary to track running processes
running_processes = {}

# Configure bioinformatics tools
HOME_PROCESSES = {
    # Data download tool
    "curl": ProcessConfig(
        name="Data Download",
        command="curl -L {url} -o {output_file}",
        input_fields=["url", "output_file"],
        required_fields=["url", "output_file"],
        description="FastQCs the files in the workingdir/data/reads directory [u]if the input field is left empty[/u]",
        success_message="Dataset downloaded successfully",
        error_message="Dataset download failed",
    ),
    "fastqc": ProcessConfig(
        name="FastQC",
        command="fastqc {input_file} -o {output_dir}",
        input_fields=["input_file", "output_dir"],
        required_fields=["input_file"],
        description="Quality control tool for high throughput sequence data",
        success_message="FastQC analysis completed successfully",
        error_message="FastQC analysis failed",
    ),
    "multiqc": ProcessConfig(
        name="MultiQC",
        command="multiqc {input_dir} -o {output_dir}",
        input_fields=["input_dir", "output_dir"],
        required_fields=["input_dir"],
        description="Aggregate results from bioinformatics analyses into a single report",
        success_message="MultiQC report generated successfully",
        error_message="MultiQC report generation failed",
    ),
    "bwa_index": ProcessConfig(
        name="BWA Index",
        command="bwa index {reference_genome}",
        input_fields=["reference_genome"],
        required_fields=["reference_genome"],
        description="Index reference genome for BWA alignment",
        success_message="Reference genome indexed successfully",
        error_message="Reference genome indexing failed",
    ),
    "bwa_mem": ProcessConfig(
        name="BWA MEM Alignment",
        command="bwa mem {reference_genome} {read1} {read2} > {output_sam}",
        input_fields=["reference_genome", "read1", "read2", "output_sam"],
        required_fields=["reference_genome", "read1", "output_sam"],
        description="Align reads to reference genome using BWA MEM algorithm",
        success_message="Read alignment completed successfully",
        error_message="Read alignment failed",
        default_output_ext=".sam",
    ),
}
SAMTOOLS_PROCESSES = {
    "samtools_view": ProcessConfig(
        name="SAM to BAM Conversion",
        command="samtools view -bS {input_sam} > {output_bam}",
        input_fields=["input_sam", "output_bam"],
        required_fields=["input_sam", "output_bam"],
        description="Convert SAM file to BAM format",
        success_message="SAM to BAM conversion completed successfully",
        error_message="SAM to BAM conversion failed",
        default_output_ext=".bam",
    ),
    "bamtools_stats": ProcessConfig(
        name="BAMTools Stats",
        command="bamtools stats -in {input_bam}",
        input_fields=["input_bam"],
        required_fields=["input_bam"],
        description="Generate statistics for BAM file",
        success_message="BAM statistics generated successfully",
        error_message="BAM statistics generation failed",
        default_output_ext=".txt",
    ),
    # Additional samtools operations
    "samtools_sort": ProcessConfig(
        name="BAM Sorting",
        command="samtools sort {input_bam} -o {output_sorted_bam}",
        input_fields=["input_bam", "output_sorted_bam", "threads"],
        required_fields=["input_bam", "output_sorted_bam"],
        description="Sort BAM file by coordinates",
        success_message="BAM file sorted successfully",
        error_message="BAM sorting failed",
        default_output_ext=".bam",
    ),
    "samtools_index": ProcessConfig(
        name="BAM Indexing",
        command="samtools index {input_bam} {output_bai}",
        input_fields=["input_bam", "output_bai"],
        required_fields=["input_bam"],
        description="Index BAM file for fast random access",
        success_message="BAM file indexed successfully",
        error_message="BAM indexing failed",
    ),
}

BCFTOOLS_PROCESS = {
    # bcftools operations
    "bcftools_mpileup": ProcessConfig(
        name="BCFtools Mpileup",
        command="bcftools mpileup -f {reference_genome} {input_bam} -o {output_vcf}",
        input_fields=[
            "reference_genome",
            "input_bam",
            "output_vcf",
            "additional_params",
        ],
        required_fields=["reference_genome", "input_bam", "output_vcf"],
        description="Generate VCF file containing genotype likelihoods for variants",
        success_message="Variant calling (mpileup) completed successfully",
        error_message="Variant calling (mpileup) failed",
        default_output_ext=".vcf",
    ),
    "bcftools_call": ProcessConfig(
        name="BCFtools Call",
        command="bcftools call -m -v {input_vcf} -o {output_vcf}",
        input_fields=["input_vcf", "output_vcf", "additional_params"],
        required_fields=["input_vcf", "output_vcf"],
        description="Call variants from VCF file containing genotype likelihoods",
        success_message="Variant calling completed successfully",
        error_message="Variant calling failed",
        default_output_ext=".vcf",
    ),
    "bcftools_filter": ProcessConfig(
        name="BCFtools Filter",
        command="bcftools filter -i'{filter_expression}' {input_vcf} -o {output_vcf}",
        input_fields=["input_vcf", "output_vcf", "filter_expression"],
        required_fields=["input_vcf", "output_vcf", "filter_expression"],
        description="Filter variants in VCF file based on specified criteria",
        success_message="Variant filtering completed successfully",
        error_message="Variant filtering failed",
        default_output_ext=".vcf",
    ),
}

GATK_PROCESS = {
    # GATK tools
    "gatk_haplotypecaller": ProcessConfig(
        name="GATK HaplotypeCaller",
        command="gatk HaplotypeCaller -R {reference_genome} -I {input_bam} -O {output_vcf}",
        input_fields=[
            "reference_genome",
            "input_bam",
            "output_vcf",
            "additional_params",
        ],
        required_fields=["reference_genome", "input_bam", "output_vcf"],
        description="Call germline SNPs and indels via local assembly of haplotypes",
        success_message="Variant calling with HaplotypeCaller completed successfully",
        error_message="Variant calling with HaplotypeCaller failed",
        default_output_ext=".vcf",
    ),
    "gatk_baserecalibrator": ProcessConfig(
        name="GATK BaseRecalibrator",
        command="gatk BaseRecalibrator -R {reference_genome} -I {input_bam} --known-sites {known_sites} -O {output_table}",
        input_fields=["reference_genome", "input_bam", "known_sites", "output_table"],
        required_fields=[
            "reference_genome",
            "input_bam",
            "known_sites",
            "output_table",
        ],
        description="Generate recalibration table for base quality score recalibration",
        success_message="Base recalibration table generated successfully",
        error_message="Base recalibration table generation failed",
    ),
    "gatk_applybqsr": ProcessConfig(
        name="GATK ApplyBQSR",
        command="gatk ApplyBQSR -R {reference_genome} -I {input_bam} --bqsr-recal-file {recalibration_table} -O {output_bam}",
        input_fields=[
            "reference_genome",
            "input_bam",
            "recalibration_table",
            "output_bam",
        ],
        required_fields=[
            "reference_genome",
            "input_bam",
            "recalibration_table",
            "output_bam",
        ],
        description="Apply base quality score recalibration",
        success_message="Base quality score recalibration applied successfully",
        error_message="Base quality score recalibration failed",
        default_output_ext=".bam",
    ),
    "gatk_variantfiltration": ProcessConfig(
        name="GATK VariantFiltration",
        command="gatk VariantFiltration -R {reference_genome} -V {input_vcf} --filter-expression '{filter_expression}' --filter-name '{filter_name}' -O {output_vcf}",
        input_fields=[
            "reference_genome",
            "input_vcf",
            "output_vcf",
            "filter_expression",
            "filter_name",
        ],
        required_fields=[
            "reference_genome",
            "input_vcf",
            "output_vcf",
            "filter_expression",
            "filter_name",
        ],
        description="Filter variants based on specified criteria",
        success_message="Variant filtration completed successfully",
        error_message="Variant filtration failed",
        default_output_ext=".vcf",
    ),
    "gatk_selectvariants": ProcessConfig(
        name="GATK SelectVariants",
        command="gatk SelectVariants -R {reference_genome} -V {input_vcf} --select-type-to-include {variant_type} -O {output_vcf}",
        input_fields=[
            "reference_genome",
            "input_vcf",
            "output_vcf",
            "variant_type",
            "additional_params",
        ],
        required_fields=["reference_genome", "input_vcf", "output_vcf"],
        description="Select a subset of variants from a VCF file",
        success_message="Variant selection completed successfully",
        error_message="Variant selection failed",
        default_output_ext=".vcf",
    ),
}

PIPLELINES = {
    "mpliup_pipeline": ProcessConfig(
        name = "Mpliup Pipeline (Mock)",
        command="echo {read_1} > {read_2}",
        input_fields=[
            "reference_genome",
            "read_1",
            "read_2",
            "additional_params",
        ],
        required_fields=["reference_genome", "read_1", "read_2" ],
        description="Mpileup Pipeline",
        success_message="Mpileup Pipeline completed successfully",
        error_message="An error occured during the execution of the Mpileup Pipeline",
        default_output_ext=".vcf",
    )
}


MASTER_CONFIG: dict[str, dict[str, ProcessConfig]] = {
    "home": HOME_PROCESSES,
    "samtools": SAMTOOLS_PROCESSES,
    "bcftools": BCFTOOLS_PROCESS,
    "gatk": GATK_PROCESS,
    "pipeline": PIPLELINES,
}
