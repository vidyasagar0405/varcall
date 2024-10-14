import subprocess
import threading


# Function to run a command using subprocess
def run_command(cmd):
    """
    Runs a shell command using subprocess.
    Executes the provided command in a new process, captures the output,
    and prints the results or errors if any occur.
    """
    # Create a new process to run the command
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    # Communicate with the process to capture stdout and stderr
    stdout, stderr = process.communicate()
    # Check if the command was successful
    if process.returncode != 0:
        print(f"Error running command: {cmd}\n{stderr.decode()}")
    else:
        print(stdout.decode())


# Thread wrapper to run commands in parallel
def run_command_threaded(cmd):
    """
    Runs a command in a separate thread.
    This allows commands to be executed in parallel, improving efficiency
    for tasks that can be run independently.
    """
    # Create a new thread to run the command
    thread = threading.Thread(target=run_command, args=(cmd,))
    # Start the thread, which will run the command asynchronously
    thread.start()
    return thread


# Example usage of threading with bcftools commands


# 1. bcftools mpileup command
def bcftools_mpileup(bam_file, reference, output_file):
    """
    Runs the bcftools mpileup command to generate the pileup format
    from the BAM file against a reference genome.
    This is used to gather information about the sequence data for variant calling.
    """
    # Construct the command string for mpileup
    cmd = f"bcftools mpileup -f {reference} {bam_file} -o {output_file}"
    # Run the command in a separate thread
    return run_command_threaded(cmd)


# 2. bcftools call command
def bcftools_call(input_file, output_file):
    """
    Runs the bcftools call command to perform variant calling
    on the input file (typically an output from mpileup).
    This step identifies variants such as SNPs and indels.
    """
    # Construct the command string for variant calling
    cmd = f"bcftools call -mv -Oz -o {output_file} {input_file}"
    # Run the command in a separate thread
    return run_command_threaded(cmd)


# 3. bcftools filter command
def bcftools_filter(input_file, output_file, filter_expression):
    """
    Runs the bcftools filter command to apply quality and content filters
    on the variant call file. Filters help in removing low-quality or
    unwanted variants based on the specified criteria.
    """
    # Construct the command string for filtering variants
    cmd = f"bcftools filter -s LOWQUAL -e '{filter_expression}' {input_file} -o {output_file}"
    # Run the command in a separate thread
    return run_command_threaded(cmd)


# 4. bcftools norm command
def bcftools_norm(input_file, output_file):
    """
    Runs the bcftools norm command to normalize variant calls.
    This step adjusts representation of indels, left-aligns variants,
    splits multi-allelic sites into bi-allelic records, and checks
    consistency with the reference genome.

    Args:
        input_file (str): Path to the input VCF or BCF file that contains
                          variant calls.
        output_file (str): Path to the output file where the normalized
                           variant data will be saved.

    The normalization process helps standardize variant calls and is an
    essential step before downstream analysis like annotation or merging
    multiple VCF files.
    """
    # Construct the command string for normalizing variants
    cmd = f"bcftools norm -f {input_file} {output_file}"
    # Run the command using subprocess to perform the normalization
    run_command(cmd)


# 5. bcftools stats command
def bcftools_stats(input_file, output_file):
    """
    Runs the bcftools stats command to generate statistical summaries
    of the variant data in the input file. The statistics include counts
    of SNPs, indels, transitions, transversions, and other key metrics
    that provide an overview of the variant calling results.

    Args:
        input_file (str): Path to the input VCF or BCF file that contains
                          the variant calls to be analyzed.
        output_file (str): Path to the output file where the statistics
                           will be written.

    This command is useful for assessing the quality and characteristics
    of the variant calls produced by your pipeline.
    """
    # Construct the command string for generating statistics
    cmd = f"bcftools stats {input_file} > {output_file}"
    # Run the command using subprocess to generate the stats
    run_command(cmd)
