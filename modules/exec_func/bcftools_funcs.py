# import subprocess
# import threading
#
#
# # Function to run a command using subprocess and threading for parallel execution
# def run_command(cmd):
#     process = subprocess.Popen(
#         cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
#     )
#     stdout, stderr = process.communicate()
#     if process.returncode != 0:
#         print(f"Error running command: {cmd}\n{stderr.decode()}")
#     else:
#         print(stdout.decode())
#
#
# # 1. bcftools mpileup command
# def bcftools_mpileup(bam_file, reference, output_file):
#     cmd = f"bcftools mpileup -f {reference} {bam_file} -o {output_file}"
#     run_command(cmd)
#
#
# # 2. bcftools call command
# def bcftools_call(input_file, output_file):
#     cmd = f"bcftools call -mv -Oz -o {output_file} {input_file}"
#     run_command(cmd)
#
#
# # 3. bcftools filter command
# def bcftools_filter(input_file, output_file, filter_expression):
#     cmd = f"bcftools filter -s LOWQUAL -e '{filter_expression}' {input_file} -o {output_file}"
#     run_command(cmd)
#
#
# # 4. bcftools norm command
# def bcftools_norm(input_file, output_file):
#     cmd = f"bcftools norm -f {input_file} {output_file}"
#     run_command(cmd)
#
#
# # 5. bcftools stats command
# def bcftools_stats(input_file, output_file):
#     cmd = f"bcftools stats {input_file} > {output_file}"
#     run_command(cmd)
#
#
#
#
import subprocess
import threading


# Function to run a command using subprocess
def run_command(cmd):
    """Runs a shell command using subprocess."""
    process = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error running command: {cmd}\n{stderr.decode()}")
    else:
        print(stdout.decode())


# Thread wrapper to run commands in parallel
def run_command_threaded(cmd):
    """Runs a command in a separate thread."""
    thread = threading.Thread(target=run_command, args=(cmd,))
    thread.start()
    return thread


# Example usage of threading with bcftools commands
def bcftools_mpileup(bam_file, reference, output_file):
    cmd = f"bcftools mpileup -f {reference} {bam_file} -o {output_file}"
    return run_command_threaded(cmd)


def bcftools_call(input_file, output_file):
    cmd = f"bcftools call -mv -Oz -o {output_file} {input_file}"
    return run_command_threaded(cmd)


def bcftools_filter(input_file, output_file, filter_expression):
    cmd = f"bcftools filter -s LOWQUAL -e '{filter_expression}' {input_file} -o {output_file}"
    return run_command_threaded(cmd)


# Example of running commands in parallel
# if __name__ == "__main__":
#     # Define some example inputs
#     bam_file = "example.bam"
#     reference = "reference.fasta"
#     output_mpileup = "output_mpileup.bcf"
#     output_call = "output_call.vcf.gz"
#     filter_expression = "QUAL<20"
#
#     # Run commands in parallel threads
#     mpileup_thread = bcftools_mpileup(bam_file, reference, output_mpileup)
#     call_thread = bcftools_call(output_mpileup, output_call)
#     filter_thread = bcftools_filter(output_call, "filtered_output.vcf.gz", filter_expression)
#
#     # Wait for all threads to finish
#     mpileup_thread.join()
#     call_thread.join()
#     filter_thread.join()
