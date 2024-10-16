import os.path
import pysam
from pysam.utils import SamtoolsError
import asyncio
import threading
import logging
from modules.logging import setup_logging
from pathlib import Path
from modules.exec.common  import get_input

setup_logging()

def get_basename_and_ext(path) -> tuple:
    filename = os.path.splitext(os.path.basename(path))
    return filename

async def run_samtools_sort(self):
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "sam_sort_input_input")
    output_path = get_input(self, "sam_sort_output_input")

    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools sort")
        return

    if not output_path:
        input_path_tuple = get_basename_and_ext(input_path)
        last = len(input_path_tuple)-1
        output_path = f"{input_path_tuple[0]}.sorted{input_path_tuple[last]}"

    self.notify(f"Samtools sort {str(input_path)} to {str(output_path)} ...", title="Samtools sort")
    logging.info(f"Samtools sort {str(input_path)} to {str(output_path)} ...")
    self.query_one("#sam_sort_horizontal").add_class("running")
    
    await _run_samtools_sort(self, input_path, output_path)

async def _run_samtools_sort(self, input_file, output_file):
    try:
        samtools_sort_cmd = f"samtools sort {input_file} -o {output_file}"
        self.notify(samtools_sort_cmd, title="Samtools sort")
        logging.info("Running command: " + samtools_sort_cmd)
        
        # Run pysam.sort in a separate thread to avoid blocking the event loop
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(None, pysam.sort, "-o", output_file, input_file)
        
        logging.info(f"Samtools sort completed for {str(input_file)}")
        self.notify(f"Samtools sort completed for {str(input_file)}", title="Samtools sort")

    except SamtoolsError as e:
        logging.error(f"An error occurred during Samtools sort: {e}")
        self.notify(
            f"An error occurred during Samtools sort: {e}",
            severity="error",
            timeout=10.0,
            title="Samtools sort",
        )

    finally:
        self.query_one("#sam_sort_horizontal").remove_class("running")


def run_samtools_view(self):
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "sam_view_input_input")
    output_path = get_input(self, "sam_view_output_input")
    view_region = get_input(self, "sam_view_region_input")

    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return

    # WARN:
    # outputs output_file in workingdir
    # change to add working dir or find a way to get input prefix before output_file name
    if not output_path:
        input_path_tuple = get_basename_and_ext(input_path)
        last = len(input_path_tuple)-1
        output_path = f"{input_path_tuple[0]}.view{input_path_tuple[last]}"

    self.notify(f"Samtools view {str(input_path)} to {str(input_path)} ...", title="Samtools view")
    logging.info(f"Samtools view {str(input_path)} to {str(input_path)} ...")
    self.query_one("#sam_view_horizontal").add_class("running")
    threading.Thread(
        target=_run_samtools_view,
        args=(
            self,
            input_path,
            output_path,
            view_region,
        ),
    ).start()



def _run_samtools_view(self, input_file, output_file, view_region):

    try:
        samtools_view_cmd = f"samtools view {input_file} -o {output_file} region={view_region}"
        self.notify(samtools_view_cmd, title="Samtools view")
        logging.info("Running command: " + samtools_view_cmd)

        with pysam.AlignmentFile(input_file, "rb") as infile, pysam.AlignmentFile(
            output_file, "wb", header=infile.header
        ) as outfile:
            for read in infile.fetch(region=view_region) if view_region else infile:
                outfile.write(read)

        logging.info(f"Samtools view completed for {str(input_file)}")
        self.notify(f"Samtools view completed for {str(input_file)}", title="Samtools view")

    except SamtoolsError as e:
        logging.error(f"An error occurred during Samtools view: {e}")
        self.notify(
            f"An error occurred during Samtools view: {e}",
            severity="error",
            timeout=10.0,
            title="Samtools view",
        )

    finally:
        self.query_one("#sam_view_horizontal").remove_class("running")


def run_samtools_index(self):
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "sam_index_input_input")

    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools index")
        return

    self.notify(f"Samtools index {str(input_path)} ...", title="Samtools index")
    logging.info(f"Samtools index {str(input_path)} ...")
    self.query_one("#sam_index_horizontal").add_class("running")
    threading.Thread(
        target=_run_samtools_index,
        args=(
            self,
            input_path,
        ),
    ).start()

def _run_samtools_index(self, input_file):


    try:
        samtools_index_cmd = f"samtools index {input_file}"
        self.notify(samtools_index_cmd, title="Samtools index")
        logging.info("Running command: " + samtools_index_cmd)

        pysam.index(input_file)

        logging.info(f"Samtools index completed for {str(input_file)}")
        self.notify(f"Samtools index completed for {str(input_file)}", title="Samtools index")

    except SamtoolsError as e:
        logging.error(f"An error occurred during Samtools index: {e}")
        self.notify(
            f"An error occurred during Samtools index: {e}",
            severity="error",
            timeout=10.0,
            title="Samtools index",
        )

    finally:
        self.query_one("#sam_index_horizontal").remove_class("running")



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

