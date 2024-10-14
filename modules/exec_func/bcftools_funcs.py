import logging
import subprocess
import threading

from modules.exec_func.common_funcs import get_input
from modules.logging import setup_logging

setup_logging()

# 1. bcftools mpileup command
def run_bcf_mpileup(self):
    """

    """
    if not self.workingDir:
        self.workingDir = __file__
    ref = get_input(self, "bcf_mpileup_reference_input")
    input_path = get_input(self, "bcf_mpileup_input_input")
    output_path = get_input(self, "bcf_mpileup_output_input")
    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return
    if not output_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return
    if not ref:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return

    self.notify(f"bcf mpileup {str(input_path)}...", title="bcf mpileup")
    logging.info(f"bcf mpileup {str(input_path)}...")
    self.query_one("#bcf_mpileup_horizontal").add_class("running")
    threading.Thread(
        target=_run_bcf_mpileup,
        args=(
            self,
            ref,
            input_path,
            output_path,
        ),
    ).start()


def _run_bcf_mpileup(self, ref: str, input_path: str, output_path: str) -> None:
    """

    """
    try:
        bcf_mpileup_cmd  = f"bcftools mpileup -f {ref} {input_path} -o {output_path}"
        self.notify(bcf_mpileup_cmd, title="bcf mpileup")
        logging.info("Running command: " + bcf_mpileup_cmd)
        subprocess.run(
            [bcf_mpileup_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"bcf mpileup completed for {str(input_path)}")
        self.notify(f"bcf mpileup completed for {str(input_path)}", title="bcf mpileup")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during bcf mpileup: {e.stderr}")
        self.notify(
            f"An error occurred during FastQc: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="bcf mpileup",
        )

    finally:
        self.query_one("#bcf_mpileup_horizontal").remove_class("running")


# 2. bcftools call command

def run_bcf_call(self):
    """

    """
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "bcf_call_input_input")
    output_path = get_input(self, "bcf_call_output_input")
    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return
    if not output_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return

    self.notify(f"bcf call {str(input_path)}...", title="bcf call")
    logging.info(f"bcf call {str(input_path)}...")
    self.query_one("#bcf_call_horizontal").add_class("running")
    threading.Thread(
        target=_run_bcf_call,
        args=(
            self,
            input_path,
            output_path,
        ),
    ).start()


def _run_bcf_call(self, input_path: str, output_path: str) -> None:
    """

    """
    try:
        bcf_call_cmd = f"bcftools call -mv -Oz -o {output_path} {input_path}"
        self.notify(bcf_call_cmd, title="bcf call")
        logging.info("Running command: " + bcf_call_cmd)
        subprocess.run(
            [bcf_call_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"bcf call completed for {str(input_path)}")
        self.notify(f"bcf call completed for {str(input_path)}", title="bcf call")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during bcf call: {e.stderr}")
        self.notify(
            f"An error occurred during FastQc: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="bcf call",
        )

    finally:
        self.query_one("#bcf_call_horizontal").remove_class("running")

# 3. bcftools filter command
def run_bcf_filter(self):
    """

    """
    if not self.workingDir:
        self.workingDir = __file__
    filter = get_input(self, "bcf_filter_filter_input")
    input_path = get_input(self, "bcf_filter_input_input")
    output_path = get_input(self, "bcf_filter_output_input")
    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return
    if not output_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return

    self.notify(f"bcf filter {str(input_path)}...", title="bcf filter")
    logging.info(f"bcf filter {str(input_path)}...")
    self.query_one("#bcf_filter_horizontal").add_class("running")
    threading.Thread(
        target=_run_bcf_filter,
        args=(
            self,
            filter,
            input_path,
            output_path,
        ),
    ).start()


def _run_bcf_filter(self,filter: str, input_path: str, output_path: str) -> None:
    """

    """
    try:
        bcf_filter_cmd = f"bcftools filter -s LOWQUAL -e '{filter}' {input_path} -o {output_path}"
        self.notify(bcf_filter_cmd, title="bcf filter")
        logging.info("Running command: " + bcf_filter_cmd)
        subprocess.run(
            [bcf_filter_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"bcf filter completed for {str(input_path)}")
        self.notify(f"bcf filter completed for {str(input_path)}", title="bcf filter")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during bcf filter: {e.stderr}")
        self.notify(
            f"An error occurred during FastQc: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="bcf filter",
        )

    finally:
        self.query_one("#bcf_filter_horizontal").remove_class("running")



# 4. bcftools norm command

def run_bcf_norm(self):
    """

    """
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "bcf_norm_input_input")
    output_path = get_input(self, "bcf_norm_output_input")
    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return
    if not output_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return

    self.notify(f"bcf norm {str(input_path)}...", title="bcf norm")
    logging.info(f"bcf norm {str(input_path)}...")
    self.query_one("#bcf_norm_horizontal").add_class("running")
    threading.Thread(
        target=_run_bcf_norm,
        args=(
            self,
            input_path,
            output_path,
        ),
    ).start()


def _run_bcf_norm(self, input_path: str, output_path: str) -> None:
    """

    """
    try:
        bcf_norm_cmd = f"bcftools norm -f {input_path} {output_path}"
        self.notify(bcf_norm_cmd, title="bcf norm")
        logging.info("Running command: " + bcf_norm_cmd)
        subprocess.run(
            [bcf_norm_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"bcf norm completed for {str(input_path)}")
        self.notify(f"bcf norm completed for {str(input_path)}", title="bcf norm")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during bcf norm: {e.stderr}")
        self.notify(
            f"An error occurred during FastQc: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="bcf norm",
        )

    finally:
        self.query_one("#bcf_norm_horizontal").remove_class("running")


# 5. bcftools stats command

def run_bcf_stats(self):
    """

    """
    if not self.workingDir:
        self.workingDir = __file__
    input_path = get_input(self, "bcf_stats_input_input")
    output_path = get_input(self, "bcf_stats_output_input")
    if not input_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return
    if not output_path:
        self.notify("Please provide a valid path", severity="warning", title="Samtools view")
        return

    self.notify(f"bcf stats {str(input_path)}...", title="bcf stats")
    logging.info(f"bcf stats {str(input_path)}...")
    self.query_one("#bcf_stats_horizontal").add_class("running")
    threading.Thread(
        target=_run_bcf_stats,
        args=(
            self,
            input_path,
            output_path,
        ),
    ).start()


def _run_bcf_stats(self, input_path: str, output_path: str) -> None:
    """

    """
    try:
        bcf_stats_cmd = f"bcftools stats {input_path} > {output_path}"
        self.notify(bcf_stats_cmd, title="bcf stats")
        logging.info("Running command: " + bcf_stats_cmd)
        subprocess.run(
            [bcf_stats_cmd],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=True,
        )
        logging.info(f"bcf stats completed for {str(input_path)}")
        self.notify(f"bcf stats completed for {str(input_path)}", title="bcf stats")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred during bcf stats: {e.stderr}")
        self.notify(
            f"An error occurred during FastQc: {e.stderr}",
            severity="error",
            timeout=10.0,
            title="bcf stats",
        )

    finally:
        self.query_one("#bcf_stats_horizontal").remove_class("running")
