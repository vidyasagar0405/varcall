from textual.app import ComposeResult
from textual.containers import Container, Horizontal
from textual.validation import Number
from textual.widgets import Button, Input, Label, LoadingIndicator, Select, Static

from modules.exec.common  import file_suggester


class HomeWidgets(Static):
    selection_list = [
        "reference",
        "reads",
        "vcf",
        "bed",
        "gtf",
        "gff",
    ]

    def compose(self) -> ComposeResult:
        with Horizontal(id="InputProjectName"):
            yield Label( "Enter project name (Absolute path is preferred)", id="Project_name")
            yield Input(id="Input_Project_name")
            yield Label( "[bold #01B0DC]This will be your working directory & is necessary to run all the processes[/bold #01B0DC]", id="Project_location")

        with Container(id="Download_widget"):
            yield Label("Download (uses curl)", id="download_title")

            with Horizontal(id="Download_inputs"):
                yield Input(placeholder="Input url", id="Input_url")
                yield Input( placeholder="outputfile name", id="Input_outputfile_name", suggester=file_suggester)

            with Horizontal(id="Download_options"):
                yield Select.from_values( self.selection_list, prompt="The downloaded file is for", id="Select_outputdir")
                yield Label("Saved in: ", id="output_path_label")
                yield Button("Download", id="Download_button", classes="action_buttons")
                yield LoadingIndicator(id="Download_loading")

        with Container(id="FastQC_widget"):
            yield Label("FastQC", id="FastQC_title")
            with Horizontal(id="MultiQC_horizontal"):
                yield Input( placeholder="Input file name", id="FastQC_input_input", suggester=file_suggester)
                yield Input( placeholder="Output file name", id="FastQC_output_input", suggester=file_suggester)

            with Horizontal(id="FastQC_Horizontal"):
                yield Label( "FastQCs the files in the workingdir/data/reads directory [u]if the input field is left empty[/u]", id="bwa_description")
                yield Button("Fastqc", id="FastQC_Button", classes="action_buttons")
                yield Button("View Results", id="view_FastQC_res")
                yield LoadingIndicator(id="FastQC_Loading")

        with Container(id="MultiQC_widget"):
            yield Label("MultiQC", id="MultiQC_title")
            with Horizontal(id="MultiQC_horizontal"):
                yield Input( placeholder="Input file name", id="MultiQC_input_input", suggester=file_suggester)
                yield Input( placeholder="Output file name", id="MultiQC_output_input", suggester=file_suggester)

            with Horizontal(id="MultiQC_Horizontal"):
                yield Label(
                    "MultiQCs the files in the workingdir/results/fastqc directory [u]if the input field is left empty[/u] \n[bold red]NOTE[/bold red]: MultiQC can be used to visualise the output of various tools like samtools flagstat.\nFor more info check this",
                    id="bwa_description",
                )
                yield Button("MultiQC", id="MultiQC_Button", classes="action_buttons")
                yield Button("View Results", id="view_MultiQC_res")
                yield LoadingIndicator(id="MultiQC_Loading")

        with Container(id="bwa_widget"):
            yield Label("BWA (alignment uses mem)", id="bwa_title")

            with Horizontal(id="bwa_input_horizontal_1"):
                yield Input( placeholder="Input reference genome file name", id="bwa_ref_input", suggester=file_suggester)
                yield Input( placeholder="No. of threads (default=4)", id="bwa_threads_input", validators=[Number(minimum=1, maximum=500)])

            with Horizontal(id="bwa_input_horizontal_2"):
                yield Input( placeholder="Input read file", id="bwa_reads_input_1", suggester=file_suggester)
                yield Input( placeholder="Input read file", id="bwa_reads_input_2", suggester=file_suggester)

            yield Input( placeholder="Output file name", id="bwa_output_input", suggester=file_suggester)

            with Horizontal(id="bwa_Horizontal"):
                yield Label(
                    "[bold red]NOTE[/bold red]: The reference genome must be inedx before mapping\nIndexes all files in workingdir/data/reference directory if the input field is left empty\nMake sure to input two reads (1 sample)",
                    id="bwa_description",
                )
                yield Button("Index", id="bwa_index_Button", classes="action_buttons")
                yield Button("Align", id="bwa_align_Button", classes="action_buttons")
                yield Button("View Results", id="view_bwa_res")
                yield LoadingIndicator(id="bwa_Loading")
