from textual.widgets import Label, LoadingIndicator, Button, Input, Select, Static
from textual.containers import Container, Horizontal
from modules.exec_func.hometab_funcs import FileSuggester
from textual.validation import Number


class HomeWidgets(Static):
    selection_list = [
        "reference",
        "reads",
        "vcf",
        "bed",
        "gtf",
        "gff",
    ]

    def compose(self):
        with Horizontal(id="InputProjectName"):
            yield Label("Enter project name (Enter absolute path)", id="Project_name")
            yield Input(id="Input_Project_name")
            yield Label("This will be your working directory", id="Project_location")
        with Container(id="Download_widget"):
            yield Label("Download (uses curl)", id="download_title")
            with Horizontal(id="Download_inputs"):
                yield Input(placeholder="Input url", id="Input_url")
                yield Input(
                    placeholder="outputfile name",
                    id="Input_outputfile_name",
                    suggester=FileSuggester(use_cache=False, case_sensitive=True),
                )
            with Horizontal(id="Download_options"):
                yield Select.from_values(
                    self.selection_list,
                    prompt="The downloaded file is for",
                    id="Select_outputdir",
                )
                yield Label("Saved in: ", id="output_path_label")
                yield Button("Download", id="Download_button", classes="action_buttons")
                yield LoadingIndicator(id="Download_loading")
        with Container(id="FastQC_widget"):
            yield Label("FastQC", id="FastQC_title")
            yield Input(
                placeholder="Input file name",
                id="FastQC_Input",
                suggester=FileSuggester(use_cache=False, case_sensitive=True),
            )
            with Horizontal(id="FastQC_Horizontal"):
                yield Label(
                    "FastQCs the files in the workingdir/data/reads directory if the input field is left empty",
                    id="bwa_description",
                )
                yield Button("Fastqc", id="FastQC_Button", classes="action_buttons")
                yield Button("View Results", id="view_FastQC_res")
                yield LoadingIndicator(id="FastQC_Loading")
        with Container(id="MultiQC_widget"):
            yield Label("MultiQC", id="MultiQC_title")
            yield Input(
                placeholder="Input file name",
                id="MultiQC_Input",
                suggester=FileSuggester(use_cache=False, case_sensitive=True),
            )
            with Horizontal(id="MultiQC_Horizontal"):
                yield Label(
                    "MultiQCs the files in the workingdir/results/fastqc directory if the input field is left empty \nNOTE: MultiQC can be used to visualise the output of various tools like samtools flagstat.\nFor more info check this",
                    id="bwa_description",
                )
                yield Button("MultiQC", id="MultiQC_Button", classes="action_buttons")
                yield Button("View Results", id="view_MultiQC_res")
                yield LoadingIndicator(id="MultiQC_Loading")
        with Container(id="bwa_widget"):
            yield Label("BWA (alignment uses mem)", id="bwa_title")
            with Horizontal(id="bwa_input_horizontal"):
                yield Input(
                    placeholder="Input reference genome file name",
                    id="bwa_ref_Input",
                    suggester=FileSuggester(use_cache=False, case_sensitive=True),
                )
                yield Input(
                    placeholder="No. of threads (default=4)",
                    id="bwa_threads_Input",
                    validators=[Number(minimum=1, maximum=200)],
                )
            yield Input(
                placeholder="Input reads file name (2 files)",
                id="bwa_reads_Input",
                suggester=FileSuggester(use_cache=False, case_sensitive=True),
            )

            with Horizontal(id="bwa_Horizontal"):
                yield Label(
                    "NOTE: The reference genome must be inedx before mapping\nIndexes all files in workingdir/data/reference directory if the input field is left empty\nMake sure to input two reads (1 sample)",
                    id="bwa_description",
                )
                yield Button("Index", id="bwa_index_Button", classes="action_buttons")
                yield Button("Align", id="bwa_align_Button", classes="action_buttons")
                yield Button("View Results", id="view_bwa_res")
                yield LoadingIndicator(id="bwa_Loading")
