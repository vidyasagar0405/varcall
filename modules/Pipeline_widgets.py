from textual.app import ComposeResult
from textual.widgets import Button, Label, Static, Input


class PipelineWidgets(Static):
    def compose(self) -> ComposeResult:
        yield Label(
            "It is recommended to do fastqc and make sure the reads are of sufficiant quality before starting the pipeline"
        )
        yield Input(placeholder="Enter project name", id="pipeline_project_name_input")
        yield Input(placeholder="Enter reference file", id="pipeline_ref_input")
        yield Input(placeholder="Enter a read file", id="pipeline_read1_input")
        yield Input(placeholder="Enter the other read file", id="pipeline_read2_input")
        yield Button("Start Pipeline")
