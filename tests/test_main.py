import os
import filecmp
import unittest
import pytest
from src.varcall.main import Varcall
from textual.pilot import Pilot

app = Varcall()

# WARN: Valid only when samtools is focused on startup
def test_samtoolstab(snap_compare):
    assert snap_compare("../src/main.py")

def test_hometab(snap_compare):
    assert snap_compare("../src/main.py", press=["tab", "left"])

def test_bcftools(snap_compare):
    assert snap_compare("../src/main.py", press=["tab", "right"])

def test_helptab(snap_compare):
    assert snap_compare("../src/main.py", press=["tab", "right", "right"])

def test_homebuttons(snap_compare):

    @pytest.mark.asyncio
    async def home_download_button(pilot: Pilot):
        await pilot.resize_terminal(168, 150)
        await pilot.press("tab", "left", "tab", "t", "r", "i", "a", "l", "enter")
        await pilot.click("#Download_button")
    assert snap_compare(Varcall(), run_before=home_download_button,  terminal_size=(168, 75))

    @pytest.mark.asyncio
    async def home_fastqc_button(pilot: Pilot):
        await pilot.resize_terminal(168, 75)
        await pilot.press("tab", "left", "tab", "t", "r", "i", "a", "l", "enter")
        await pilot.click("#FastQC_Button")
    assert snap_compare(Varcall(), run_before=home_fastqc_button,  terminal_size=(168, 75))

    @pytest.mark.asyncio
    async def home_multiqc_button(pilot: Pilot):
        await pilot.resize_terminal(168, 75)
        await pilot.press("tab", "left", "tab", "t", "r", "i", "a", "l", "enter")
        await pilot.click("#MultiQC_Button")
    assert snap_compare(Varcall(), run_before=home_multiqc_button,  terminal_size=(168, 75))

    @pytest.mark.asyncio
    async def home_bwa_index_button(pilot: Pilot):
        await pilot.resize_terminal(168, 75)
        await pilot.press("tab", "left", "tab", "t", "r", "i", "a", "l", "enter")
        await pilot.click("#bwa_index_Button")
    assert snap_compare(Varcall(), run_before=home_bwa_index_button,  terminal_size=(168, 75))

    @pytest.mark.asyncio
    async def home_bwa_mem_button(pilot: Pilot):
        await pilot.resize_terminal(168, 75)
        await pilot.press("tab", "left", "tab", "t", "r", "i", "a", "l", "enter")
        await pilot.click("#bwa_mem_Button")
    assert snap_compare(Varcall(), run_before=home_bwa_mem_button,  terminal_size=(168, 75))

def test_pipeline_buttons(snap_compare):

    @pytest.mark.asyncio
    async def pipeline_button(pilot: Pilot):
        await pilot.resize_terminal(168, 40)
        await pilot.press("tab", "right", "right", "end")
        await pilot.click("#start_pipeline_button")
    assert snap_compare(Varcall(), run_before=pipeline_button,  terminal_size=(168, 40))

class TestFindMatchingFiles(unittest.TestCase):
    def compare_file(self):
        test_files_dir = "../trial/"
        template_files_dir = "./test_trial"
        _, mismatch, _ = filecmp.cmpfiles(test_files_dir, template_files_dir, shallow=False, common=["t"])
        assert mismatch == []

if __name__ == "__main__":
    test_samtoolstab(Varcall)
    test_hometab(Varcall)
