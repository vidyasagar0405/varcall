import unittest
# from varcall.main import Varcall
from varcall.modules.exec.command import Command 
from varcall.modules.exec.utils import get_output_name

class TestCommand(unittest.TestCase):
    def test_download_cmd(self):
        # Varcall().run()
        cmd = Command(name="download", workingDir = "trial", first_path="https://example.com/file.txt", second_path="file.txt")
        commands = cmd.get_command()
        expected = "curl -L https://example.com/file.txt -o file.txt"
        self.assertEqual(commands, expected)

    def test_get_output_name(self):
        cmd = Command(name="fastqc", workingDir = "trial", first_path="", second_path="t")
        commands = cmd.get_command()
        expected = "fastqc trial/data/reads/* -o trial/results/fastqc/"
        self.assertEqual(commands, expected)


if __name__ == "__main__":
    unittest.main()
