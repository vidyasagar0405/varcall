import threading
import subprocess
import logging
from typing import Dict, List, Any
from dataclasses import dataclass


# Import the Process and ProcessConfig classes
@dataclass
class ProcessConfig:
    name: str
    command: str
    input_fields: List[str]
    required_fields: List[str]
    description: str = ""
    success_message: str = ""
    error_message: str = ""


class Process:
    def __init__(self, app_instance: Any, config: ProcessConfig):
        self.app = app_instance
        self.config = config
        self.inputs: Dict[str, str] = {}

    def get_inputs(self, form_data) -> bool:
        """Collect all inputs from form data"""
        for field in self.config.input_fields:
            value = form_data.get(field, "")
            if not value and field in self.config.required_fields:
                return False, f"Please provide a value for {field}"
            self.inputs[field] = value or ""
        return True, ""

    def format_command(self) -> str:
        """Format command string with input values"""
        return self.config.command.format(**self.inputs)

    def run(self, form_data):
        """Main entry point to run the process"""
        success, message = self.get_inputs(form_data)
        if not success:
            return {"status": "error", "message": message}

        logging.info(f"Running {self.config.name}...")

        # Start process thread
        thread = threading.Thread(target=self._run_process, args=())
        thread.daemon = True
        thread.start()

        return {"status": "running", "message": f"Running {self.config.name}..."}

    def _run_process(self):
        """Execute the actual process in a thread"""
        try:
            cmd = self.format_command()
            logging.info(f"Running command: {cmd}")

            result = subprocess.run(
                cmd.split(),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            msg = (
                self.config.success_message
                or f"{self.config.name} completed successfully"
            )
            logging.info(msg)
            logging.info(f"Output: {result.stdout}")

            # Store the result in a file for later retrieval
            with open(f"results/{self.config.name.lower()}_result.txt", "w") as f:
                f.write(result.stdout)

        except subprocess.CalledProcessError as e:
            error_msg = (self.config.error_message or f"Error in {self.config.name}: {e.stderr}")
            logging.error(error_msg)

            # Store the error in a file for later retrieval
            with open(f"results/{self.config.name.lower()}_error.txt", "w") as f:
                f.write(e.stderr)
