import os
import logging
from pathlib import Path
import tempfile
import sys
import traceback

def setup_logging(log_file=None, log_level=logging.INFO):
    """
    Set up logging in an OS-independent way.

    Args:
        log_file (str or Path, optional): Custom log file path. If None, uses temp directory.
        log_level (int, optional): Logging level. Defaults to logging.INFO.
    """
    try:
        # Print to stderr to avoid disrupting TUI
        stderr_print = lambda msg: print(msg, file=sys.stderr)

        if log_file is None:
            # Specify absolute path in home directory instead of temp
            home_dir = Path.home()
            log_dir = home_dir / ".varcall"
            log_dir.mkdir(exist_ok=True)
            log_file = log_dir / "varcall.log"
            stderr_print(f"Using log file at: {log_file}")
        else:
            log_file = Path(log_file)
            stderr_print(f"Using custom log file at: {log_file}")

        # Ensure parent directory exists
        log_file.parent.mkdir(parents=True, exist_ok=True)

        # Make sure we can write to the file
        try:
            # Open the file directly to ensure it's writable
            with open(log_file, 'a') as f:
                pass
        except Exception as e:
            stderr_print(f"Cannot write to log file {log_file}: {e}")
            # Try alternative location as fallback
            log_file = Path.home() / "varcall.log"
            stderr_print(f"Trying alternative log location: {log_file}")

        # Configure logging - file only, no console handler
        # Clear any existing handlers
        root_logger = logging.getLogger()
        root_logger.handlers = []  # Remove all handlers
        root_logger.setLevel(log_level)

        # File handler for writing to log file
        try:
            handler = logging.FileHandler(str(log_file))
            formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
            handler.setFormatter(formatter)
            root_logger.addHandler(handler)

            # Force a write to the log file to verify it works
            logging.info(f"Logging initialized. Log file: {log_file}")
            stderr_print(f"Logging successfully configured to write to: {log_file}")
            return True
        except Exception as e:
            stderr_print(f"ERROR: Could not create log file handler: {e}")
            stderr_print(f"Error details: {traceback.format_exc()}")
            return False

    except Exception as e:
        stderr_print(f"CRITICAL ERROR in setup_logging: {e}")
        stderr_print(f"Error details: {traceback.format_exc()}")
        return False
