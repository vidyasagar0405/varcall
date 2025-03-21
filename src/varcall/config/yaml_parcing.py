from pathlib import Path
import yaml
from varcall.process.process_class import ProcessConfig
from varcall.config.process_dict import DEFAULT_MASTER_CONFIG
import os
import logging

def compile_master_config(yaml_file: Path | None = None) -> dict[str, dict[str, ProcessConfig]]:
    """
    Reads a YAML configuration file and converts its content into a nested dictionary
    of ProcessConfig instances, structured as:
    {
        "home": { "process_name": ProcessConfig, ... },
        "samtools": { "process_name": ProcessConfig, ... },
        "bcftools": { "process_name": ProcessConfig, ... },
        "gatk": { "process_name": ProcessConfig, ... },
        "pipeline": { "process_name": ProcessConfig, ... },
    }

    Args:
        yaml_file (Path, optional): Path to the YAML file containing the configuration.
                                   If None, attempts to load from default locations.

    Returns:
        dict[str, dict[str, ProcessConfig]]: The compiled master configuration.
    """
    # Set up logging
    logger = logging.getLogger(__name__)

    # Define potential config file locations in order of preference
    if yaml_file is None:
        possible_paths = [
            Path("~/.config/varcall/config.yaml").expanduser(),
            Path(os.environ.get("VARCALL_CONFIG", "~/.varcall/config.yaml")).expanduser(),
            Path("./config.yaml"),
            Path("./varcall_config.yaml"),
        ]

        # Find the first existing config file
        for path in possible_paths:
            if path.exists():
                yaml_file = path
                logger.info(f"Using configuration file: {yaml_file}")
                break

        # If no config file found, use default or raise error
        if yaml_file is None:
            logger.warning("No configuration file found. Using default configuration.")
            return create_default_config()

    # Ensure the file exists
    if not yaml_file.exists():
        logger.error(f"Configuration file not found: {yaml_file}")
        logger.info("Using default configuration.")
        return create_default_config()

    # Load YAML configuration file into a dictionary
    try:
        with open(yaml_file, "r") as file:
            config_yaml = yaml.safe_load(file)

        if not config_yaml:
            logger.warning(f"Empty or invalid configuration file: {yaml_file}")
            return create_default_config()

    except (yaml.YAMLError, IOError) as e:
        logger.error(f"Error reading configuration file {yaml_file}: {e}")
        return create_default_config()

    master_config = {}
    # Loop over each tab and its processes
    for tab_name, processes in config_yaml.items():
        tab_dict = {}
        if not isinstance(processes, dict):
            logger.warning(f"Invalid format for tab '{tab_name}': expected dict but got {type(processes)}. Skipping.")
            continue

        for process_name, process_config_dict in processes.items():
            try:
                # Convert each process config dictionary into a ProcessConfig instance
                process_config = ProcessConfig(**process_config_dict)
                tab_dict[process_name] = process_config
            except Exception as e:
                logger.error(f"Error creating ProcessConfig for '{process_name}' in tab '{tab_name}': {e}")
                # Consider adding default ProcessConfig for this entry

        if tab_dict:  # Only add the tab if it has valid processes
            master_config[tab_name] = tab_dict

    # Ensure master_config is not empty
    if not master_config:
        logger.warning("No valid configuration found. Using default configuration.")
        return create_default_config()

    return master_config

def create_default_config() -> dict[str, dict[str, ProcessConfig]]:
    """
    Creates a default configuration when no valid configuration file is found.

    Returns:
        dict[str, dict[str, ProcessConfig]]: Default configuration structure.
    """
    # Define minimal default configuration
    default_config = DEFAULT_MASTER_CONFIG

    return default_config

# Compile the master configuration from the YAML file
MASTER_CONFIG: dict[str, dict[str, ProcessConfig]] = compile_master_config()
