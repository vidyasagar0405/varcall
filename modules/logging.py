import logging


def setup_logging():
    logging.basicConfig(
        filename="master_log.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
