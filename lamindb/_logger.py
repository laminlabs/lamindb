import sys

from loguru import logger

default_handler = dict(
    sink=sys.stdout,
    format="{message}",
)

logger.configure(handlers=[default_handler])
