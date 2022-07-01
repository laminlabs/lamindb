import sys

from loguru import logger

logger.configure(
    handlers=[
        dict(
            sink=sys.stdout,
            format="{message}",
        ),
    ],
)
