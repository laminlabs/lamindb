import sys
from loguru import logger

# all times in UTC
logger.configure(
    handlers=[
        dict(
            sink=sys.stdout,
            format="{time:YYYY-MM-DD HH:mm:ss!UTC} | {message}",
        ),
    ],
)
