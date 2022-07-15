import sys

from loguru import logger

# ANSI color code: https://gist.github.com/iansan5653/c4a0b9f5c30d74258c5f132084b78db9
ANSI_COLORS = dict(
    bold="\x1b[1m",
    green="\x1b[1;92m",
    blue="\x1b[1;94m",
    purple="\x1b[1;95m",
    yellow="\x1b[1;93m",
    reset="\x1b[0m",
)

default_handler = dict(
    sink=sys.stdout,
    format="{level.icon} | {message}",
)

logger.configure(handlers=[default_handler])
logger.level("INGEST", no=15, icon="✅")
