"""
Debugging configuration for CobraMod

The logs are saved in the directory 'logs' created in the working directory.
Cobramod saves the logs by date, i.e, running the commands in different days
will results in different files. The default level is set to INFO.

The name of logger variable is `debug_log`
"""

import datetime as dt
import logging
from pathlib import Path

import colorlog

debug_log = logging.getLogger("debug_log")
debug_log.setLevel(logging.INFO)

format_str = "%(log_color)s%(message)s"
TIME_STR = "%H:%M:%S"

colors = {
    "DEBUG": "white",
    "INFO": "cyan",
    "WARNING": "yellow",
    "ERROR": "red",
    "CRITICAL": "purple",
}
formatter = colorlog.ColoredFormatter(format_str, TIME_STR, log_colors=colors)

# In direct stream
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
debug_log.addHandler(stream_handler)

# File logs
format_str = "[%(asctime)s] %(levelname)s %(message)s"
formatter_file = logging.Formatter(format_str, TIME_STR)

log_dir = Path("logs").absolute()
log_dir.mkdir(exist_ok=True)
log_path = log_dir.joinpath(
    f"cobramod_{dt.date.today().strftime('%Y%m%d')}.log"
)

debug_handler = logging.FileHandler(log_path, mode="a+")
debug_handler.setFormatter(formatter_file)

debug_log.addHandler(debug_handler)


def change_to_debug() -> None:
    """
    Changes the logging level to DEBUG and changes the formatting of it.
    This function is intented to be used only for the tests
    """
    debug_log.setLevel(logging.DEBUG)
    format_str = "[%(asctime)s] %(log_color)s%(levelname)s %(message)s"
    formatter = colorlog.ColoredFormatter(
        format_str, TIME_STR, log_colors=colors
    )
    stream_handler.setFormatter(formatter)
