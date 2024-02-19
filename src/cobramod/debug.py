"""Debugging configuration for CobraMod

This module configures the debug logging tool. The format of the syntax:
`(asctime) (levelname) (message)`.

The FileHandler is set to 'debug.log' with the mode "+a".
The default level is set to INFO.

The name of logger variable is `debug_log`
"""
import datetime as dt
import logging
from pathlib import Path

import colorlog

# Log
debug_log = logging.getLogger("debug_log")
debug_log.setLevel(logging.INFO)

format_str = "[%(asctime)s] %(log_color)s%(levelname)s %(message)s"
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

format_str = "[%(asctime)s] %(levelname)s %(message)s"
# File logs
formatter = logging.Formatter(format_str, TIME_STR)

log_dir = Path("logs").absolute()
log_dir.mkdir(exist_ok=True)
log_path = log_dir.joinpath(f"cobramod_{dt.date.today().strftime('%Y%m%d')}")

debug_handler = logging.FileHandler(log_path, mode="a+")
debug_handler.setFormatter(formatter)

debug_log.addHandler(debug_handler)
