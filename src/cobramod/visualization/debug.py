"""
Debugging configuration for sub-package visualization.

Configures the debug logging tool. The format follows the syntax: `(asctime)
(levelname) (message)`. The logs are saved in the logs directory and includes
'visualization' in the file and its date. The default level is set to DEBUG.

The name of logger variable is `visualization`
"""

import datetime as dt
import logging
from pathlib import Path


debug_log = logging.getLogger("visualization")
debug_log.setLevel(logging.DEBUG)

FORMAT_STR = "[%(asctime)s] %(levelname)s %(message)s"
TIME_STR = "%H:%M:%S"

# Format
debug_formatter = logging.Formatter(FORMAT_STR, TIME_STR)

log_dir = Path("logs").absolute()
log_dir.mkdir(exist_ok=True)

log_path = log_dir.joinpath(
    f"visualization_{dt.date.today().strftime('%Y%m%d')}.log"
)
debug_handler = logging.FileHandler(log_path, mode="a+")
debug_handler.setFormatter(debug_formatter)

debug_log.addHandler(debug_handler)
