#!/usr/bin/env python3
"""Debugging configuration for sub-package visualization.

Configures the debug logging tool. The format follows the syntax:
`(asctime) (levelname) (message)`.

The FileHandler is set to 'visualization.log' with the mode "+a".
The default level is set to DEBUG.

The name of logger variable is `visualization`

"""
from logging import Formatter, FileHandler, getLogger, DEBUG

# Format
debug_formatter = Formatter("%(asctime)s %(levelname)s %(message)s")
# Handler
debug_handler = FileHandler("visualization.log", mode="a+")
debug_handler.setFormatter(debug_formatter)
# Log
debug_log = getLogger("visualization")
debug_log.setLevel(DEBUG)
# GenLog.ad
debug_log.addHandler(debug_handler)
