#!/usr/bin/env python3
"""Debugging configuration for CobraMod

This module configures the debug logging tool. The format of the syntax:
`(asctime) (levelname) (message)`.

The FileHandler is set to 'debug.log' with the mode "+a".
The default level is set to INFO.

The name of logger variable is `debug_log`
"""
from logging import Formatter, FileHandler, getLogger, INFO

# Creating corresponding Logs
# Format
debug_formatter = Formatter("%(asctime)s %(levelname)s %(message)s")
# Handler
debug_handler = FileHandler("debug.log", mode="a+")
debug_handler.setFormatter(debug_formatter)
# Log
debug_log = getLogger("debug_log")
debug_log.setLevel(INFO)
# GenLog.ad
debug_log.addHandler(debug_handler)
