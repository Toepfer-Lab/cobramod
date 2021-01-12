#!/usr/bin/env python3
# Creating corresponding log for the visualization sub module
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
