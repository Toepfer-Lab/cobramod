#!/usr/bin/env python3
import logging

# Creating corresponding Logs
# Format
debug_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# Handler
debug_handler = logging.FileHandler("debug.log", mode="a+")
debug_handler.setFormatter(debug_formatter)
# Log
debug_log = logging.getLogger("debug_log")
debug_log.setLevel(logging.DEBUG)
# GenLog.ad
debug_log.addHandler(debug_handler)
