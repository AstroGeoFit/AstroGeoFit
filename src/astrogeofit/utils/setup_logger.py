"""
Copyright (C) 2025  CNRS
Lead author : J. Laskar
    
This file is part of the AstroGeoFit software.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import logging


def setup_logger(logs: bool):
    logger = logging.getLogger("ToolLogger")

    # Avoid adding handlers multiple times
    if not logger.hasHandlers():
        # Set the minimum log level
        logger.setLevel(logging.DEBUG)

        # Formatter for log messages
        formatter = logging.Formatter("%(levelname)s - %(message)s")

        # Console Handler (for printing to the console)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_handler.setFormatter(formatter)
        # Add handlers to the logger
        logger.addHandler(console_handler)
        if logs:
            # File Handler (for saving logs to a file)
            file_handler = logging.FileHandler("tool.log")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
    return logger
