#!/usr/bin/env python3
"""Abstract class Base

This module contains an abstract class, which is used as a base for
the retrieval of all the data from different databases. The class
:class:`cobramod.parsing.base.BaseParser` defines important methods such as:

- _retrieve_data: Retrieve from database or locally.
- _parse: Pase information and return it as a dictionary.
- _return_database: Check method.
- _read_file: Read file and return the information.
"""
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Union


class BaseParser(ABC):
    @staticmethod
    def _check_transport(data_dict: dict) -> bool:
        """
        Check if given reaction includes same participant in both sides,
        meaning that there is a tranport reaction. Prefix is not taken into
        consideration.
        """
        sequence = [item[2:] for item in data_dict.keys()]
        return len(sequence) != len(set(sequence))

    @staticmethod
    @abstractmethod
    def _retrieve_data(
        directory: Path,
        identifier: str,
        database: str,
        debug_level: int,
        **kwargs
    ) -> dict:
        """
        Method to retrieve data from server or local directory.
        """
        pass

    @staticmethod
    @abstractmethod
    def _parse(root: Union[Any, dict], directory: Path) -> dict:
        """
        Basic method to parse information and return it as dictionary.
        """
        pass

    @staticmethod
    @abstractmethod
    def _check_database(database: str):
        """
        Basic check method.
        """
        pass

    @staticmethod
    @abstractmethod
    def _read_file(filename: Path) -> Any:
        """
        Basic method to read file.
        """
        pass
