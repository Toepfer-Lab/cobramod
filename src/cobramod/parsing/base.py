#!/usr/bin/env python3
from pathlib import Path
from abc import ABC, abstractmethod
from typing import Any, Union


class BaseParser(ABC):
    @staticmethod
    def _define_base_dir(directory: Path, database: str) -> Path:
        """
        Returns Path object for given database. If directory does not exist.
        It will be created.

        Args:
            directory (Path): Parent directory.
            database (str): Name of database. Options: "META", "ARA". "KEGG"

        Returns:
            Path: Path object for database.
        """
        if directory.joinpath(database).exists():
            return directory.joinpath(database)
        else:
            directory.joinpath(database).mkdir()
            return directory.joinpath(database)

    @staticmethod
    def _check_transport(data_dict: dict) -> bool:
        """
        check if data_dict includes repetitions, meaning that there is a
        tranport reactions. Prefix is not taken into consideration.
        """
        sequence = [item[2:] for item in data_dict.keys()]
        return len(sequence) != len(set(sequence))

    @staticmethod
    @abstractmethod
    def _retrieve_data(
        directory: Path, identifier: str, database: str, debug_level: int
    ) -> dict:
        pass

    @staticmethod
    @abstractmethod
    def _parse(root: Union[Any, str]) -> dict:
        pass

    @staticmethod
    @abstractmethod
    def _return_database(database: str) -> str:
        pass
