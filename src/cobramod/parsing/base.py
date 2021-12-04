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
from typing import Any, Union, Optional

import pandas as pd
from warnings import warn

from cobramod.error import UserInterruption


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

    ignore_db_versions: bool = False
    database_version: Optional[pd.DataFrame] = None

    @staticmethod
    def _get_database_version(directory: Path) -> pd.DataFrame:
        """
        Loads and returns the database versioning file.

        Args:
            directory (Path): The folder used for storing data.

        Returns:
            (pd.DataFrame): A DataFrame containing the orgid and version, at
                the time of the first retrieval, for all databases used so far.
        """
        if isinstance(directory, str):
            directory = Path(directory)

        if BaseParser.database_version is not None:
            return BaseParser.database_version
        try:
            file = directory / "DatabaseVersions.csv"
            BaseParser.database_version = pd.read_csv(file)
        except FileNotFoundError:
            BaseParser.database_version = pd.DataFrame(
                columns=["orgid", "version"]
            )

        return BaseParser.database_version

    @staticmethod
    def _get_local_databases(directory: Path) -> pd.Series:
        databases = BaseParser._get_database_version(directory)
        databases = databases["orgid"]

        return databases

    @staticmethod
    def check_database_version(directory: Path, database: str, version: str):
        """
        Function to compare the saved database version with the one
        of the retrieved data.

        Args:
            directory (Path): The folder used for storing data.
            database (str): Identifier of the database.
            version (str): The version of the database.
        """
        if BaseParser.database_version is None:
            BaseParser._get_database_version(directory=directory)

        assert isinstance(BaseParser.database_version, pd.DataFrame)
        row = BaseParser.database_version.loc[
            BaseParser.database_version["orgid"] == database
        ]

        if row.shape[0] == 0:
            return BaseParser._set_database_version(
                directory, database, version
            )

        if row.shape[0] > 1:
            exit()

        expected_version = row["version"].tolist()[0]

        # row can not be more than one value but python does not know that
        # => any

        if expected_version != version:
            msg = (
                "Versions of {} do not match. Remote has version {} "
                "and local version is {}.".format(
                    database, version, expected_version
                )
            )
            warn(
                message=msg,
                category=UserWarning,
            )

            if BaseParser.ignore_db_versions:
                return

            while True:
                choice = input("Ignore version mismatch? (Y)es (N)o")

                if choice.lower() == "y":
                    break

                if choice.lower() == "n":
                    raise UserInterruption("Interrupted by user input")

    @staticmethod
    def _set_database_version(
        directory: Path, database: str, version: str
    ) -> bool:
        """
        Adds the version of a database to the local data versioning file.

        Args:
            directory (Path): The folder used for storing data.
            database (str): Identifier of the database.
            version (str): The version of the database.

        Returns:
             (bool): Returns True if the addition was successful.

        """
        if BaseParser.database_version is None:
            BaseParser._get_database_version(directory=directory)

        assert isinstance(BaseParser.database_version, pd.DataFrame)
        BaseParser.database_version = BaseParser.database_version.append(
            {"orgid": database, "version": version}, ignore_index=True
        )

        BaseParser.database_version.to_csv(
            path_or_buf=str(directory / "DatabaseVersions.csv"), index=False
        )

        return True

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
    def _check_database(directory: Path, database: str):
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
