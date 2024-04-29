"""
Database versioning configurator

This modules includes the class DataVersionConfigurator which is in charge of
obtaining and comparing the version from the obtained metabolic data.
It uses the same structure of using a Singleton for the configuration in
COBRApy
"""

from pathlib import Path
from typing import Optional, Union

import pandas as pd
from cobra.core.singleton import Singleton

import cobramod.error as cmod_error
from cobramod.debug import debug_log


class DataVersionConfigurator(metaclass=Singleton):
    def __init__(self):
        self.ignore_db_versions: bool = False
        self.database_version: Optional[pd.DataFrame] = None
        self.force_same_version: bool = False

    def get_database_version(self, directory: Union[str, Path]) -> pd.DataFrame:
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

        if self.database_version is not None:
            return self.database_version
        try:
            file = directory / "DatabaseVersions.csv"
            self.database_version = pd.read_csv(file)
        except FileNotFoundError:
            self.database_version = pd.DataFrame(columns=["orgid", "version"])

        return self.database_version

    def get_local_databases(self, directory: Union[Path, str]) -> pd.Series:
        if isinstance(directory, str):
            directory = Path(directory).absolute()
        databases = self.get_database_version(directory)
        databases = databases["orgid"]

        return databases

    def check_database_version(
        self, directory: Union[str, Path], database: str, version: Optional[str]
    ):
        """
        Function to compare the saved database version with the one
        of the retrieved data.

        Args:
            directory (Path): The folder used for storing data.
            database (str): Identifier of the database.
            version (str): The version of the database.
        """
        if isinstance(directory, str):
            directory = Path(directory).absolute()

        if self.database_version is None:
            self.get_database_version(directory=directory)

        assert isinstance(self.database_version, pd.DataFrame)
        row = self.database_version.loc[
            self.database_version["orgid"] == database
        ]

        if row.shape[0] == 0:
            return self.set_database_version(directory, database, version)

        if row.shape[0] > 1:
            exit()

        # NOTE: not all version use the same versioning
        expected_version = str(row["version"].tolist()[0])

        # row can not be more than one value but python does not know that
        # => any

        if expected_version != version:
            msg = (
                "Versions of {} do not match. Remote has version {} "
                "and local version is {}.".format(
                    database, version, expected_version
                )
            )
            debug_log.warning(msg)

            if self.ignore_db_versions:
                return

            if self.force_same_version:
                raise Exception(
                    "DB Configuration: force_same_version is True\n"
                    "Database differs from local dabatase and retrieved data!"
                )
            else:
                while True:
                    choice = input("Ignore version mismatch? (Y)es (N)o: ")

                    if choice.lower() == "y":
                        self.ignore_db_versions = True
                        break

                    if choice.lower() == "n":
                        raise cmod_error.UserInterruption(
                            "Interrupted by user input"
                        )

    def set_database_version(
        self, directory: Union[str, Path], database: str, version: Optional[str]
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
        if isinstance(directory, str):
            directory = Path(directory).absolute()

        if self.database_version is None:
            self.get_database_version(directory=directory)

        assert isinstance(self.database_version, pd.DataFrame)
        self.database_version = pd.concat(
            [
                self.database_version,
                pd.DataFrame({"orgid": [database], "version": [version]}),
            ],
            ignore_index=True,
        )

        self.database_version.to_csv(
            path_or_buf=str(directory / "DatabaseVersions.csv"), index=False
        )

        return True
