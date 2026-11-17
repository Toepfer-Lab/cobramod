import logging
import os
from pathlib import Path
from threading import Lock
from typing import Optional

import platformdirs
import requests
from pyrate_limiter import Duration, Limiter, Rate
from requests.adapters import HTTPAdapter, Retry

import cobramod.utils as cmod_utils

logger = logging.getLogger("cobramod.Settings")
logger.propagate = True


class SingletonMeta(type):
    """
    Threadsafe Singleton metalass.
    """

    _instances = {}  # type: ignore[var-annotated]
    _lock: Lock = Lock()

    def __call__(cls, *args, **kwargs):
        with cls._lock:
            if cls not in cls._instances:
                instance = super().__call__(*args, **kwargs)
                cls._instances[cls] = instance
        return cls._instances[cls]


class Settings(metaclass=SingletonMeta):
    def __init__(self):
        # One session per database, e.g. "BIOCYC" or "KEGG", see _session
        self.__sessions: dict[str, requests.Session] = {}

        self.__set_biocyc_login_from_env()

        if self.__biocyc_name is None or self.__biocyc_password is None:
            self.__set_biocyc_login_from_file()

        self.__cache_dir: Path = Path(
            platformdirs.user_cache_dir(
                appname="cobramod",
                appauthor="cobramod",
                ensure_exists=True,
            )
        )

        self.limiter = Limiter(
            Rate(1, Duration.SECOND),
        )

    __biocyc_password: Optional[str] = None
    __biocyc_name: Optional[str] = None
    __biocyc_login: bool = False
    # Message from BioCyc if it rejected the login information
    __biocyc_login_error: Optional[str] = None
    add_smiles_as_cross_reference: bool = False
    """
    SMILES strings are not registered in indentifiers.org therefore they are
    not considered as valid identifiers and will not be added by default to metabolites.
    If you want to add SMILES strings as cross-references, set this to True.

    .. versionadded:: 1.4.0

    .. note::
        This setting only affects cross-references that are retrieved from
        MetaNetX after it is changed. Cross-references that are already
        cached in the data directory are used as they are, with or without
        SMILES.
    """

    add_metanetx_specific_prefixes: bool = False
    """
    MetaNetX returns most identifiers two times: once with the prefix
    registered in identifiers.org and once with a MetaNetX specific prefix,
    e.g. 'kegg.compound:C00031' and 'keggC:C00031'. Cross-references with a
    MetaNetX specific prefix, such as 'keggC' or 'biggM', will not be added
    by default. If you want to add them as cross-references, set this to True.

    .. versionadded:: 1.4.0
    """

    auto_open_close_biocyc_session = True

    def __set__biocyc_password(self, value: str):
        if value != self.__biocyc_password:
            self.__biocyc_password = value
            # The next request logs in again with the new credentials
            self._close_session("BIOCYC")

    biocyc_password = property(fset=__set__biocyc_password)
    del __set__biocyc_password

    def __set__biocyc_name(self, value: str):
        if value != self.__biocyc_name:
            self.__biocyc_name = value
            # The next request logs in again with the new credentials
            self._close_session("BIOCYC")

    biocyc_name = property(fset=__set__biocyc_name)
    del __set__biocyc_name

    def set_biocyc_login(self, file: Path):
        user, pwd = cmod_utils.get_credentials(file)
        self.biocyc_password = pwd
        self.biocyc_name = user

    def __set_biocyc_login_from_file(self):
        try:
            self.set_biocyc_login(Path.cwd().joinpath("credentials.txt"))
        except FileNotFoundError:
            logger.debug("No credentials.txt found.")

    def __set_biocyc_login_from_env(self):
        biocyc_name = os.getenv("BIOCYC_NAME")
        biocyc_pwd = os.getenv("BIOCYC_PASSWORD")

        if biocyc_name and biocyc_pwd:
            self.__biocyc_name = biocyc_name
            self.__biocyc_password = biocyc_pwd
            logger.debug("BioCyc credentials set from environment variables.")

    def _session(self, database: str) -> requests.Session:
        """
        Returns the session for a database and creates it if necessary. Each
        database has its own session, so that login information is only sent
        to the database it belongs to.

        Args:
            database (str): The database as named in
                cobramod.retrieval.databases, i.e. "BIOCYC", "PMN", "SOL",
                "KEGG" or "BIGG".

        Returns:
            requests.Session: The session. For "BIOCYC", it is logged in if
                login information is available.
        """
        session = self.__sessions.get(database)

        if session is None:
            logger.debug(f"Creating new {database} session ...")
            session = requests.Session()
            retries = HTTPAdapter(
                max_retries=Retry(
                    total=5,
                    backoff_factor=0.5,
                    allowed_methods=frozenset(["GET", "POST"]),
                )
            )
            session.mount("https://", retries)
            session.mount("http://", retries)
            self.__sessions[database] = session

        if database == "BIOCYC":
            self.__login_biocyc(session)

        return session

    def __login_biocyc(self, session: requests.Session):
        if (
            not self.__biocyc_login
            and self.__biocyc_login_error is None
            and self._has_biocyc_credentials
        ):
            logger.debug("Using BioCyc credentials to login ...")
            response = session.post(
                "https://websvc.biocyc.org/credentials/login/",
                data={
                    "email": self.__biocyc_name,
                    "password": self.__biocyc_password,
                },
            )

            if response.ok:
                self.__biocyc_login = True
            elif response.status_code == 401:
                # Wrong login information. It is not sent again until new
                # login information is set
                try:
                    self.__biocyc_login_error = response.json()["error"]
                except (ValueError, KeyError, TypeError):
                    self.__biocyc_login_error = response.text[:200]
                logger.warning(
                    f"BioCyc login failed: {self.__biocyc_login_error}"
                )
            else:
                logger.warning(
                    f"BioCyc login failed with status {response.status_code}."
                    " It is tried again with the next request."
                )

    @property
    def _has_biocyc_credentials(self) -> bool:
        return bool(self.__biocyc_name and self.__biocyc_password)

    @property
    def _biocyc_login_error(self) -> Optional[str]:
        """
        Message from BioCyc if it rejected the login information, otherwise
        None.
        """
        return self.__biocyc_login_error

    @property
    def biocyc_logged_in(self) -> bool:
        """
        Check if the BioCyc session is logged in.

        Returns:
            bool: True if logged in, False otherwise.
        """
        return self.__biocyc_login

    @property
    def cache_dir(self) -> Path:
        return self.__cache_dir

    def _close_session(self, database: str):
        """
        Closes the session of a database. The next request creates a new one.
        """
        session = self.__sessions.pop(database, None)
        if session is not None:
            session.close()

        if database == "BIOCYC":
            # A new session tries to log in again
            self.__biocyc_login = False
            self.__biocyc_login_error = None

    def __del__(self):
        for database in list(self.__sessions):
            self._close_session(database)
