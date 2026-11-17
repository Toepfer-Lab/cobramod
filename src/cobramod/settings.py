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
    __biocyc_session: Optional[requests.Session] = None
    __biocyc_login: bool = False
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
            self._close_biocyc_session()

    biocyc_password = property(fset=__set__biocyc_password)
    del __set__biocyc_password

    def __set__biocyc_name(self, value: str):
        if value != self.__biocyc_name:
            self.__biocyc_name = value
            # The next request logs in again with the new credentials
            self._close_biocyc_session()

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

    @property
    def _biocyc_session(self) -> requests.Session:
        if self.__biocyc_session is None:
            logger.debug("Creating new BioCyc Session ...")
            self.__biocyc_session = requests.Session()
            self.__biocyc_session.mount(
                "https://",
                HTTPAdapter(
                    max_retries=Retry(
                        total=5,
                        backoff_factor=0.5,
                        allowed_methods=frozenset(["GET", "POST"]),
                    )
                ),
            )
        if (
            not self.__biocyc_login
            and self.__biocyc_name
            and self.__biocyc_password
        ):
            logger.debug("Using BioCyc credentials to login ...")
            self.__biocyc_session.post(
                "https://websvc.biocyc.org/credentials/login/",
                data={
                    "email": self.__biocyc_name,
                    "password": self.__biocyc_password,
                },
            )

            self.__biocyc_login = True

        return self.__biocyc_session

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

    def _close_biocyc_session(self):
        if self.__biocyc_session:
            self.__biocyc_session.close()
            self.__biocyc_session = None
            self.__biocyc_login = False

    def __del__(self):
        self._close_biocyc_session()
