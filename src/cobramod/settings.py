import logging
import os
import sys
from pathlib import Path
from threading import Lock
from typing import Optional

import platformdirs
import requests
from pyrate_limiter import Limiter, Rate, Duration
from requests.adapters import HTTPAdapter

import cobramod.utils as cmod_utils

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)]
)

logger = logging.getLogger("cobramod.Settings")
logger.setLevel(logging.DEBUG)
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
        self.__setBioCycLoginFromEnv()

        if self.__biocyc_name is None or self.__biocyc_password is None:
            self.__setBioCycLoginFromFile()

        self.__cache_dir: Path = Path(
            platformdirs.user_cache_dir(
                appname="cobramod",
                appauthor="cobramod",
                ensure_exists=True,
            )
        )

        self.limiter = Limiter(Rate(1, Duration.SECOND),raise_when_fail=False, max_delay=3600)

    __biocyc_password: Optional[str] = None
    __biocyc_name: Optional[str] = None
    __biocyc_session: Optional[requests.Session] = None
    __biocyc_login: bool = False
    add_smiles_as_cross_reference: bool = False
    """
    SMILES strings are not registered in indentifiers.org therefore they are
    not considered as valid identifiers and will not be added by default to metabolites.
    If you want to add SMILES strings as cross-references, set this to True.
    """

    autoOpenCloseBioCycSession = True

    def __set__biocyc_password(self, value: str):
        self.__biocyc_password = value

    biocyc_password = property(fset=__set__biocyc_password)
    del __set__biocyc_password

    def __set__biocyc_name(self, value: str):
        self.__biocyc_name = value

    biocyc_name = property(fset=__set__biocyc_name)
    del __set__biocyc_name

    def SetBioCycLogin(self, file: Path):
        user, pwd = cmod_utils.get_credentials(file)
        self.__biocyc_password = pwd
        self.__biocyc_name = user

    def __setBioCycLoginFromFile(self):
        try:
            self.SetBioCycLogin(Path.cwd().joinpath("credentials.txt"))
        except FileNotFoundError:
            logger.debug("No credentials.txt found.")

    def __setBioCycLoginFromEnv(self):
        biocyc_name = os.getenv("BIOCYC_NAME")
        biocyc_pwd = os.getenv("BIOCYC_PASSWORD")

        if biocyc_name and biocyc_pwd:
            self.__biocyc_name = biocyc_name
            self.__biocyc_password = biocyc_pwd
            logger.debug("BioCyc credentials set from environment variables.")

    @property
    def _biocycSession(self) -> requests.Session:
        if self.__biocyc_session is None:
            logger.debug("Creating new BioCyc Session ...")
            self.__biocyc_session = requests.Session()
            self.__biocyc_session.mount('https://', HTTPAdapter(max_retries=5))

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
    def BioCycLoggedIn(self) -> bool:
        """
        Check if the BioCyc session is logged in.

        Returns:
            bool: True if logged in, False otherwise.
        """
        return self.__biocyc_login

    @property
    def cacheDir(self) -> Path:
        return self.__cache_dir

    def _closeBiocycSession(self):
        if self.__biocyc_session:
            self.__biocyc_session.close()
            self.__biocyc_session = None
            self.__biocyc_login = False

    def __del__(self):
        self._closeBiocycSession()
