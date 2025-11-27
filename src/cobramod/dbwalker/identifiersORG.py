import logging
import re
from typing import Dict, Optional

import requests

from cobramod.settings import SingletonMeta

logger = logging.getLogger("cobramod.DBWalker.IdentifiersORG")
logger.propagate = True

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


class Validator(metaclass=SingletonMeta):
    def __init__(self):
        super().__init__()

        self.lookup: Dict[str, re.Pattern] = self.get_registry()

    def get_registry(self) -> Dict[str, re.Pattern]:
        """
        Getting the registry information from identifiers.org and constructing the lookup for the regex patterns of
        the IDs as defined by identifiers.org.
        Returns:

        """

        url = "https://registry.api.identifiers.org/resolutionApi/getResolverDataset?rewriteForEmbeddedPrefixes=true"

        logger.debug(f"Getting registry info from {url}")

        response = requests.get(url)
        response.raise_for_status()

        data = response.json()["payload"]["namespaces"]
        lookup: Dict[str, re.Pattern] = {}

        for entry in data:
            prefix = entry["prefix"]
            pattern = re.compile(entry["pattern"])

            lookup[prefix] = pattern

        logger.info(
            "Sucessfully retrieved registry info and created regex pattern."
        )
        return lookup

    def validate_id_pattern(self, prefix: str, identifier: str) -> bool:
        logger.debug(f"Validating {prefix}:{identifier}")

        try:
            pattern = self.lookup[prefix]
        except KeyError:
            logger.debug(f"No entry found for {prefix}")
            return False

        found = False
        if pattern.match(f"{identifier}") is not None:
            found = True
            logger.info(f"Structure of {prefix}:{identifier} is valid")
        else:
            logger.info(f"Structure of {prefix}:{identifier} is invalid")

        return found
