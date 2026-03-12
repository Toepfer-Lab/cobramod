import logging
import tarfile
from pathlib import Path
from typing import Optional, Union

import pandas as pd
import requests

from cobramod import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import (
    GenerellIdentifiers,
    Unavailable,
    UnavailableType,
    Uncertain,
)

logger = logging.getLogger("cobramod.DBWalker.MetaNetX")
logger.propagate = True


class MetaNetX(Database):
    """
    Database walker for MetaNetX. Uses the TSV flat files from MetaNetX's FTP
    server (chem_prop.tsv for chemical properties, chem_xref.tsv for
    cross-references) to resolve identifiers locally.
    """

    MNX_FTP = "https://www.metanetx.org/ftp/latest"

    def __init__(self):
        super().__init__()
        self.__settings = Settings()
        self.__cache_dir = self.__settings.cacheDir / self.name
        self.__cache_dir.mkdir(parents=True, exist_ok=True)

        self.__download_files()
        self.__load_files()

    @property
    def name(self) -> str:
        return "MetaNetX"

    @property
    def AnnotationPrefix(self) -> str:
        return "metanetx.chemical"

    def save_cache(self):
        pass

    # ── File management ──────────────────────────────────────────────

    def __download_files(self):
        """Download mnxref_tsv.tar.gz and extract the needed TSV files."""
        # Check if the needed files are already extracted
        chem_prop = self.__cache_dir / "chem_prop.tsv"
        chem_xref = self.__cache_dir / "chem_xref.tsv"
        if chem_prop.exists() and chem_xref.exists():
            logger.debug("Found local MetaNetX TSV files")
            return

        tar_path = self.__cache_dir / "mnxref_tsv.tar.gz"

        if not tar_path.exists():
            url = f"{self.MNX_FTP}/mnxref_tsv.tar.gz"
            logger.info("Downloading mnxref_tsv.tar.gz from MetaNetX (~200MB)...")

            self.__settings.limiter.try_acquire("metanetx")
            response = requests.get(url, timeout=300, stream=True)
            response.raise_for_status()

            with open(tar_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            logger.info(f"Saved mnxref_tsv.tar.gz at {tar_path}")

        logger.info("Extracting TSV files from mnxref_tsv.tar.gz...")
        with tarfile.open(tar_path, "r:gz") as tar:
            for member in tar.getmembers():
                basename = Path(member.name).name
                if basename in ("chem_prop.tsv", "chem_xref.tsv"):
                    member.name = basename
                    tar.extract(member, path=self.__cache_dir)

        logger.info("Extracted files for MetaNetX")

    @staticmethod
    def __parse_header(filepath: Path) -> tuple[int, list[str]]:
        """
        Parse leading comment lines of a MetaNetX TSV file.

        Returns (skip_rows, column_names) where skip_rows is the number of
        comment lines and column_names are parsed from the last comment line
        (which serves as the header).
        """
        count = 0
        last_comment = ""
        with open(filepath) as f:
            for line in f:
                if line.startswith("#"):
                    last_comment = line
                    count += 1
                else:
                    break
        columns = last_comment.lstrip("#").strip().split("\t")
        return count, columns

    def __load_files(self):
        # Cannot use comment="#" because '#' appears inside SMILES (triple bonds)
        prop_path = self.__cache_dir / "chem_prop.tsv"
        skip, cols = self.__parse_header(prop_path)
        self._chem_prop = pd.read_csv(
            prop_path,
            sep="\t",
            skiprows=skip,
            header=None,
            names=cols,
            index_col=False,
        )

        xref_path = self.__cache_dir / "chem_xref.tsv"
        skip, cols = self.__parse_header(xref_path)
        self._chem_xref = pd.read_csv(
            xref_path,
            sep="\t",
            skiprows=skip,
            header=None,
            names=cols,
            index_col=False,
        )

    # ── Core lookups ─────────────────────────────────────────────────

    def getGenerellIdentifier(
        self, dbIdentifier: str, **kwargs
    ) -> GenerellIdentifiers:
        """
        Get InChI, InChIKey, and SMILES for a MetaNetX compound ID (e.g. MNXM3).
        """
        row = self._chem_prop.loc[self._chem_prop["ID"] == dbIdentifier]

        if len(row) == 0:
            logger.warning(f"MetaNetX ID '{dbIdentifier}' does not exist in chem_prop.tsv")
            return Unavailable

        row = row.iloc[0]

        return GenerellIdentifiers(
            inchi=row["InChI"] if pd.notna(row["InChI"]) else None,
            inchi_key=row["InChIKey"] if pd.notna(row["InChIKey"]) else None,
            smiles=row["SMILES"] if pd.notna(row["SMILES"]) else None,
        )

    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        if isinstance(smiles, GenerellIdentifiers):
            smiles = smiles.smiles

        result = self._chem_prop.loc[
            self._chem_prop["SMILES"] == smiles, "ID"
        ]
        return self._single_or_uncertain_or_unavailable(result, "SMILES", smiles)

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        if isinstance(inchi, GenerellIdentifiers):
            inchi = inchi.inchi

        result = self._chem_prop.loc[
            self._chem_prop["InChI"] == inchi, "ID"
        ]
        return self._single_or_uncertain_or_unavailable(result, "InChI", inchi)

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> Union[str, UnavailableType]:
        if isinstance(inchikey, GenerellIdentifiers):
            inchikey = inchikey.inchi_key

        result = self._chem_prop.loc[
            self._chem_prop["InChIKey"] == inchikey, "ID"
        ]
        return self._single_or_uncertain_or_unavailable(result, "InChIKey", inchikey)

    # ── Cross-reference lookups ──────────────────────────────────────

    def get_xrefs(self, mnx_id: str) -> list[str]:
        """
        Get all cross-references for a MetaNetX compound ID.

        Returns a list of strings like ["kegg.compound:C00002", "chebi:15422", ...].
        """
        rows = self._chem_xref.loc[self._chem_xref["ID"] == mnx_id, "source"]
        return rows.tolist()

    def get_mnx_id_from_external(self, external_id: str) -> Optional[str]:
        """
        Resolve an external ID (e.g. "kegg.compound:C00002") to a MetaNetX ID.

        Returns the MNX ID, Uncertain if multiple, or Unavailable if none.
        """
        result = self._chem_xref.loc[
            self._chem_xref["source"] == external_id, "ID"
        ]

        unique = result.unique()

        if len(unique) == 1:
            return unique[0]
        elif len(unique) > 1:
            logger.warning(
                f"Multiple MNX IDs found for external ID ({external_id}): {unique.tolist()}"
            )
            return Uncertain(type="DBID", possibilities=unique.tolist())
        else:
            return Unavailable

    # ── Helpers ───────────────────────────────────────────────────────

    def _single_or_uncertain_or_unavailable(self, result, field: str, query: str):
        if len(result) == 1:
            return str(result.iloc[0])
        elif len(result) > 1:
            logger.warning(
                f"Multiple ({len(result)}) entries for {field} ({query}) in MetaNetX"
            )
            return Uncertain(
                type="DBID",
                possibilities=result.astype(str).tolist(),
            )
        else:
            logger.debug(f"No match for {field} ({query}) in MetaNetX")
            return Unavailable
