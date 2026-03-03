import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas
from cobra import Metabolite, Model, Reaction
from rdkit.VLib.NodeLib.SmartsRemover import biggerTest
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

import cobramod
from cobramod.core.crossreferences import add2dict_unique, validate_id
from cobramod.dbwalker.bigg import Bigg
from cobramod.dbwalker.BioCyc import BioCyc
from cobramod.dbwalker.chebi import Chebi
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable
from cobramod.dbwalker.kegg import Kegg
from cobramod.dbwalker.pubchem import PubChem

general_identifiers = ["inchikey", "inchi", "smiles"]

logger = logging.getLogger("cobramod.DBWalker.CrossReferences")
logger.propagate = True

settings = cobramod.Settings()

supportedDatabases: List[Database] = [
    BioCyc(),
    Kegg(),
    # Bigg(),
    # Chebi(),
    PubChem(),
]


class AnnotationLogger:
    dict = {}

    def add_entry(
        self, metabolite_id, db_name, GenerellIdentifiers: GenerellIdentifiers
    ):
        self.dict[(metabolite_id, db_name)] = GenerellIdentifiers.to_dict()

    def to_dataframe(self):
        return pandas.DataFrame(self.dict)


def add_crossreferences2metabolite(
    metabolite: Metabolite, annotation_logger: Optional[AnnotationLogger] = None
):
    general_identifiers = GenerellIdentifiers.fromAnnotation(
        metabolite.annotation
    )

    logger.debug(
        f"For Metabolite ({metabolite.id}) the following GID was already present {general_identifiers}"
    )

    for database in supportedDatabases:
        databaseID = metabolite.annotation.get(database.AnnotationPrefix, None)
        if databaseID is None:
            continue

        general_identifier = database.getGenerellIdentifier(databaseID)

        if general_identifier is Unavailable:
            logger.debug(
                f"No GIDs available in {database.name} for ID {databaseID}"
            )
            continue

        logger.debug(
            f"Found in {database.name} the following GIDs {general_identifier}."
        )

        if annotation_logger:
            annotation_logger.add_entry(
                metabolite_id=metabolite.id,
                db_name=database.name,
                GenerellIdentifiers=general_identifiers,
            )

        if general_identifiers.empty():
            general_identifiers = general_identifier
        elif general_identifiers.weakEq(general_identifier):
            general_identifiers += general_identifier
        else:
            logger.error(
                f"While adding general identifiers from {database.name} for ID {databaseID}, a missmatch occurred."
                f" Previous identifiers are: {general_identifiers}, Database IDs are {general_identifier}."
            )
            return

    for database in supportedDatabases:
        logger.debug(
            f"Getting DB ID for database {database.name} and object {metabolite.id}"
        )

        dbID = database.getDBIdentifier(general_identifiers)

        if dbID is Unavailable:
            continue

        assert dbID != "Unavailable"
        assert dbID is not None
        assert dbID is not Unavailable

        add2dict_unique(
            key=database.AnnotationPrefix,
            value=dbID,
            dictionary=metabolite.annotation,
        )

    if "inchi" not in metabolite.annotation and isinstance(
        general_identifiers.inchi, str
    ):
        metabolite.annotation["inchi"] = general_identifiers.inchi
    if "inchikey" not in metabolite.annotation and isinstance(
        general_identifiers.inchi_key, str
    ):
        metabolite.annotation["inchikey"] = general_identifiers.inchi_key
    if "smiles" not in metabolite.annotation and isinstance(
        general_identifiers.smiles, str
    ):
        metabolite.annotation["smiles"] = general_identifiers.smiles

    for database in supportedDatabases:
        database.save_cache()


def id2annotation(object: Model):
    for metabolite in object.metabolites:
        ID = metabolite.id
        metabolite.annotation["kegg.compound"] = ID


def add_crossreferences(
    object: Union[Model, Reaction, Metabolite],
    consider_subobjects: bool = True,
    save_log: Union[None, str, Path] = None,
):
    """
    Add cross-references to the given object. This function works with a cobrapy
    Model, Reaction, or Metabolite.

    Args:
        object (Union[Model, Reaction, Metabolite]): The COBRApy object to which
            cross-references should be added.
    """
    autoOpenCloseBioCycSession = settings.autoOpenCloseBioCycSession
    settings.autoOpenCloseBioCycSession = False
    annotation_logger = AnnotationLogger()

    if isinstance(object, Model):
        with logging_redirect_tqdm():
            if not consider_subobjects:
                return

            for reaction in object.reactions:
                pass

            for metabolite in object.metabolites:  # tqdm(object.metabolites):
                add_crossreferences2metabolite(
                    metabolite=metabolite, annotation_logger=annotation_logger
                )

    elif isinstance(object, Reaction):
        pass

    elif isinstance(object, Metabolite):
        add_crossreferences2metabolite(metabolite=object)
    else:
        raise TypeError("Object must be a Model, Reaction, or Metabolite.")

    # close the session and restore original behaviour
    settings._closeBiocycSession()
    settings.autoOpenCloseBioCycSession = autoOpenCloseBioCycSession

    if save_log:
        df = annotation_logger.to_dataframe()
        df.to_csv(
            path_or_buf=save_log,
        )
