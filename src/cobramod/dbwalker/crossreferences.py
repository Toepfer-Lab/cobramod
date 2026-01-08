import logging
from typing import Union, Dict, List, Any

from cobra import Model, Reaction, Metabolite
from rdkit.VLib.NodeLib.SmartsRemover import biggerTest

import cobramod
from cobramod.core.crossreferences import validate_id, add2dict_unique
from cobramod.dbwalker.BioCyc import BioCyc
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.bigg import Bigg
from cobramod.dbwalker.chebi import Chebi
from cobramod.dbwalker.dataclasses import GenerellIdentifiers
from cobramod.dbwalker.kegg import Kegg

from cobramod.dbwalker.pubchem import PubChem

general_identifiers = ["inchikey", "inchi", "smiles"]

logger = logging.getLogger("cobramod.DBWalker.CrossReferences")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

settings = cobramod.Settings()

supportedDatabases: List[Database]=[
    BioCyc(),
    Kegg(),
    #Bigg(),
    #Chebi(),
    PubChem(),
]

def add_crossreferences2metabolite(metabolite: Metabolite):
    annotations = metabolite.annotation
    general_identifiers = GenerellIdentifiers.fromAnnotation(metabolite.annotation)

    for database in supportedDatabases:
        databaseID = metabolite.annotation.get(database.AnnotationPrefix, None)
        if databaseID is None:
            continue

        general_identifier = database.getGenerellIdentifier(databaseID)

        if general_identifiers.empty():
            general_identifiers = general_identifier
        elif general_identifiers.weakEq(general_identifier):
            general_identifiers += general_identifier
        else:
            logger.error(
                f"While adding general identifiers from {database.name}, a missmatch occurred."
                f"Previous identifiers are: {general_identifiers}, Database IDs are {general_identifier}."
            )


    for database in supportedDatabases:
        dbID = database.getDBIdentifier(general_identifiers)
        add2dict_unique(
            key=database.AnnotationPrefix,
            value=dbID,
            dictionary=metabolite.annotation,
        )

def add_crossreferences(
    object: Union[Model, Reaction, Metabolite], consider_subobjects: bool = True
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

    if isinstance(object, Model):
        if not consider_subobjects:
            return

        for reaction in object.reactions:
            pass

        for metabolite in object.metabolites:
            add_crossreferences2metabolite(metabolite=metabolite)

    elif isinstance(object, Reaction):
        pass

    elif isinstance(object, Metabolite):
        add_crossreferences2metabolite(metabolite=object)
    else:
        raise TypeError("Object must be a Model, Reaction, or Metabolite.")

    # close the session and restore original behaviour
    settings._closeBiocycSession()
    settings.autoOpenCloseBioCycSession = autoOpenCloseBioCycSession
