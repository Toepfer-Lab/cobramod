import logging
from pathlib import Path
from typing import Union, Dict, List, Any, Optional

import pandas
from cobra import Model, Reaction, Metabolite
from rdkit.VLib.NodeLib.SmartsRemover import biggerTest
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

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

settings = cobramod.Settings()

supportedDatabases: List[Database]=[
    BioCyc(),
    Kegg(),
    #Bigg(),
    #Chebi(),
    PubChem(),
]

class AnnotationLogger():
    dict = {}

    def add_entry(self, metabolite_id,db_name,GenerellIdentifiers: GenerellIdentifiers):
        self.dict[(metabolite_id, db_name)] = GenerellIdentifiers


    def to_dataframe(self):

        return pandas.DataFrame(self.dict)




def add_crossreferences2metabolite(metabolite: Metabolite, annotation_logger: Optional[AnnotationLogger] = None):
    general_identifiers = GenerellIdentifiers.fromAnnotation(metabolite.annotation)

    for database in supportedDatabases:
        databaseID = metabolite.annotation.get(database.AnnotationPrefix, None)
        if databaseID is None:
            continue

        general_identifier = database.getGenerellIdentifier(databaseID)
        if annotation_logger:
            annotation_logger.add_entry(
                metabolite_id = metabolite.id,
                db_name= database.name,
                GenerellIdentifiers=general_identifiers,
            )

        if general_identifiers.empty():
            general_identifiers = general_identifier
        elif general_identifiers.weakEq(general_identifier):
            general_identifiers += general_identifier
        else:
            logger.error(
                f"While adding general identifiers from {database.name}, a missmatch occurred."
                f"Previous identifiers are: {general_identifiers}, Database IDs are {general_identifier}."
            )
            return


    for database in supportedDatabases:

        logger.debug(
            f"Getting DB ID for database {database.name} and object {metabolite.name}"
        )

        dbID = database.getDBIdentifier(general_identifiers)

        if dbID is None:
            continue

        add2dict_unique(
            key=database.AnnotationPrefix,
            value=dbID,
            dictionary=metabolite.annotation,
        )

def id2annotation(object: Model):

    for metabolite in object.metabolites:
        ID = metabolite.id
        metabolite.annotation["kegg.compound"] = ID

def add_crossreferences(
        object: Union[Model, Reaction, Metabolite],
        consider_subobjects: bool = True,
        save_log:Optional[Path] = None
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

            for metabolite in object.metabolites: #tqdm(object.metabolites):
                add_crossreferences2metabolite(metabolite=metabolite, annotation_logger= annotation_logger)

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
            path_or_buf= save_log,
        )