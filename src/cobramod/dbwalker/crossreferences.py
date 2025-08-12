import logging
from typing import Union

from cobra import Model, Reaction, Metabolite

import cobramod
from cobramod.core.crossreferences import validate_id, add2dict_unique
from cobramod.dbwalker.BioCyc import get_compound_info_by_biocyc_id, smiles2BioCyc, InChI2BioCyc
from cobramod.dbwalker.chebi import get_chebi_compound_info, getChebiID
from cobramod.dbwalker.modelSeed import get_compound_info_by_modelseed_id, smiles2ModelSeed, inchikey2modelseed

general_identifiers = ["inchikey", "inchi", "smiles"]

logger = logging.getLogger("cobramod.DBWalker.CrossReferences")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

settings = cobramod.Settings()

def add_crossreferences2metabolite(metabolite:Metabolite):

    annotations = metabolite.annotation
    present_identifiers = {}
    verified_identifiers = {}
    nonverified_identifiers = {}

    for identifier in general_identifiers:
        if identifier in annotations:
            value = annotations[identifier]
            if isinstance(value, str):
                value = [value]
            present_identifiers[identifier] = value


    #biocyc
    if "biocyc" in annotations:
        IDs = annotations["biocyc"]


        if isinstance(IDs, str):
            IDs = [IDs]

        for ID in IDs:
            db, identifier = ID.split(":")
            biocyc_result = get_compound_info_by_biocyc_id(identifier, db=db)

            if biocyc_result.smiles is not None:
                query = smiles2BioCyc(biocyc_result.smiles, orgid=db)
                identifiersORG_validation = validate_id(f"biocyc:{db}:{query}")

                if present_identifiers.get("smiles") is not None:
                    if biocyc_result.smiles != present_identifiers["smiles"]:
                        logger.error(f"SMILES identifier mismatch for {identifier} in BioCyc: {biocyc_result.smiles} vs existing {present_identifiers['smiles']}")
                        logger.error(f"Please check the BioCyc ID ({ID}) and the present SMILES ({present_identifiers['smiles']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the SMILES match, we can ignore it as it exists already
                        continue

                if query == identifier and identifiersORG_validation:
                    logger.info(
                        f"Found SMILES {biocyc_result.smiles} in BioCyc ({db}) for {metabolite.id}"
                    )
                    verified_identifiers = add2dict_unique(
                        key= "smiles",
                        value= biocyc_result.smiles,
                        dictionary=verified_identifiers
                    )

            if biocyc_result.inchi is not None:
                query = InChI2BioCyc(biocyc_result.inchi, orgid=db)
                identifiersORG_validation = validate_id(f"biocyc:{db}:{query}")

                if present_identifiers.get("inchi") is not None:
                    if biocyc_result.inchi != present_identifiers["inchi"]:
                        logger.error(f"InChI identifier mismatch for {identifier} in BioCyc: {biocyc_result.inchi} vs existing {present_identifiers['inchi']}")
                        logger.error(f"Please check the BioCyc ID ({ID}) and the present inchi ({present_identifiers['inchi']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the InChI match, we can ignore it as it exists already
                        continue

                if query == identifier and identifiersORG_validation:
                    logger.info(
                        f"Found InChI {biocyc_result.inchi} in BioCyc ({db}) for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="inchi",
                        value=biocyc_result.inchi,
                        dictionary=verified_identifiers
                    )
            if biocyc_result.inchi_key is not None:
                if present_identifiers.get("inchikey") is not None:
                    if present_identifiers.get("inchikey") == biocyc_result.inchi_key:
                        # If the InChIKey match, we can ignore it as it exists already
                        continue
                    else:
                        logger.error(f"InChIKey identifier mismatch for {identifier} in BioCyc: {biocyc_result.inchi_key} vs existing {present_identifiers['inchikey']}")
                        logger.error(f"Please check the BioCyc ID ({ID}) and the present InChIKey ({present_identifiers['inchikey']}) for {metabolite.id}.")
                        continue

                logger.info(
                    f"Found unverified InChIKey {biocyc_result.inchi_key} in BioCyc ({db}) for {metabolite.id}"
                )
                nonverified_identifiers["inchikey"] = biocyc_result.inchi_key

    if "metacyc.compound" in annotations:
        db = "META"
        IDs = annotations["metacyc.compound"]

        if isinstance(IDs, str):
            IDs = [IDs]

        for ID in IDs:
            biocyc_result = get_compound_info_by_biocyc_id(ID, db=db)

            if biocyc_result.smiles is not None:
                query = smiles2BioCyc(biocyc_result.smiles, orgid=db)
                identifiersORG_validation = validate_id(f"biocyc:{db}:{query}")

                if present_identifiers.get("smiles") is not None:
                    if biocyc_result.smiles != present_identifiers["smiles"]:
                        logger.error(f"SMILES identifier mismatch for {ID} in BioCyc: {biocyc_result.smiles} vs existing {present_identifiers['smiles']}")
                        logger.error(f"Please check the BioCyc ID ({db}:{ID}) and the present SMILES ({present_identifiers['smiles']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the SMILES match, we can ignore it as it exists already
                        continue

                if query == ID and identifiersORG_validation:
                    logger.info(
                        f"Found SMILES {biocyc_result.smiles} in BioCyc ({db}) for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="smiles",
                        value=biocyc_result.smiles,
                        dictionary=verified_identifiers
                    )

            if biocyc_result.inchi is not None:
                query = InChI2BioCyc(biocyc_result.inchi, orgid=db)
                identifiersORG_validation = validate_id(f"biocyc:{db}:{query}")

                if present_identifiers.get("inchi") is not None:
                    if biocyc_result.inchi != present_identifiers["inchi"]:
                        logger.error(f"InChI identifier mismatch for {ID} in BioCyc: {biocyc_result.inchi} vs existing {present_identifiers['inchi']}")
                        logger.error(f"Please check the BioCyc ID ({db}:{ID}) and the present inchi ({present_identifiers['inchi']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the InChI match, we can ignore it as it exists already
                        continue

                if query == ID and identifiersORG_validation:
                    logger.info(
                        f"Found InChI {biocyc_result.inchi} in BioCyc ({db}) for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="inchi",
                        value=biocyc_result.inchi,
                        dictionary=verified_identifiers
                    )

            if biocyc_result.inchi_key is not None:
                if present_identifiers.get("inchikey") is not None:
                    if present_identifiers.get("inchikey") == biocyc_result.inchi_key:
                        # If the InChIKey match, we can ignore it as it exists already
                        continue
                    else:
                        logger.error(f"InChIKey identifier mismatch for {ID} in BioCyc: {biocyc_result.inchi_key} vs existing {present_identifiers['inchikey']}")
                        logger.error(f"Please check the BioCyc ID ({ID}) and the present InChIKey ({present_identifiers['inchikey']}) for {metabolite.id}.")
                        continue

                logger.info(
                    f"Found unverified InChIKey {biocyc_result.inchi_key} in BioCyc ({db}) for {metabolite.id}"
                )
                nonverified_identifiers["inchikey"] = biocyc_result.inchi_key

    if "chebi" in annotations:
        IDs = annotations["chebi"]

        if isinstance(IDs, str):
            IDs = [IDs]

        for ID in IDs:
            chebi_result = get_chebi_compound_info(ID)

            if chebi_result.smiles is not None:
                if present_identifiers.get("smiles") is not None:
                    if chebi_result.smiles != present_identifiers["smiles"]:
                        logger.error(f"SMILES identifier mismatch for {ID} in ChEBI: {chebi_result.smiles} vs existing {present_identifiers['smiles']}")
                        logger.error(f"Please check the ChEBI ID ({ID}) and the present SMILES ({present_identifiers['smiles']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the SMILES match, we can ignore it as it exists already
                        continue

                query = getChebiID(chebi_result.smiles, identifier_type="SMILES")
                identifiersORG_validation = validate_id(f"chebi:{query}")

                if query == ID and identifiersORG_validation:
                    logger.info(
                        f"Found SMILES {chebi_result.smiles} in ChEBI for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="smiles",
                        value=chebi_result.smiles,
                        dictionary=verified_identifiers
                    )
            if chebi_result.inchi is not None:
                if present_identifiers.get("inchi") is not None:
                    if chebi_result.inchi != present_identifiers["inchi"]:
                        logger.error(f"InChI identifier mismatch for {ID} in ChEBI: {chebi_result.inchi} vs existing {present_identifiers['inchi']}")
                        logger.error(f"Please check the ChEBI ID ({ID}) and the present inchi ({present_identifiers['inchi']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the InChI match, we can ignore it as it exists already
                        continue

                query = getChebiID(chebi_result.inchi, identifier_type="InChI")
                identifiersORG_validation = validate_id(f"chebi:{query}")

                if query == ID and identifiersORG_validation:
                    logger.info(
                        f"Found InChI {chebi_result.inchi} in ChEBI for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="inchi",
                        value=chebi_result.inchi,
                        dictionary=verified_identifiers
                    )

            if chebi_result.inchi_key is not None:
                if present_identifiers.get("inchikey") is not None:
                    if present_identifiers.get("inchikey") == chebi_result.inchi_key:
                        # If the InChIKey match, we can ignore it as it exists already
                        continue
                    else:
                        logger.error(f"InChIKey identifier mismatch for {ID} in ChEBI: {chebi_result.inchi_key} vs existing {present_identifiers['inchikey']}")
                        logger.error(f"Please check the ChEBI ID ({ID}) and the present InChIKey ({present_identifiers['inchikey']}) for {metabolite.id}.")
                        continue

                query = getChebiID(chebi_result.inchi_key, identifier_type="InChIKey")
                identifiersORG_validation = validate_id(f"chebi:{query}")

                if query == ID and identifiersORG_validation:
                    logger.info(
                        f"Found InChIKey {chebi_result.inchi_key} in ChEBI for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="inchikey",
                        value=chebi_result.inchi,
                        dictionary=verified_identifiers
                    )
    if "seed.compound" in annotations:
        IDs = annotations["modelseed.compound"]

        if isinstance(IDs, str):
            IDs = [IDs]

        for ID in IDs:
            modelseed_result = get_compound_info_by_modelseed_id(ID)

            if modelseed_result.smiles is not None:
                if present_identifiers.get("smiles") is not None:
                    if modelseed_result.smiles != present_identifiers["smiles"]:
                        logger.error(f"SMILES identifier mismatch for {ID} in ModelSEED: {modelseed_result.smiles} vs existing {present_identifiers['smiles']}")
                        logger.error(f"Please check the ModelSEED ID ({ID}) and the present SMILES ({present_identifiers['smiles']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the SMILES match, we can ignore it as it exists already
                        continue

                query = smiles2ModelSeed(modelseed_result.smiles)
                identifiersORG_validation = validate_id(f"seed.compound:{query}")

                if query == ID and identifiersORG_validation:
                    logger.info(
                        f"Found SMILES {modelseed_result.smiles} in ModelSEED for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="smiles",
                        value=modelseed_result.smiles,
                        dictionary=verified_identifiers
                    )
            if modelseed_result.inchi_key is not None:
                if present_identifiers.get("inchikey") is not None:
                    if modelseed_result.inchi_key != present_identifiers["inchikey"]:
                        logger.error(f"InChIKey identifier mismatch for {ID} in ModelSEED: {modelseed_result.inchi} vs existing {present_identifiers['inchikey']}")
                        logger.error(f"Please check the ModelSEED ID ({ID}) and the present inchi ({present_identifiers['inchikey']}) for {metabolite.id}.")
                        continue
                    else:
                        # If the InChIKey match, we can ignore it as it exists already
                        continue

                query = inchikey2modelseed(modelseed_result.inchi_key)
                identifiersORG_validation = validate_id(f"seed.compound:{query}")

                if query == ID and identifiersORG_validation:
                    logger.info(
                        f"Found InChIKey {modelseed_result.inchi_key} in ModelSEED for {metabolite.id}"
                    )

                    verified_identifiers = add2dict_unique(
                        key="inchikey",
                        value=modelseed_result.inchi_key,
                        directory=verified_identifiers
                    )


    # Add verified identifiers to the metabolite's annotation
    # ToDo use set?
    for key, values in verified_identifiers.items():
        if key == "smiles" and not settings.add_smiles_as_cross_reference:
            continue

        if isinstance(values, list):
            logger.error(
                f"Multiple entries found for {key} in verified identifiers: {values}. "
                f"But {key} should only have one entry. None of those entries will be added."
                f"You can check where those entries originated in the log using the INFO level."
            )
            continue

        metabolite.annotation = add2dict_unique(
            key=key,
            value=value,
            dictionary=metabolite.annotation
        )



def add_crossreferences(object: Union[Model,Reaction,Metabolite], consider_subobjects: bool = True):
    """
    Add cross-references to the given object. This function works with a cobrapy
    Model, Reaction, or Metabolite.

    Args:
        object (Union[Model, Reaction, Metabolite]): The COBRApy object to which
            cross-references should be added.
    """
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