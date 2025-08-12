import logging

import requests
from typing import Optional, Tuple, Literal

import zeep
from zeep import Client

from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.Chebi")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# Initialize the ChEBI SOAP client
wsdl_url = "https://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl"
namespace = "https://www.ebi.ac.uk/webservices/chebi"
client = Client(wsdl=wsdl_url)

def getChebiID(identifier: str, identifier_type: Literal['SMILES', 'InChI', 'InChIKey']) -> str:
    """
    Get ChEBI ID for a given identifier.

    Args:
        identifier: The chemical identifier (SMILES, InChI, or InChIKey)
        identifier_type: Type of identifier - must be one of: 'smiles', 'inchi', 'inchikey'

    Returns:
        ChEBI ID
    """

    # inchi and inchikey use the same searchCategory in ChEBI
    if identifier_type == 'InChIKey':
        identifier_type = 'InChI'

    logger.info(f"Requesting ChEBI ID for {identifier} of type {identifier_type}.")
    response = client.service.getLiteEntity(
        search=identifier,
        searchCategory=identifier_type,
        maximumResults=50,
        stars="ALL"
    )

    if response is None:
        logger.warning(f"There is no entry for {identifier} of type {identifier_type}.")
        return response

    logger.debug(f"Got following as result: {response}.")

    if len(response) > 1:
        logger.warning(f"Multiple results found for {identifier} of type {identifier_type}. Returning first result.")

    response = response[0]
    chebi_id = response.chebiId

    return chebi_id

def get_chebi_compound_info(chebi_id: str) -> GenerellIdentifiers:
    """
    Retrieve InChI, InChIKey, and SMILES for a compound from ChEBI database.
    
    Args:
        chebi_id: ChEBI identifier (e.g., 'CHEBI:15377' or '15377')
    
    Returns:
        Tuple of (inchi, inchikey, smiles) or None if not found
    """
    

    response = client.service.getCompleteEntity(
        chebiId=chebi_id,
        )
    response = zeep.helpers.serialize_object(response)

    inchi = response.get('inchi')
    if inchi.startswith("InChI="):
        inchi = inchi[6:]

    inchikey = response.get('inchiKey')
    if inchikey.startswith("InChIKey="):
        inchikey = inchikey[9:]

    smiles = response.get('smiles', None)

    return GenerellIdentifiers(inchi = inchi, inchi_key = inchikey, smiles = smiles)