import requests
import logging
from typing import Optional

from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.ModelSeed")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def get_compound_info_by_modelseed_id(modelseed_id: str) -> GenerellIdentifiers:
    """
    Get InChI, InChI Key, and SMILES for a compound using ModelSEED ID.
    
    Args:
        modelseed_id: ModelSEED compound ID (e.g., 'cpd00131')
    
    Returns:
        GenerellIdentifiers object with chemical identifiers
    """
    url = (f"https://modelseed.org/solr/compounds/select"
           f"?q=id:{modelseed_id}"
           f"&fl=id,name,smiles,inchi,inchikey"
           f"&wt=json&rows=5"
           )

    logger.debug(f"Querying ModelSEED for compound {modelseed_id}")
    logger.debug(f"Request URL: {url}")

    result = GenerellIdentifiers()
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        logger.debug(f"Response status code: {response.status_code}")

        data = response.json()
        docs = data.get('response', {}).get('docs', [])
        
        if docs:
            doc = docs[0]
            
            # Extract SMILES
            if doc.get("smiles"):
                result.smiles = doc.get("smiles").strip()

            # Extract InChI
            if doc.get("inchi"):
                inchi = doc.get("inchi").strip()
                if inchi.startswith("InChI="):
                    result.inchi = inchi[6:]
                else:
                    result.inchi = inchi

            # Extract InChI Key
            inchikey = doc.get("inchikey") or doc.get("InChIKey")
            if inchikey:
                inchikey = inchikey.strip()
                if inchikey.startswith("InChIKey="):
                    result.inchi_key = inchikey[9:]
                else:
                    result.inchi_key = inchikey
                    
            logger.debug(f"Found compound data for {modelseed_id}")
        else:
            logger.warning(f"No compound found for ModelSEED ID: {modelseed_id}")
        
    except requests.RequestException as e:
        logger.error(f"Error fetching data from ModelSEED: {e}")
    except (KeyError, ValueError) as e:
        logger.error(f"Error parsing ModelSEED response: {e}")

    logger.debug(f"Final result: {result}")
    return result


def smiles2ModelSeed(smiles: str) -> Optional[str]:
    """
    Convert SMILES to ModelSEED identifier.
    
    Args:
        smiles: SMILES string of the compound
        
    Returns:
        ModelSEED ID if found, otherwise None
    """

    url = (f"https://modelseed.org/solr/compounds/select"
           f"?q=smiles:{smiles}"
           f"&fl=id,name,smiles,inchi,inchikey"
           f"&wt=json&rows=5"
           )
    
    logger.info(f"Requesting ModelSEED for SMILES: {smiles}")
    logger.debug(f"Request URL: {url}")

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        docs = data.get('response', {}).get('docs', [])
        
        if docs:
            modelseed_id = docs[0].get("id")
            logger.debug(f"Found ModelSEED ID: {modelseed_id}")
            return modelseed_id
        else:
            logger.warning(f"No ModelSEED ID found for SMILES: {smiles}")
            return None
        
    except requests.RequestException as e:
        logger.error(f"Error fetching data from ModelSEED: {e}")
        return None

def inchikey2modelseed(inchikey: str) -> Optional[str]:
    """
    Convert InChI Key to ModelSEED identifier.

    Args:
        inchikey: InChI Key string of the compound

    Returns:
        ModelSEED ID if found, otherwise None
    """

    url = (f"https://modelseed.org/solr/compounds/select?"
           f"q=inchikey:{inchikey}"
           f"&fl=id,name,smiles,inchi,inchikey"
           f"&wt=json"
           f"&rows=5"
           )

    logger.info(f"Requesting ModelSEED for InChI Key: {inchikey}")
    logger.debug(f"Request URL: {url}")

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        data = response.json()
        docs = data.get('response', {}).get('docs', [])

        if docs:
            modelseed_id = docs[0].get("id")
            logger.debug(f"Found ModelSEED ID: {modelseed_id}")
            return modelseed_id
        else:
            logger.warning(f"No ModelSEED ID found for InChI Key: {inchikey}")
            return None

    except requests.RequestException as e:
        logger.error(f"Error fetching data from ModelSEED: {e}")
        return None
