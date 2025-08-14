import urllib

import requests
from typing import Optional, Dict, Any

import xml.etree.ElementTree as ET
import logging

from cobramod.dbwalker.dataclasses import GenerellIdentifiers
import cobramod

logger = logging.getLogger("cobramod.DBWalker.BioCyc")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

biocycTier1DB = ["Meta","ECOLI", "HUMAN", "ARA", "YEAST"]

def get_compound_info_by_biocyc_id(biocyc_id: str, db: str = "META") -> GenerellIdentifiers:
    """
    Get InChI, InChI Key, and SMILES for a compound using BioCyc ID.
    
    Args:
        biocyc_id: BioCyc compound ID (e.g., 'CPD-123')
        db: BioCyc database name (default: 'META' for MetaCyc)
    
    Returns:
        Dictionary with 'inchi', 'inchi_key', and 'smiles' keys
    """

    url = f"https://websvc.biocyc.org/getxml"
    
    params = {
            'id': f"{db}:{biocyc_id}"
    }

    # BioCyc cuts the SMILES string short otherwise
    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "application/xml"
    }

    logger.debug(f"Querying BioCyc for compound {biocyc_id} in database {db}")
    logger.debug(f"Request URL: {url}")
    logger.debug(f"Request parameters: {params}")

    result = GenerellIdentifiers()
    try:
        cSettings = cobramod.Settings()
        session = cSettings._biocycSession
        response = session.get(url, params=params,headers=headers, timeout=30)
        response.raise_for_status()

        if cSettings.autoOpenCloseBioCycSession:
            cSettings._closeBiocycSession()

        logger.debug(f"Response status code: {response.status_code}")

        # Check if BioCyc is asking for account creation
        if response.text and isinstance(response.text, str) and "<title>Create Account</title>" in response.text:
            logger.warning("BioCyc is requesting account creation.")
            if not cSettings.BioCycLoggedIn:
                logger.error("BioCyc account information are not provided can't login.")
            return result

        # Parse XML response
        root = ET.fromstring(response.content)

        
        # Extract SMILES from XML
        smiles_element = root.find(".//string[@title='smiles']")
        if smiles_element is not None:
            result.smiles = smiles_element.text.strip()

        # Extract InChI from XML
        inchi_element = root.find(".//inchi[@datatype='string']")
        if inchi_element is not None:
            inchi = inchi_element.text.strip()
            if inchi.startswith("InChI="):
                result.inchi = inchi[6:]
            else:
                result.inchi = inchi

        # Extract InChI Key from XML
        inchikey_element = root.find(".//inchi-key[@datatype='string']")
        if inchikey_element is not None:
            inchikey = inchikey_element.text.strip()

            if inchikey.startswith("InChIKey="):
                result.inchi_key = inchikey[9:]
            else:
                result.inchi_key = inchikey
        
    except requests.RequestException as e:
        logger.error(f"Error fetching data from BioCyc: {e}")
    except ET.ParseError as e:
        logger.error(f"Error parsing XML response: {e}")

    logger.debug(f"Final result: {result}")
    return result

def smiles2BioCyc(smiles: str, orgid: str = "META") -> str:
    """
    Convert SMILES to BioCyc identifiers.
    
    Args:
        smiles: SMILES string of the compound
        
    Returns:
        GenerellIdentifiers with BioCyc ID if found, otherwise empty
    """
    try:
        encoded_smiles = urllib.parse.quote(smiles, safe='')
        url = f"https://websvc.biocyc.org/{orgid}/smiles-search?smiles={encoded_smiles}&exact=T&fmt=json"
        
        logger.info(f"Requesting BioCyc for SMILES: {smiles}")
        logger.debug(f"Request URL: {url}")

        cSettings = cobramod.Settings()
        session = cSettings._biocycSession

        response = session.get(url, timeout=30)
        response.raise_for_status()

        cSettings._closeBiocycSession()
        
        data = response.json()

        print(data)
        biocycID = data["RESULTS"][0]["OBJECT-ID"]
        return biocycID
        
        
    except requests.RequestException as e:
        logger.error(f"Error fetching data from BioCyc: {e}")
        return None

def InChI2BioCyc(inchi: str, orgid: str = "META") -> Optional[str]:
    """
    Convert InChI to BioCyc identifiers.
    
    Args:
        inchi: InChI string of the compound
        
    Returns:
        GenerellIdentifiers with BioCyc ID if found, otherwise empty
    """
    try:
        url = f"https://websvc.biocyc.org/{orgid}/inchi-search?inchi={inchi}&exact=T&fmt=json"
        
        logger.info(f"Requesting BioCyc for InChI: {inchi}")

        cSettings = cobramod.Settings()
        session = cSettings._biocycSession

        response = session.get(url, timeout=30)
        response.raise_for_status()

        cSettings._closeBiocycSession()
        
        data = response.json()

        biocycID = data["RESULTS"][0]["OBJECT-ID"]
        return biocycID
        
        
    except requests.RequestException as e:
        logger.error(f"Error fetching data from BioCyc: {e}")
        return None

#def Identifiers2BioCyc(identifier: str) -> GenerellIdentifiers:
    