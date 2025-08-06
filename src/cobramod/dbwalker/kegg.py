import logging

import requests

logger = logging.getLogger("cobramod.DBWalker.Kegg")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

def get_kegg_id_from_cid(cid):
    """
    Given a PubChem CID, return the corresponding KEGG compound ID (e.g., 'C00022').
    """
    url = f"https://rest.kegg.jp/conv/compound/pubchem:{cid}"
    response = requests.get(url)
    
    if response.status_code != 200:
        logger.error(f"Error ({response.status_code}) getting KEGG compound ID from {cid}")
        return None

    lines = response.text.strip().split("\n")
    if not lines:
        return None

    # Parse the first matching line
    try:
        kegg_id = lines[0].split("\t")[1]
    except IndexError:
        return None

    if ":" in kegg_id:
        kegg_id = kegg_id.split(":")[1]

    kegg_id = f"kegg.compound:{kegg_id}"

    return kegg_id

def get_cid_from_kegg_id(kegg_id):
    """
    Given a KEGG compound ID (e.g., 'C00022' or 'kegg.compound:C00022'),
    return the corresponding PubChem CID.
    """
    # Extract the compound ID if it's in full format
    if ":" in kegg_id:
        compound_id = kegg_id.split(":")[1]
    else:
        compound_id = kegg_id

    url = f"https://rest.kegg.jp/conv/pubchem/cpd:{compound_id}"
    response = requests.get(url)

    if response.status_code != 200:
        logger.error(f"Error ({response.status_code}) while getting CID for {compound_id}")
        return None

    lines = response.text.strip().split("\n")
    if not lines:
        logger.info(f"No entry found for KEGG ID {kegg_id}")
        return None

    # Parse the first matching line
    pubchem_entry = lines[0].split("\t")[1]

    if ":" in pubchem_entry:
        cid = pubchem_entry.split(":")[1]
    else:
        cid = pubchem_entry

    return f"pubchem.compound:{cid}"