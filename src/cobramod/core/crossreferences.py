import time
from functools import lru_cache
from typing import Set, Union, List

import pandas as pd
import requests
import tqdm
from cobra import Model, Reaction, Metabolite
from cobra.core import Group
from requests import HTTPError


@lru_cache(maxsize=512)
def get_namespace_id(database: str) -> int:
    url = "https://registry.api.identifiers.org/restApi/namespaces/search/findByPrefix"
    payload = {
        "prefix": database
    }

    response = requests.get(url=url, params=payload)
    response.raise_for_status()
    response_json = response.json()

    potential_namespaces: str = response_json["_links"]["namespace"]["href"]

    if len(potential_namespaces) > 1:
        pass  # raise uncertain

    id = potential_namespaces[potential_namespaces.rindex('/') + 1:]
    return id


def inchikey2pubchem_cid(inchikey: str) -> Union[str, List[str]]:
    if isinstance(inchikey, list):
        result: [int] = []
        for key in inchikey:
            result.append(inchikey2pubchem_cid(key))
        return result

    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/" \
          + inchikey + \
          "/cids/txt"

    response = requests.get(url=url)
    response.raise_for_status()
    return response.text


@lru_cache(maxsize=512)
def namespaceid2provider_code(id: int) -> [str]:
    url = "https://registry.api.identifiers.org/restApi/resources/search/findAllByNamespaceId?id=1691"
    payload = {
        "id": id
    }

    response = requests.get(url=url, params=payload)
    response.raise_for_status()
    response_json = response.json()
    providers = response_json["_embedded"]["resources"]

    provider_codes: [str] = []

    for provider in providers:
        provider_codes.append(provider["providerCode"])

    return provider_codes


def get_crossreferences(sort, querys) -> Set[str]:
    query_list: str

    if isinstance(querys, list):
        query_list = " ".join(querys)
    else:
        query_list = querys

    url = "https://www.metanetx.org/cgi-bin/mnxweb/id-mapper"
    data = {
        "query_list": query_list,
        "query_index": sort,
        "output_format": "json"
    }

    response = requests.post(url=url, data=data)
    response.raise_for_status()

    response_json = response.json()

    crossreferences = set()

    for query in querys:
        try:
            answer = response_json[query]
        except KeyError:
            continue

        xrefs: [str] = []
        try:
            xrefs = answer["xrefs"]
        except KeyError:
            pass
        try:
            xrefs.append("inchi:" + answer["InChI"])
        except KeyError:
            pass
        try:
            inchikey = answer["InChIkey"]
            xrefs.append("inchikey:" + inchikey)
        except KeyError:
            pass

        crossreferences.update(xrefs)

    return crossreferences


def metanetx2ec(id: str) -> Union[str, List[str]]:
    data = get_reac_prop_with_ec()
    result = data.loc[data['ID'] == id, "classifs"]

    if len(result) == 0:
        raise KeyError

    result = result.iloc[0]

    if ";" in result:
        result = result.split(";")

    return result


@lru_cache(maxsize=1)
def get_reac_prop_with_ec() -> pd.DataFrame:
    url = "https://www.metanetx.org/ftp/latest/reac_prop.tsv"

    dataFrame = pd.read_csv(url,
                            sep="\t",
                            comment="#",
                            names=["ID",
                                   "mnx_equation",
                                   "reference",
                                   "classifs",
                                   "is_balanced",
                                   "is_transport"])

    dataFrame = dataFrame.loc[dataFrame["classifs"].notna()]
    return dataFrame


def add2dict_unique(key, value, dictionary):
    if key not in dictionary.keys():
        dictionary[key] = value
    else:
        ref = dictionary[key]
        if isinstance(ref, str):
            ref = [ref]
        if isinstance(value, str):
            value = [value]

        unique_pub = list(set(ref + value))
        if len(unique_pub) == 1:
            dictionary[key] = unique_pub[0]
        else:
            dictionary[key] = unique_pub

    return dictionary


def add_crossreferences(object: Union[Model, Group, Reaction, Metabolite]) -> None:
    """Extends the passed object by cross references. Here, only the
    cross references of reactions or metabolites are expanded. There
    must be at least one MetaNetX resolvable identifier in the annotation
    of the object in order to get as many cross references as possible.
    Overall, an attempt is made to add all crossreferences required by Memote.
    The only exception are reactom reactions. These are not added due to the
    current lack of query capabilities on the part of reactom. (Depending on
    the number of objects this function can take some time).

    Args:
        object: The CobraPy object to be extended.

    """
    if isinstance(object, Model):
        size = len(object.reactions) + len(object.metabolites)
        with tqdm(total=size) as pbar:
            for reaction in object.reactions:
                add_crossreferences(reaction)
                pbar.update(1)

            for metabolite in object.metabolites:
                add_crossreferences(metabolite)
                pbar.update(1)

    elif isinstance(object, Group):

        for member in tqdm(object.members):
            add_crossreferences(member)

    else:
        sort = None

        if isinstance(object, Reaction):
            sort = "reac"

        elif isinstance(object, Metabolite):
            sort = "chem"

        if sort is not None:

            # metanetx

            ids = [object.id]
            xrefs = object.annotation

            for key, value in object.annotation.items():

                if isinstance(value, list):
                    for id in value:
                        ids.append(key + ":" + id)
                else:
                    ids.append(key + ":" + value)

            potential_xrefs = get_crossreferences(sort, ids)

            for potential_xref in potential_xrefs:

                try:
                    prefix, new_id = potential_xref.split(":")
                except ValueError:
                    continue
                prefix = prefix.lower()

                # MetaNetX does not provide correct CHEBI IDs. This is because
                # they lack the correct prefix. "CHEBI:CHEBI:0000" is correct
                # but "CHEBI:0000" is delivered. Therefore the following if
                # statement is used. This will also work if this condition is
                # fixed. (But then it is redundant)
                if prefix == "chebi" and not "CHEBI:" in new_id:
                    new_id = "CHEBI:" + new_id

                # obsolete ids are ignored
                elif prefix == "deprecated":
                    continue

                xrefs = add2dict_unique(prefix, new_id, xrefs)

            # pubchem.compound

            try:
                inchikey = xrefs["inchikey"]
                pubchem_compound = inchikey2pubchem_cid(inchikey)
            except (KeyError, HTTPError):
                pass
            else:

                xrefs = add2dict_unique("pubchem.compound", pubchem_compound, xrefs)

            # ec_number
            try:
                metanetXID = xrefs["metanetx.reaction"]
                ec_number = metanetx2ec(metanetXID)
            except KeyError:
                pass
            else:
                xrefs = add2dict_unique("ec-code", ec_number, xrefs)

            # brenda
            try:
                ec_number = xrefs["ec-code"]
            except KeyError:
                pass
            else:
                xrefs = add2dict_unique("brenda", ec_number, xrefs)

            object.annotation = xrefs
