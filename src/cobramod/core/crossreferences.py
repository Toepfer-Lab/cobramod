from functools import lru_cache
from typing import Set

import requests
from cobra import Model, Reaction, Metabolite
from cobra.core import Group
from requests import HTTPError

from cobramod.test import textbook

@lru_cache(maxsize=512)
def get_namespace_id(database:str) -> int:

    url ="https://registry.api.identifiers.org/restApi/namespaces/search/findByPrefix"
    payload ={
        "prefix":database
    }

    response = requests.get(url=url,params=payload)
    response.raise_for_status()
    response_json = response.json()

    potential_namespaces:str = response_json["_links"]["namespace"]["href"]

    if len(potential_namespaces) > 1:
        pass # raise uncertain

    id = potential_namespaces[potential_namespaces.rindex('/')+1:]
    return id


@lru_cache(maxsize=512)
def namespaceId2providerCode(id:int)-> [str]:

    url ="https://registry.api.identifiers.org/restApi/resources/search/findAllByNamespaceId?id=1691"
    payload = {
        "id":id
    }

    response = requests.get(url=url, params=payload)
    response.raise_for_status()
    response_json = response.json()
    providers = response_json["_embedded"]["resources"]

    provider_codes : [str] = []

    for provider in providers:
        provider_codes.append(provider["providerCode"])


    return  provider_codes


def get_crossreferences(sort, querys) -> Set[str]:

    query_list: str

    if isinstance(querys,list):
        query_list = " ".join(querys)
    else:
        query_list = querys

    url = "https://www.metanetx.org/cgi-bin/mnxweb/id-mapper"
    data = {
        "query_list":query_list,
        "query_index":sort,
        "output_format":"json"
    }


    response = requests.post(url=url,data=data)
    response.raise_for_status()

    response_json = response.json()

    #print(response_json)

    crossreferences = set()

    for query in querys:
        try:
            xrefs = response_json[query]["xrefs"]
            crossreferences.update(xrefs)
        except KeyError:
            pass

    return crossreferences


def add_crossreferences(object):

    if isinstance(object,Model):
        for reaction in object.reactions:
            add_crossreferences(reaction)

        for metabolite in object.metabolites:
            add_crossreferences(metabolite)

    elif isinstance(object,Group):

        for member in object.members:
            add_crossreferences(member)

    else:
        sort = None

        if isinstance(object,Reaction):
            sort = "reac"

        elif isinstance(object, Metabolite):
            sort = "chem"

        if sort is not None:

            ids = [object.id]
            xrefs = object.annotation

            for key, value in object.annotation.items():

                if isinstance(value,list):
                    for id in value:
                        ids.append(key+":"+id)
                else:
                    ids.append(key+":"+value)

            potential_xrefs = get_crossreferences(sort,ids)
            for potential_xref in potential_xrefs:
                database, id = potential_xref.split(":")

                try:
                    namespaceID = get_namespace_id(database=database)
                    providerCodes = namespaceId2providerCode(namespaceID)

                    for providerCode in providerCodes:
                        if providerCode not in xrefs.keys():
                            xrefs[providerCode] = potential_xref
                        else:
                            database_id = xrefs[providerCode]
                            if isinstance(database_id,str):
                                if database_id != potential_xref:
                                    xrefs[providerCode] = [database_id,potential_xref]
                            elif isinstance(database_id,list):
                                if potential_xref not in database_id:
                                    database_id.append(potential_xref)
                                    xrefs[providerCode] = database_id
                            else:
                                pass # raise warning unknown data in dict for providercode
                except HTTPError:
                    pass

            object.annotation = xrefs
