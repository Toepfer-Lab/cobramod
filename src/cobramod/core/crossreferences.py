import requests
from cobra import Model, Reaction, Metabolite
from cobra.core import Group

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

get_namespace_id("metacyc.compound")


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

namespaceId2providerCode(1691)

def get_crossreferences(sort, querys, xrefs = {}):

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

    try:

        response = requests.post(url=url,data=data)
        response.raise_for_status()

        response_json = response.json()

        print(response_json)

        for query in querys:
            try:
                buffer = response_json[query]["xrefs"]
                for ref in buffer:
                    key, value = ref.split(":")

                    if key in xrefs:
                        if isinstance(xrefs[key],str):
                            xrefs[key] = [xrefs[key],value]
                        else:
                            if value not in xrefs[key]:
                                xrefs[key].append(value)

                    else:
                        xrefs[key] = value
            except KeyError:
                pass
    except requests.exceptions.ConnectionError:
        pass

    return xrefs

get_crossreferences("chem",["metacyc.compound:ACETALD","metacycM:ACETALD"])

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

            for key, value in object.annotation.items():

                if isinstance(value,list):
                    for id in value:
                        ids.append(key+":"+id)
                else:
                    ids.append(key+":"+value)

            new_annotations = get_crossreferences(sort,ids,object.annotation)
            object.annotation = new_annotations
