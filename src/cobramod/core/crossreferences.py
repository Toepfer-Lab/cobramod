import re
from functools import lru_cache
from pathlib import Path
from typing import Set, Union, List, Any, Literal, Optional

import pandas as pd
import requests
from cobra import Model, Reaction, Metabolite
from cobra.core import Group
from pandas.core.interchange.dataframe_protocol import DataFrame
from requests import HTTPError
from tqdm import tqdm

from cobramod.dbwalker.identifiersORG import Validator
from cobramod.debug import debug_log


def findXRefs(id: str, cache: DataFrame) -> Optional[Union[List[str], str]]:
    """ """

    value = cache.loc[cache["ID"] == id]["XRefs"]

    if len(value) > 0:
        element = value.iloc[0]

        # single element return the containing str
        if len(element) == 1:
            return element[0]

        # return list of str
        return element

    return None


def inchikey2pubchem_cid(
    inchikey: Union[str, List[str]], directory: Path
) -> Union[str, List[str]]:
    """
    This function returns the corresponding PubChem compound ID for an
    InChIKey. A local cache is used, which is located in the specified
    directory under the folder XREF.
    Args:
        inchikey: The InChIKey for which the PubChem compound ID is to be
            searched.
        directory: The directory for storing the data. This is where
            the cache is stored in a folder called XRef.

    Returns:
        The PubChem compound ID(s) found as a string if it is one or as a
            list of strings otherwise.
    """
    if isinstance(inchikey, list):
        result: List[str] = []
        for key in inchikey:
            # following can return strings
            result.append(inchikey2pubchem_cid(key, directory))  # type: ignore
        return result

    cache = load_cache_from_disk("pubchem", directory)

    value = findXRefs(id=inchikey, cache=cache)

    if value is not None:
        return value

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        + inchikey
        + "/cids/txt"
    )

    response = requests.get(url=url)
    response.raise_for_status()

    value = response.text.rstrip()

    # PubChem returns one cid per line if multiple valid
    if "\n" in value:
        value = value.split("\n")

    # cache = cache.append({"ID": inchikey, "XRefs": value}, ignore_index=True)
    cache = pd.concat(
        [cache, pd.DataFrame([{"ID": inchikey, "XRefs": value}])],
        ignore_index=True,
    )
    path = directory / "XRef" / str("pubchem" + ".feather")
    path.parent.mkdir(parents=True, exist_ok=True)

    # create consistent content here list of str instead of lists[str] & str
    cache["XRefs"] = cache["XRefs"].apply(
        lambda x: [x] if isinstance(x, str) else x
    )

    cache.to_feather(path)
    load_cache_from_disk.cache_clear()
    return value


@lru_cache(maxsize=10)
def load_cache_from_disk(sort: str, directory: Path) -> pd.DataFrame:
    """
    The function loads the locally stored cache and returns it as a
    pandas DataFrame. If there is no cache file yet an empty DataFrame
    with columns 'ID' and 'XRefs' is returned.
    Args:
        sort: The cache name.
        directory: The directory for storing the data. This is where
            the cache is stored in a folder called XRef.

    """

    try:
        df = pd.read_feather(directory / "XRef" / str(sort + ".feather"))
    except FileNotFoundError:
        df = pd.DataFrame({"ID": [], "XRefs": []})

    return df


def get_crossreferences(  # noqa: C901
    sort: str, querys: Union[str, List[str]], directory: Path
) -> Set[str]:
    """
    Searches for IDs, other IDs from other databases. MetaNetX is
    used for this purpose. Results are stored locally in the specified
    directory and if they are found in it, they are loaded from it.
    Args:
        sort: Type of IDs' possible specifications are "chem" for
            metabolites or "reac" for reactions.
        querys: The IDs of a metabolite or a reaction. Can be either a string
            of the form "database:identifier" or a list of such strings.
            The list should only consist of identifiers for an object,
            as these will be merged.
        directory: The directory for storing the data. This is where
            the cache is stored in a folder called XRef.

    Returns:
        All references retrieved as a set of strings of the
        structure: "database:ID".
    """

    query_list: str
    cache = load_cache_from_disk(sort, directory)
    crossreferences = set()

    if isinstance(querys, list):
        all_querys = querys.copy()
        for query in querys:
            result = cache.loc[cache["ID"] == query]["XRefs"]

            result = findXRefs(id=query, cache=cache)

            if result is not None:
                all_querys.remove(query)
                crossreferences.update(set(result))

        if len(all_querys) == 0:
            return crossreferences
        query_list = " ".join(all_querys)
    else:
        all_querys = [querys]

        result = findXRefs(id=querys, cache=cache)
        if result is not None:
            return set(result)
        query_list = querys

    url = "https://www.metanetx.org/cgi-bin/mnxweb/id-mapper"
    data = {
        "query_list": query_list,
        "query_index": sort,
        "output_format": "json",
    }

    response = requests.post(url=url, data=data)
    try:
        response.raise_for_status()
    except HTTPError as error:
        if error.response.status_code == 413:
            size = len(querys)
            half = round(size / 2)
            chunk1 = querys[1:half]
            chunk2 = querys[half + 1 : size]

            result1 = get_crossreferences(sort, chunk1, directory)
            result2 = get_crossreferences(sort, chunk2, directory)

            return set.union(result1, result2)
        else:
            return set()
    response_json = response.json()

    for query in all_querys:
        try:
            answer = response_json[query]
        except KeyError:
            continue

        xrefs: List[str] = []
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

        replace = {"chem": "metanetx.chemical", "reac": "metanetx.reaction"}

        try:
            mnx_id = answer["mnx_id"]
            xrefs.append(replace[sort] + ":" + mnx_id)
        except KeyError:
            pass

        cache = pd.concat(
            [cache, pd.DataFrame([{"ID": query, "XRefs": xrefs}])],
            ignore_index=True,
        )
        # cache = cache.append({"ID": query, "XRefs": xrefs}, ignore_index=True)
        path = directory / "XRef" / str(sort + ".feather")
        path.parent.mkdir(parents=True, exist_ok=True)

        cache["XRefs"] = cache["XRefs"].apply(
            lambda x: [x] if isinstance(x, str) else x
        )
        cache.to_feather(path)
        load_cache_from_disk.cache_clear()
        crossreferences.update(xrefs)

    return crossreferences


def metanetx2ec(
    id: Union[str, List[str]],
    directory: Path,
    include_metanetx_specific_ec: bool = False,
) -> Union[str, List[str]]:
    """
    Returns the corresponding EC number for a specific MetaNetX ID.
    Args:
        id: The MetaNetX ID.
        directory: The directory for storing the data. This is where
            the cache is stored in a folder called XRef.
        include_metanetx_specific_ec: Determines whether MetaNetX specific
            EC numbers should be taken over. These are generally not found in
            other databases. The default value is False.
    Raises
        KeyError: If no EC number can be assigned to this ID.
    Returns:
        All found EC numbers.
    """
    if isinstance(id, List):
        result = set()
        for single_id in id:
            try:
                ec_numbers = metanetx2ec(
                    single_id, directory, include_metanetx_specific_ec
                )
            except KeyError:
                continue

            if isinstance(ec_numbers, str):
                ec_numbers = [ec_numbers]
            result.update(ec_numbers)

        size = len(result)
        if size == 0:
            raise KeyError
        elif size == 1:
            return result.pop()
        else:
            return list(result)

    data = get_reac_prop_with_ec(directory)
    found_ec_numbers = data.loc[data["ID"] == id, "classifs"]
    pattern = re.compile(
        r"^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\."
        r"\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$"
    )

    if len(found_ec_numbers) == 0:
        raise KeyError

    found_ec_numbers = found_ec_numbers.iloc[0]

    if ";" in found_ec_numbers:
        found_ec_numbers = found_ec_numbers.split(";")

        for found_ec_number in found_ec_numbers:
            # id can only be a string at this point
            if not include_metanetx_specific_ec and not pattern.match(
                found_ec_number
            ):
                found_ec_numbers.remove(found_ec_number)

        if len(found_ec_numbers) == 0:
            raise KeyError

    else:
        if not include_metanetx_specific_ec and not pattern.match(
            found_ec_numbers
        ):
            raise KeyError

    return found_ec_numbers


@lru_cache(maxsize=1)
def get_reac_prop_with_ec(directory: Path) -> pd.DataFrame:
    """
    This function loads the file reac_prop from MetaNetX and stores it in
    memory using the lru_cache. The returned DataFrame contains only the
    rows that have an EC number.

    Args:
        directory: The directory for storing the data. This is where
            the cache is stored in a folder called XRef.
    """

    path = directory / "XRef" / str("reac_prop" + ".feather")

    try:
        return pd.read_feather(path)
    except FileNotFoundError:
        pass

    url = "https://www.metanetx.org/ftp/latest/reac_prop.tsv"

    df = pd.read_csv(
        url,
        sep="\t",
        comment="#",
        names=[
            "ID",
            "mnx_equation",
            "reference",
            "classifs",
            "is_balanced",
            "is_transport",
        ],
    )

    df = df[~df["classifs"].isna()].reset_index(drop=True)

    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_feather(path)

    return df


def add2dict_unique(
    key, value: Union[Any, List[Any]], dictionary: dict
) -> dict:
    """
    Adds key-value pairs to a dictionary. It is expected that the values
    are either single objects or in a list. If the key does not yet exist
    in the dictionary, the value is added directly with it. If the key
    already exists, the values are joined together like this. That a list
    is assigned to the key, which contains unique values. If only one value
    remains, it will be added directly and not in a list.
    Args:
        key: The key to be used.
        value: The Values to be added. This can be a single object or
            objects in a list.
        dictionary: The dictionary to which the objects should be added.

    Returns:
        The new dictionary.
    """

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


def add_crossreferences(  # noqa: C901
    object: Union[Model, Group, Reaction, Metabolite],
    directory: Union[Path, str],
    consider_sub_elements: bool = True,
    include_metanetx_specific_ec: bool = False,
    validate: Literal["all", "only_new", "none"] = "only_new",
    use_metanetx: bool = False,
) -> None:
    """Extends the passed object by cross-references. Here, only the
    cross-references of reactions or metabolites are expanded. There
    must be at least one MetaNetX resolvable identifier in the annotation
    of the object in order to get as many cross-references as possible.
    Overall, an attempt is made to add all cross-references required by Memote.
    The only exception are reactom reactions. These are not added due to the
    current lack of query capabilities on the part of reactom. (Depending on
    the number of objects this function can take some time).

    Args:
        object: The CobraPy object to be extended.
        directory: The directory for storing the data. This is where
            the cache is stored in a folder called XRef.
        consider_sub_elements: Specifies whether additional cross-references
            should also be added to the subelements. For example, you can
            specify whether only the reaction or also its metabolites
            should be expanded.
        include_metanetx_specific_ec: Determines whether MetaNetX specific
            EC numbers should be taken over. These are generally not found in
            other databases. Furthermore, this could result in non-existing
            Brenda IDs being created. The default value is False.
        validate:

    """

    validate_all = False
    validate_new = False

    if validate == "all":
        validate_all = True
    elif validate == "only_new":
        validate_new = True
    elif validate == "none":
        pass
    else:
        raise ValueError("Invalid value for validate")

    if isinstance(directory, str):
        directory = Path(directory).absolute()

    if isinstance(object, Model):
        if not consider_sub_elements:
            debug_log.warn(
                "A model was passed but no sub-elements are to be annotated. "
                "Nothing changes."
            )
            return
        size = len(object.reactions) + len(object.metabolites)
        debug_log.debug(
            f"A model was passed. Trying to find further references for "
            f"{size} elements."
        )
        with tqdm(total=size) as pbar:
            for reaction in object.reactions:
                add_crossreferences(
                    reaction, directory=directory, consider_sub_elements=False
                )
                pbar.update(1)

            for metabolite in object.metabolites:
                add_crossreferences(metabolite, directory)
                pbar.update(1)

    elif isinstance(object, Group):
        if not consider_sub_elements:
            debug_log.warn(
                "A group was passed but no sub-elements are to be annotated. "
                "Nothing changes."
            )
            return

        debug_log.debug(
            f"A group was passed. Trying to find further references for "
            f"{len(object.members)} elements."
        )
        for member in tqdm(object.members):
            add_crossreferences(
                member, directory=directory, consider_sub_elements=False
            )

    else:
        sort = None

        if isinstance(object, Reaction):
            if consider_sub_elements:
                for metabolite in object.metabolites:
                    add_crossreferences(metabolite, directory)

            sort = "reac"

        elif isinstance(object, Metabolite):
            sort = "chem"

        if sort is None:
            debug_log.error(
                "An object was passed that is not of type Model, Group, "
                "Reaction or Metabolite! Check the 'object' value."
            )
            raise ValueError

        if object.id is not None:
            id = object.id
            ids = [object.id]
        else:
            id = "Undefined"
            ids = []
        xrefs = object.annotation
        original_xrefs = object.annotation.copy()
        new_xrefs: dict[str, Union[list[str], str]] = {}
        size = len(xrefs)
        total_found = 0

        # metanetx
        # metanetx can only be disabled for reactions
        if use_metanetx is True or isinstance(object, Metabolite):
            for key, value in object.annotation.items():
                if isinstance(value, list):
                    for id in value:
                        ids.append(key + ":" + id)
                else:
                    ids.append(key + ":" + value)

            potential_xrefs = get_crossreferences(sort, ids, directory)

            for potential_xref in potential_xrefs:
                try:
                    prefix, new_id = potential_xref.split(":")
                except ValueError:
                    continue
                prefix = prefix.lower()

                # MetaNetX does not provide correct CHEBI IDs. This is because
                # they lack the correct prefix. "chebi:CHEBI:0000" is correct
                # but "CHEBI:0000" is delivered. Therefore the following if
                # statement is used. This will also work if this condition is
                # fixed. (But then it is redundant)
                if prefix == "chebi" and "CHEBI:" not in new_id:
                    new_id = "CHEBI:" + new_id

                # obsolete ids are ignored
                elif prefix == "deprecated":
                    continue

                elif prefix == "mnx":
                    total_found += 1

                    if sort == "reac":
                        prefix_ = "metanetx.reaction"
                        xrefs = add2dict_unique(prefix_, new_id, xrefs)
                    else:
                        prefix_ = "metanetx.chemical"
                        xrefs = add2dict_unique(prefix_, new_id, xrefs)

                else:
                    total_found += 1
                    xrefs = add2dict_unique(prefix, new_id, xrefs)

                if validate_new:
                    new_xrefs = add2dict_unique(prefix, new_id, new_xrefs)

        # pubchem.compound
        try:
            inchikey = xrefs["inchikey"]
            pubchem_compound = inchikey2pubchem_cid(inchikey, directory)
        except (KeyError, HTTPError):
            pass
        else:
            xrefs = add2dict_unique("pubchem.compound", pubchem_compound, xrefs)

            if validate_new:
                new_xrefs = add2dict_unique(
                    "pubchem.compound", pubchem_compound, new_xrefs
                )

            if isinstance(pubchem_compound, str):
                pubchem_compound = [pubchem_compound]

            total_found += len(pubchem_compound)

        # ec_number
        try:
            metanet_xid = xrefs["metanetx.reaction"]
            ec_number = metanetx2ec(
                metanet_xid,
                directory,
                include_metanetx_specific_ec=include_metanetx_specific_ec,
            )
        except KeyError:
            pass
        else:
            xrefs = add2dict_unique("ec-code", ec_number, xrefs)

            if validate_new:
                new_xrefs = add2dict_unique("ec-code", ec_number, new_xrefs)

            if isinstance(ec_number, str):
                ec_number = [ec_number]

            total_found += len(ec_number)

        # brenda
        try:
            ec_number = xrefs["ec-code"]
        except KeyError:
            pass
        else:
            xrefs = add2dict_unique("brenda", ec_number, xrefs)
            if validate_new:
                new_xrefs = add2dict_unique("brenda", ec_number, new_xrefs)

            if isinstance(ec_number, str):
                ec_number = [ec_number]

            total_found += len(ec_number)

        total_added = len(xrefs) - size

        debug_log.debug(
            f'For the object with the ID {id} a total of "{total_found}" '
            f"references were found. Of these, {total_added} "
            f"were missing and have been added."
        )

        if validate_all:
            object.annotation = validate_id_dict(xrefs=xrefs)

        elif validate_new:
            new_xrefs = validate_id_dict(xrefs=new_xrefs)

            for prefix, identifier in original_xrefs.items():
                xrefs = add2dict_unique(prefix, identifier, new_xrefs)

            object.annotation = xrefs

        else:
            object.annotation = xrefs

validator = Validator()

def validate_id_dict(
    xrefs: dict[str, Union[str, list[str]]],
) -> dict[str, Union[str, list[str]]]:
    validated_xrefs: dict[str, Union[str, list[str]]] = {}

    for prefix, identifiers in xrefs.items():
        if isinstance(identifiers, str):
            if validator.validate_id_pattern(prefix=prefix, identifier=identifiers):
                validated_xrefs[prefix] = identifiers
        else:
            valid_list = []
            # Has to be list of identifiers
            for identifier in identifiers:
                if validator.validate_id_pattern(prefix=prefix, identifier=identifier):
                    valid_list.append(identifier)

            # if at least one entry remains add those either as list
            # if more than one or as string if exactly one
            if len(valid_list) > 1:
                validated_xrefs[prefix] = valid_list
            elif len(valid_list) == 1:
                validated_xrefs[prefix] = valid_list[0]

    print(validated_xrefs)
    return validated_xrefs

def validate_id(identifier: str) -> bool:
    debug_log.debug(f"Checking validity for identifier {identifier}")
    if ":" not in identifier:
        return False
    # Chebi indetifier look like this: chebi:CHEBI:123456
    # => for the look-up we need to remove one chebi
    # both URLS need the capital variant
    elif "CHEBI" in identifier:
        provider, prefix, ID = identifier.split(":")

        # capitalisation of provider is irrelevant
        if provider.lower() != "chebi":
            return False

        actual_id = f"{prefix}:{ID}"

    # biocyc identifiers look like this: biocyc:orgid:ID
    elif "biocyc" in identifier:
        provider, prefix, ID = identifier.split(":")
        ID = f"{prefix}:{ID}"
        prefix = f"{provider}"
        actual_id = ID
    else:
        prefix, ID = identifier.split(":")
        actual_id = ID

    url_identifiers = f"https://resolver.api.identifiers.org/{prefix}:{ID}"

    response = requests.get(url_identifiers)

    try:
        response.raise_for_status()
    except HTTPError as e:
        debug_log.debug(
            f"Validation request results in HTTP error: {e}"
        )
        return False


    answer = response.json()
    valid = True if answer["errorMessage"] is None else False

    if valid:
        debug_log.debug(
            f"Identifier ({identifier}) prefix & structure is verified: {valid}"
        )
    else:
        debug_log.debug(
            f"Validation request (identifiers.org) results in error message: {answer['errorMessage']}"
        )

    url_in_sbml = f"https://identifiers.org/{prefix}/{actual_id}"

    #response = requests.get(url_in_sbml)

    #try:
    #    response.raise_for_status()
    #except HTTPError as e:
    #    debug_log.debug(
    #        f"Validation request (Database) results in HTTP error: {e}"
    #    )
    #    return False

    return valid
