#!/usr/bin/env python3
from contextlib import suppress
from itertools import chain
from pathlib import Path
from typing import Generator, Union, Dict, NamedTuple

from requests import get

from cobramod.debug import debug_log
from cobramod.parsing.base import BaseParser
from cobramod.error import WrongParserError


class MetaboliteTuple(NamedTuple):
    """
    Simple named tuple with information of the metabolite.
    """

    identifier: str
    coefficient: float


def _get_keys(raw: str) -> Generator:
    """
    Returns as Generator the keys for KEGGs data. Input is raw text
    """
    lines = (line for line in raw.split(sep="\n"))
    for line in lines:
        segment = line.split(" ")
        if segment[0] == "///":
            continue
        if segment[0] != "":
            yield segment[0]


def _find_key(line: str, keys: set) -> Union[str, None]:
    """
    Return key if found in given line. Else, None is returned.
    """
    for key in keys:
        if key in line:
            return key
    return None


def _create_dict(raw: str) -> dict:
    """
    Formats most of the keys for Keggs data and returns a dictionary.
    """
    lines = (line for line in raw.split("\n"))
    keys = set(_get_keys(raw=raw))
    actual_key = str()
    kegg_dict: Dict[str, list] = {"": []}
    # TODO: find a proper way to parse it
    for line in lines:
        if "///" in line:
            break
        else:
            key = _find_key(line=line, keys=keys)
            if key is None:
                # Append to older key.
                kegg_dict[actual_key].append(line.strip().rstrip())
            else:
                actual_key = key
                line = line[line.find(key) + len(key) :].strip().rstrip()
                kegg_dict[actual_key] = [line]
    del kegg_dict[""]
    return kegg_dict


def _get_coefficient(line: str, SIDE: int, prefix: str) -> MetaboliteTuple:
    """
    Returns a NamedTuple with the metabolite identifier and coefficient
    that appears in given line.

    Args:
        line (str): string with information
        SIDE (int): Constant to multiply depending on the side
        prefix (str): prefix for the metabolite

    Returns:
        NamedTuple: tupple with identifier and coefficient
    """
    # " 4 H+ "
    line_list = line.strip().rstrip().split(" ")
    try:
        return MetaboliteTuple(
            identifier=prefix + line_list[1],
            coefficient=float(line_list[0]) * SIDE,
        )
    except IndexError:
        return MetaboliteTuple(
            identifier=prefix + line_list[0], coefficient=1.0 * SIDE
        )


def _build_dict_meta(line: str) -> dict:
    """
    Builds dictionary with metabolites identifiers and corresponding
    coeffients.
    """
    middle = line.find("=")
    temp_dict = dict()
    reactants = line[: middle - 1]
    products = line[middle + 2 :]
    for item in reactants.split(" + "):
        MetaTuple = _get_coefficient(line=item, SIDE=-1, prefix="l_")
        temp_dict[MetaTuple.identifier] = MetaTuple.coefficient
    for item in products.split(" + "):
        MetaTuple = _get_coefficient(line=item, SIDE=1, prefix="r_")
        temp_dict[MetaTuple.identifier] = MetaTuple.coefficient
    return temp_dict


def _get_reversibility(line: str) -> tuple:
    """
    Returns the bounds for reaction depending of the string
    """
    # FIXME: Direction depends also from extra keys
    middle = line.find("=")
    line = line[middle - 1 : middle + 2]
    if line == "<=>":
        bounds = (-1000, 1000)
    elif line == "==>":
        bounds = (0, 1000)
    elif line == "<==":
        bounds = (-1000, 0)
    else:
        raise Warning(f'"Reversibility for "{line}" could not be found.')
    return bounds


def _build_reference(data_dict: dict) -> Union[dict, None]:
    """
    From the dictionary with KEGG raw information, return a dictionary with
    where the keys are the names of cross-references and keys their
    identifiers.
    """
    try:
        return {
            item[0].strip().rstrip(): item[1].strip().rstrip()
            for item in [name.split(":") for name in data_dict["DBLINKS"]]
        }
    except KeyError:
        return None


def _p_reaction(kegg_dict: dict) -> dict:
    """
    Parses the KEGG dictionary and returns a dictionary with the most important
    attributes about the compound.
    """
    try:
        object_type = list(
            chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
        )[-1]
        if object_type != "Reaction":
            raise WrongParserError(
                "Given dictionary does not belong to a reaction."
            )
        return {
            "TYPE": object_type,
            "ENTRY": list(
                chain.from_iterable(
                    item.split() for item in kegg_dict["ENTRY"]
                )
            )[-2],
            "NAME": kegg_dict["NAME"][0],
            "DATABASE": "KEGG",
            "EQUATION": _build_dict_meta(line=kegg_dict["EQUATION"][0]),
            "BOUNDS": _get_reversibility(line=kegg_dict["EQUATION"][0]),
            "TRANSPORT": BaseParser._check_transport(
                data_dict=_build_dict_meta(line=kegg_dict["EQUATION"][0])
            ),
            "XREF": _build_reference(data_dict=kegg_dict),
        }
    except KeyError:
        raise WrongParserError


def _p_compound(kegg_dict: dict) -> dict:
    """
    Parses the KEGG dictionary and returns a dictionary with the most important
    attributes about the reaction.
    """
    try:
        object_type = list(
            chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
        )[-1]
        if object_type != "Compound":
            raise WrongParserError(
                "Given dictionary does not belong to a metabolite."
            )
        identifier = list(
            chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
        )[-2]
        try:
            if "(" in kegg_dict["FORMULA"][0]:
                formula = "X"
                debug_log.warning(
                    f'KEGG ID "{identifier}" could not find a chemical '
                    f'formula. Formula set to "X" and charge to 0'
                )
            else:
                formula = "".join(kegg_dict["FORMULA"]).strip().rstrip()
        except KeyError:
            # Some general compounds do not come with formulas
            formula = "X"
            debug_log.warning(
                f'KEGG ID "{identifier}" could not find a chemical formuala. '
                f'Formula set to "X" and charge to 0'
            )
        return {
            "TYPE": list(
                chain.from_iterable(
                    item.split() for item in kegg_dict["ENTRY"]
                )
            )[-1],
            "ENTRY": identifier,
            "NAME": kegg_dict["NAME"][0],
            "DATABASE": "KEGG",
            "FORMULA": formula,
            # FIXME: search information about charges
            "CHARGE": 0,
            "XREF": _build_reference(data_dict=kegg_dict),
        }
    except KeyError:
        raise WrongParserError


def _build_pathway(kegg_dict: dict) -> tuple:
    """
    Returns dictionary with sequences for a graph, where the key the prior
    reaction is and the value, its succesor; and a set with all participant
    reaction (vertex)
    """
    sequences = [item.split(sep=" ")[0] for item in kegg_dict["REACTION"]]
    whole_set = set(
        chain.from_iterable(line.split(sep=",") for line in sequences)
    )
    sequences = [line.split(sep=",")[0] for line in sequences]
    pathway = dict()
    for index, reaction in enumerate(iterable=sequences):
        with suppress(IndexError):
            # If index not found, then ignore it
            pathway[sequences[index + 1]] = reaction
    return pathway, whole_set


def _p_pathway(kegg_dict: dict) -> dict:
    """
    Parses the KEGG dictionary and returns a new dictionary with the
    information of the pathway.
    """
    try:
        if "Pathway" not in kegg_dict["ENTRY"][0]:
            raise WrongParserError(
                "Given dictionary does not belong to a pathway."
            )
        pathway, whole_set = _build_pathway(kegg_dict=kegg_dict)
        return {
            "TYPE": "Pathway",
            "ENTRY": list(
                chain.from_iterable(
                    item.split() for item in kegg_dict["ENTRY"]
                )
            )[0],
            "NAME": kegg_dict["NAME"][0],
            "DATABASE": "KEGG",
            "PATHWAY": pathway,
            "SET": whole_set,
            "XREF": _build_reference(data_dict=kegg_dict),
        }
    except KeyError:
        raise WrongParserError


def _p_enzyme(kegg_dict: dict):
    """
    If enzyme information is found, raise NotImplementedError.
    """
    try:
        if "Enzyme" in kegg_dict["ENTRY"][0]:
            raise NotImplementedError("Enzymes are currently not implemented")
        else:
            raise WrongParserError(
                "Given dictionary does not belong to an Enzyme."
            )
    except KeyError:
        raise WrongParserError


class KeggParser(BaseParser):
    @staticmethod
    def _parse(root: dict) -> dict:
        """
        Parses raw text from KEGG database and transforms it in a dictionary
        that can be used later as commom data type for cobramod. It identifies
        the type of object depending of it information.

        Args:
            root (dict): dictionary with information for an object

        Returns:
            dict: data from KEGG
        """
        for parse_method in (_p_enzyme, _p_compound, _p_reaction, _p_pathway):
            with suppress(WrongParserError):
                return parse_method(kegg_dict=root)
        raise NotImplementedError(
            "Given identifier could not be parsed properly. "
            "Contact maintainers."
        )

    @staticmethod
    def _retrieve_data(
        directory: Path,
        identifier: str,
        database: str,
        debug_level: int,
        **kwargs,
    ) -> dict:
        """
        Retrieves data from KEGG database and parses the most important
        attributes into a dictionary

        Args:
            directory (Path): Directory to store and retrieve local data.
            identifier (str): original identifier
            database (str): name of database.
        debug_level (int): Level of debugging. Read package logging
            for more info.

        Returns:
            dict: relevant data for given identifier
        """
        raw = _get_unformatted_kegg(directory=directory, identifier=identifier)
        debug_log.log(
            level=debug_level, msg=f'Data for "{identifier}" retrieved.'
        )
        return KeggParser._parse(root=raw)

    @staticmethod
    def _return_database(database: str) -> str:
        """
        Returns the name of the database. This method is used to compare with
        given database name. It will raise a warning if both names are not
        equal or belong to the list of proper names.
        """
        if database == "KEGG":
            return database
        else:
            raise WrongParserError

    @staticmethod
    def _read_file(filename: Path) -> dict:
        """
        Reads the given file a returns a KEGG dictionary with most important
        information from it.
        """
        try:
            with open(file=filename, mode="r") as f:
                unformatted_data = f.read()
            return _create_dict(raw=unformatted_data)
        except Exception:
            # TODO find exception type
            raise Warning("Wrong filetype")


def _get_unformatted_kegg(directory: Path, identifier: str) -> dict:
    """
    Retrieves and stores the data of given identifier for the KEGG database.

    Args:
        directory (Path): Path for data storage.
        identifier (str): official name for given object in KEGGs.

    Raises:
        NotImplementedError: If identifier contains a semicolon (:)
        Warning: If given identifier is not found in KEGGs.
        NotADirectoryError:  If parent directory is not found.

    Returns:
        dict: dictioanry from  KEGGs identifier with basic information
    """
    # NOTE: As KEGG only support xml for pathways, Metabolites and Reactions
    # have to be parsed differently. Only specific keys are necesary
    database = "KEGG"
    if ":" in identifier:
        raise NotImplementedError(
            "Identifiers with semicolons cannot be parse yet"
        )
    if directory.exists():
        data_dir = KeggParser._define_base_dir(
            directory=directory, database=database
        )
        filename = data_dir.joinpath(f"{identifier}.txt")
        debug_log.debug(f'Searching "{identifier}" in directory "{database}".')
        try:
            debug_log.debug(f"Identifier '{identifier}' found.")
            return KeggParser._read_file(filename=filename)
        except FileNotFoundError:
            debug_log.debug(
                f'"{identifier}" not found in directory "{database}".'
            )
            # Retrieve from URL
            url_text = f"http://rest.kegg.jp/get/{identifier}/"
            debug_log.debug(f"Searching in {url_text}")
            r = get(url_text)
            if r.status_code >= 400:
                msg = f'"{identifier}" not available in "{database}".'
                debug_log.error(msg)
                raise Warning(msg)
            else:
                unformatted_data = r.text
                if len(unformatted_data) == 0:
                    raise Warning(
                        f'Object "{identifier}" returned empty string'
                    )
                else:
                    debug_log.info(
                        f'Object "{identifier}" found in database. Saving in '
                        f'directory "{database}".'
                    )
                    with open(file=filename, mode="w+") as f:
                        f.write(unformatted_data)
                    return _create_dict(raw=unformatted_data)
    else:
        msg = "Directory not found"
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
