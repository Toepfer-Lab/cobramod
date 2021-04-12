#!/usr/bin/env python3
"""Data parsing for KEGG

This module handles the retrieval of data from KEGG into a local directory.
The posible type of data that can be download:

- Metabolite: Identifiers that start with the letter C, e.g C00001
- Reactions: Identifiers that start with the letter R, e.g R00001. The
gene information for reactions is also included if the specie is specified
- Module Pathways: Identifiers that start with the letter M, e.g M00001

Contact maintainers if other types should be added.

Important class of the module:
- KeggParser: Child of the abstract class
:func:`cobramod.parsing.base.BaseParser`.
"""
from contextlib import suppress
from itertools import chain, islice
from pathlib import Path
from typing import Generator, Union, Dict, NamedTuple
from warnings import warn

from requests import get, HTTPError

from cobramod.debug import debug_log
from cobramod.parsing.base import BaseParser
from cobramod.error import WrongParserError, FalseAbbreviation


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


def _get_ko(identifier: str) -> list:
    """
    Returns a list with the corresponding KO-entries for given identifier
    """
    url_text = f"http://rest.kegg.jp/link/ko/{identifier}"
    try:
        r = get(url_text)
        r.raise_for_status()
        return [
            single for single in r.text.split() if identifier not in single
        ]
    except HTTPError:
        raise HTTPError(
            f'No K-Identifiers could be retrieved for "{identifier}"'
        )


def _parse_ko(string: str, reaction: str, genome: str = None) -> list:
    """
    Parses given kegg text and returns with a list with the corresponding genes
    for given genome or if no genome is specified, the genes for the first 5
    genes. Argument reaction is just for logging information.
    """
    # TODO: add logging?
    # Retrive the KO-identifiers
    text = string.replace("\t", " ").splitlines()
    genes: Dict[str, list] = dict()
    # Separate the lines and append the gene identifier to the corresponding
    # genome
    for line in text:
        ko, single = line.split()
        specie, identifier = single.split(":")
        try:
            genes[specie].append(identifier)
        except KeyError:
            genes[specie] = [identifier]
    # Return the corresponding gene or truncate to 5 genes
    if genome:
        try:
            return genes[genome]
        except KeyError:
            raise FalseAbbreviation(string=genome, reaction=reaction)
    if len(genes) > 5:
        warn(
            message=(
                f'High amount of genes for reaction "{reaction}". '
                + "Species truncated to 5"
            ),
            category=UserWarning,
        )
        genes = dict(islice(genes.items(), 5))
    return list(chain.from_iterable(genes.values()))


def _get_genes(directory: Path, identifier: str):
    """
    Stores the genes for given reaction in given directory. Function will call
    a HTTPError if nothing is found.
    """
    # GENES directory will depend from the subdatabase
    directory = directory.joinpath("KEGG").joinpath("GENES")
    if not directory.exists():
        directory.mkdir()
    # Retrieval of the Gene information
    filename = directory.joinpath(f"{identifier}_genes.txt")
    if not filename.exists():
        # Obtain list of KO
        # This URL will not necessarily raise exception
        string = "+".join(_get_ko(identifier=identifier))
        url_text = f"http://rest.kegg.jp/link/genes/{string}"
        try:
            r = get(url_text)
            r.raise_for_status()
            with open(file=filename, mode="w+") as f:
                f.write(r.text)
        except HTTPError:
            raise HTTPError(f'Gene information for "{identifier}" unavailable')


def _p_entry_genes(
    kegg_dict: dict, directory: Path, identifier: str, genome: str
) -> dict:
    """
    From given KEGG dictionary returns a dictionary with the key
    "genes" which include a dictionary with the identifier and name of the
    gene; and the key "rule" for the COBRApy representation of the
    gene-reaction rule
    """
    genes: Dict[str, str] = dict()
    rule = str()
    filename = (
        directory.joinpath("KEGG")
        .joinpath("GENES")
        .joinpath(f"{identifier}_genes.txt")
    )
    with suppress(FileNotFoundError):
        with open(file=filename, mode="r") as f:
            genes_list = _parse_ko(
                string=str(f.read()), reaction=identifier, genome=genome
            )
            for gene in genes_list:
                genes[gene] = ""
            rule = " or ".join(genes.keys())
    # FIXME: Temporal OR rule
    if genes:
        rule = " or ".join(genes.keys())
    return {"genes": genes, "rule": rule}


def _p_reaction(kegg_dict: dict, directory: Path, genome: str) -> dict:
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
        identifier = list(
            chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
        )[-2]
        return {
            "TYPE": object_type,
            "ENTRY": identifier,
            "NAME": kegg_dict["NAME"][0],
            "DATABASE": "KEGG",
            "EQUATION": _build_dict_meta(line=kegg_dict["EQUATION"][0]),
            "BOUNDS": _get_reversibility(line=kegg_dict["EQUATION"][0]),
            "TRANSPORT": BaseParser._check_transport(
                data_dict=_build_dict_meta(line=kegg_dict["EQUATION"][0])
            ),
            "XREF": _build_reference(data_dict=kegg_dict),
            "GENES": _p_entry_genes(
                kegg_dict=kegg_dict,
                identifier=identifier,
                directory=directory,
                genome=genome,
            ),
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


def get_graph(kegg_dict: dict) -> dict:
    """
    Returns dictionary with sequences for a graph, where the key the prior
    reaction is and the value, its succesor; and a set with all participant
    reaction (vertex)
    """
    # Obtain sequences as a set and regular reactions
    sequences = [item.split(sep=" ")[0] for item in kegg_dict["REACTION"]]
    whole_set = set(
        chain.from_iterable(line.split(sep=",") for line in sequences)
    )
    reactions = [line.split(sep=",") for line in sequences]
    graph = dict()  # type: ignore
    line: list
    for index_j, line in enumerate(iterable=reactions):
        reaction: str
        for index_i, reaction in enumerate(iterable=line):
            # If index not found, then ignore it
            with suppress(IndexError):
                # In KEGG the order of the reactions is defined by its order in
                # the sequence.
                parent = reaction
                children = reactions[index_j + 1]
                if len(children) == 1:
                    # Take as a string
                    children = children[0]
                else:
                    # All elements as a tuple
                    children = tuple(children)
                # Expand keys
                try:
                    # If at least 2
                    if isinstance(graph[parent], tuple):
                        graph[parent] += children
                        continue
                    # Else transform single string to tuple
                    graph[parent] = (graph[parent], children)
                except KeyError:
                    # Create new one
                    graph[parent] = children
    for reaction in whole_set:
        try:
            graph[reaction]
        except KeyError:
            graph[reaction] = None
    if not graph:
        raise WrongParserError("Given root does not belong to a Pathway")
    return graph


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
        graph = get_graph(kegg_dict=kegg_dict)
        return {
            "TYPE": "Pathway",
            "ENTRY": list(
                chain.from_iterable(
                    item.split() for item in kegg_dict["ENTRY"]
                )
            )[0],
            "NAME": kegg_dict["NAME"][0],
            "DATABASE": "KEGG",
            "PATHWAY": graph,
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
    def _parse(root: dict, directory: Path, genome: str = None) -> dict:
        """
        Parses raw text from KEGG database and transforms it in a dictionary
        that can be used later as commom data type for cobramod. It identifies
        the type of object depending of it information.

        Args:
            root (dict): dictionary with information for an object

        Returns:
            dict: data from KEGG
        """
        for method in (_p_enzyme, _p_compound, _p_reaction, _p_pathway):
            with suppress(WrongParserError):
                try:
                    kegg_dict = method(  # type: ignore
                        kegg_dict=root, directory=directory, genome=genome
                    )
                except TypeError:
                    # TODO: Find better solution
                    try:
                        kegg_dict = method(  # type: ignore
                            kegg_dict=root, directory=directory
                        )
                    except TypeError:
                        kegg_dict = method(kegg_dict=root)  # type: ignore
                return kegg_dict
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
        try:
            genome = kwargs["genome"]
        except KeyError:
            genome = None
        return KeggParser._parse(root=raw, directory=directory, genome=genome)

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
        except TypeError:
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
        dict: dictionary from KEGGs identifier with basic information
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
        # Try to obtain the data locally or get it from KEGG and store it
        try:
            return KeggParser._read_file(filename=filename)
        except FileNotFoundError:
            debug_log.debug(
                f'"{identifier}" not found in directory "{database}".'
            )
            # Retrieve from URL
            url_text = f"http://rest.kegg.jp/get/{identifier}/"
            debug_log.debug(f"Searching in {url_text}")
            r = get(url_text)
            try:
                r.raise_for_status()
                unformatted_data = r.text
                debug_log.info(
                    f'Object "{identifier}" found in database. Saving in '
                    f'directory "{database}".'
                )
                with open(file=filename, mode="w+") as f:
                    f.write(unformatted_data)
                # Search for the gene information and store it if found.
                with suppress(HTTPError):
                    _get_genes(directory=directory, identifier=identifier)
                return _create_dict(raw=unformatted_data)
            except HTTPError:
                msg = f'"{identifier}" not available in "{database}".'
                debug_log.error(msg)
                raise HTTPError(msg)
    else:
        msg = "Directory not found. Please create the given directory."
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
