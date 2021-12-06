#!/usr/bin/env python3
"""Data parsing for KEGG

This module handles the retrieval of data from KEGG into a local directory.
The possible type of data that can be download:

- Metabolite: Identifiers that start with the letter C, e.g C00001
- Reactions: Identifiers that start with the letter R, e.g R00001. The
gene information for reactions is also included if the specie is specified
- Module Pathways: Identifiers that start with the letter M, e.g M00001

Contact maintainers if other types should be added.

Important class of the module:
- KeggParser: Child of the abstract class
:class:`cobramod.parsing.base.BaseParser`.
"""
import io
from contextlib import suppress
from itertools import chain
from pathlib import Path
from typing import Generator, Union, Dict
from warnings import warn

from requests import get, HTTPError

from cobramod.debug import debug_log
from cobramod.parsing.base import BaseParser
from cobramod.error import WrongParserError, AbbreviationWarning


def _kegg_info_to_version(info: str) -> str:
    with io.StringIO(info) as lines:
        for num, line in enumerate(lines, 1):
            if num == 2:
                return line[line.find("Release") + 8 :].rstrip()

    msg = (
        "Error determining the kegg version. "
        'Instead, "Undefined" is used as version.'
    )
    warn(
        message=msg,
        category=UserWarning,
    )
    return "Undefined"


def _get_keys(raw: str) -> Generator:
    """
    Returns as Generator the keys for KEGG data. Input is raw text
    """
    lines = (line for line in raw.split(sep="\n"))
    for line in lines:
        segment = line.split(" ")
        if segment[0] != "///" and segment[0]:
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
    Formats most of the keys for KEGG data and returns a dictionary.
    """
    lines = (line for line in raw.split("\n"))
    keys = set(_get_keys(raw=raw))
    actual_key = str()
    kegg_dict: Dict[str, list] = {}
    # Create dictionary
    for line in lines:
        if "///" in line:
            break
        key = _find_key(line=line, keys=keys)
        if not key:
            # Append to older key.
            kegg_dict[actual_key].append(line.strip().rstrip())
        else:
            actual_key = key
            line = line[line.find(key) + len(key) :].strip().rstrip()
            kegg_dict[actual_key] = [line]
    return kegg_dict


def _get_metabolites(line: str) -> tuple:
    """
    Returns the bounds of the reaction as a tuple and a dictionary with the
    participant metabolites of the reaction

    Args:
        line (str): Equation in string form

    Returns:
        Tuple: bounds (tuple), dictionary with metabolites

    Raises:
        UnboundLocalError: If reversibility cannot be found for given reaction
    """
    # Obtaining bounds
    arrows = {"<=>": (-1000, 1000), "==>": (0, 1000), "<==": (-1000, 1000)}
    for item in arrows.keys():
        if item in line:
            position = line.find(item)
            arrow = line[position : position + 3]
    try:
        arrow
    except UnboundLocalError:
        raise UnboundLocalError(
            f'Reversibility for "{line}" could not be found.'
        )
    # Getting metabolites
    metabolites = dict()
    reactants, products = [
        side.rstrip().strip().split("+") for side in line.split(arrow)
    ]
    for side in (reactants, products):
        FACTOR = 1
        prefix = "r"
        if side is reactants:
            FACTOR = -1
            prefix = "l"
        for metabolite in side:
            try:
                coefficient, identifier = [
                    item.strip().rstrip()
                    for item in metabolite.rstrip().strip().split(" ")
                ]
            except ValueError:
                # In case no coefficient is given, it must be 1
                identifier = metabolite.strip().rstrip()
                coefficient = "1"
            metabolites[f"{prefix}_{identifier}"] = float(coefficient) * FACTOR
    return arrows[arrow], metabolites


def _build_reference(data_dict: dict) -> Union[dict, None]:
    """
    Return a dictionary, where the keys are the names of cross-references
    and keys their identifiers. If nothing is found, it will return None
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
    Returns a list with the corresponding KO-entries for given identifier.
    Otherwise it will raise a HTTPError
    """
    url_text = f"http://rest.kegg.jp/link/ko/{identifier}"
    try:
        response = get(url_text)
        response.raise_for_status()
        return [
            single
            for single in response.text.split()
            # This statement removes identifier from the split. e.g
            # [rn:R00001, ko:k00001]
            if identifier not in single
        ]
    except HTTPError:
        raise HTTPError(
            f'No K-Identifiers could be retrieved for "{identifier}".'
        )


def _parse_ko(string: str, reaction: str, genome: str = None) -> list:
    """
    Returns with a list with the corresponding genes for given genome. String
    is the raw text that include the gene information. If Abbreviation or no
    genome is present then, an empty list is returned.
    """
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
    # Return the corresponding gene or empty list
    if genome:
        try:
            return genes[genome]
        except KeyError:
            msg = (
                f'Reaction "{reaction}" does not have a "{genome}" '
                + "abbreviation as a specie. No genes will be added."
            )
            warn(message=msg, category=AbbreviationWarning)
            debug_log.warning(msg=msg)
            return []
    msg = (
        f'Nothing was specified in argument "genome". Reaction "{reaction}"'
        " will not include genes. Please modify if necessary."
    )
    warn(message=msg, category=UserWarning)
    debug_log.warning(msg=msg)
    return []


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
            response = get(url_text)
            response.raise_for_status()
            with open(file=filename, mode="w+") as file:
                file.write(response.text)
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
        with open(file=filename, mode="r") as file:
            genes_list = _parse_ko(
                string=str(file.read()), reaction=identifier, genome=genome
            )
            for gene in genes_list:
                genes[gene] = ""
            rule = " or ".join(genes.keys())
    if genes:
        rule = " or ".join(genes.keys())
        debug_log.warning(
            f'Gene-reaction rule for reaction "{identifier}" from KEGG '
            + ' set to "OR".'
        )
    return {"genes": genes, "rule": rule}


def _p_reaction(kegg_dict: dict, directory: Path, genome: str) -> dict:
    """
    Returns a dictionary that is the representation of a reaction
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
        # Obtain bounds
        bounds, metabolites = _get_metabolites(line=kegg_dict["EQUATION"][0])
        with suppress(KeyError):
            if "irreversible" in ",".join(kegg_dict["COMMENT"]):
                bounds = (0, 1000)
        return {
            "TYPE": object_type,
            "ENTRY": identifier,
            "NAME": kegg_dict["NAME"][0],
            "DATABASE": "KEGG",
            "EQUATION": metabolites,
            "BOUNDS": bounds,
            "TRANSPORT": BaseParser._check_transport(data_dict=metabolites),
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
    Returns a dictionary that is the representation of a metabolite
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
            # When formulas includes 'n' e.g 'H4P2O7(HPO3)n'
            if "(" in kegg_dict["FORMULA"][0]:
                formula = "X"
                msg = (
                    f'Sum formula for metabolite "{identifier}" from KEGG '
                    'could not be found. Formula set to "X" and charge to 0. '
                    "Please modify it if necessary."
                )
                debug_log.warning(msg=msg)
                warn(msg)
            else:
                formula = "".join(kegg_dict["FORMULA"]).strip().rstrip()
        except KeyError:
            # Some general compounds do not come with formulas
            formula = "X"
            msg = (
                f'Sum formula for metabolite "{identifier}" from KEGG '
                'could not be found. Formula set to "X" and charge to 0. '
                "Please modify it if necessary."
            )
            debug_log.warning(msg=msg)
            warn(msg)
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
    reaction is and the value, its successor; and a set with all participant
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
    # Check that the graph includes all reactions
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
    Returns a dictionary, which is the representation of a pathway.
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
            raise NotImplementedError(
                "Enzymes are currently not implemented. Please contact "
                + "maintainers."
            )
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
            "Please contact maintainers."
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
        KeggParser._check_database(directory=directory, database=database)
        raw = retrieve_data(directory=directory, identifier=identifier)
        debug_log.log(
            level=debug_level,
            msg=f'Data for identifier "{identifier}" '
            + "retrieved from KEGG.",
        )
        try:
            genome = kwargs["genome"]
        except KeyError:
            genome = None
        return KeggParser._parse(root=raw, directory=directory, genome=genome)

    @staticmethod
    def _check_database(directory: Path, database: str):
        """
        Returns the name of the database. This method is used to compare with
        given database name. It will raise a warning if both names are not
        equal or belong to the list of proper names.
        """
        if database != "KEGG":
            raise WrongParserError

    @staticmethod
    def _read_file(filename: Path) -> dict:
        """
        Reads the given file a returns a KEGG dictionary with most important
        information from it.
        """
        try:
            with open(file=filename, mode="r") as file:
                unformatted_data = file.read()
            return _create_dict(raw=unformatted_data)
        except TypeError:
            raise WrongParserError("Wrong filetype")


def retrieve_data(directory: Path, identifier: str) -> dict:
    """
    Retrieves and stores the data of given identifier for the KEGG database.

    Args:
        directory (Path): Path for data storage.
        identifier (str): official name for given object in KEGG

    Raises:
        NotImplementedError: If identifier contains a semicolon (:)
        HTTPError: If given identifier is not found in KEGG.
        NotADirectoryError:  If parent directory is not found.

    Returns:
        dict: dictionary from KEGG identifier with basic information
    """
    # NOTE: As KEGG only support xml for pathways, Metabolites and Reactions
    # have to be parsed differently. Only specific keys are necessary
    if ":" in identifier:
        # Get rid of prefix
        identifier = identifier[identifier.find(":") + 1 :]
    if directory.exists():
        database = "KEGG"
        data_dir = directory.joinpath(database)
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
            debug_log.debug(f"Searching {url_text} for biochemical data.")
            r = get(url_text)
            try:
                r.raise_for_status()
                unformatted_data = r.text

                database_version = get("http://rest.kegg.jp/info/kegg")
                database_version.raise_for_status()
                database_version = database_version.text
                database_version = _kegg_info_to_version(database_version)

                BaseParser.check_database_version(
                    directory, "kegg", database_version
                )

                debug_log.info(
                    f'Object "{identifier}" found in database. Saving in '
                    f'directory "{database}".'
                )
                data_dir.mkdir(exist_ok=True)
                with open(file=filename, mode="w+") as file:
                    file.write(unformatted_data)
                # Search for the gene information and store it if found.
                with suppress(HTTPError):
                    _get_genes(directory=directory, identifier=identifier)
                return _create_dict(raw=unformatted_data)
            except HTTPError:
                msg = (
                    f'"{identifier}" not available in "{database}". '
                    + "Please verify the name of the identifier."
                )
                debug_log.critical(msg)
                raise HTTPError(msg)
    else:
        msg = "Directory not found. Please create the given directory."
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
