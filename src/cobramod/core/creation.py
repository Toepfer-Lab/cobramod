"""Object creation

This module handles the creation of COBRApy's objects
:class:`cobra.core.metabolite.Metabolite` and
:class:`cobra.core.reaction.Reaction`. Dictionaries containing the data of a
given database are used.
"""

from __future__ import annotations

from contextlib import suppress
from pathlib import Path
from typing import Optional, Union

import cobra.core as cobra_core
import requests

import cobramod.error as cmod_error
import cobramod.retrieval as cmod_retrieval
import cobramod.utils as cmod_utils
from cobramod.core import crossreferences
from cobramod.debug import debug_log


def build_metabolite(
    identifier: str, formula: str, name: str, charge: float, compartment: str
) -> cobra_core.Metabolite:
    """
    Returns a basic :class:`cobra.Metabolite`.

    LOG level: Debug

    Args:
        identifier (str): Short name for Metabolite
        formula (str): Chemical formula
        name (name): Long name for Metabolite
        charge (float): Charge of the metabolite
        compartment (str): Location

    Returns:
        Metabolite
    """
    metabolite = cobra_core.Metabolite(
        id=identifier,
        formula=formula,
        name=name,
        charge=charge,
        compartment=compartment,
    )
    debug_log.debug(f'Manually curated metabolite "{identifier}" was created.')
    return metabolite


def metabolite_from_string(
    line_string: str, replacement: dict[str, str], model: cobra_core.Model
) -> cobra_core.Metabolite:
    """
    Creates a Metabolite object based on a string. Function will try to find
    the metabolite in the model or cross-reference and return it if possible.

    The string must follow the syntax:

    :code:`formatted identifier, name , compartment, chemical_formula,
    molecular_charge`

    Args:
        line_string (str): string with information
        replacement (dict[str, str]): Dictionary with either the new identifier
            and/or the identifier of an object inside the model
        model (Model): Model to search for cross-references

    Raises:
        WrongSyntax: if format is does not follow syntax

    Returns:
        Metabolite
    """
    line = [part.strip().rstrip() for part in line_string.split(",")]

    if len(line) != 5:
        raise cmod_error.WrongSyntax(
            "Given line is invalid. Format is:\nid, name, compartment, "
            f"chemical_formula, molecular_charge.\nVerify line:\n{line_string}"
        )

    identifier = line[0]

    # Try to find replacement
    new_identifier = find_replacements(
        identifier=identifier,
        obj_type="metabolites",
        replace_dict=replacement,
        model=model,
    )
    if isinstance(new_identifier, cobra_core.Metabolite):
        return new_identifier

    with suppress(KeyError):
        identifier = replacement[identifier]

    name = line[1]
    compartment = line[2]
    formula = line[3]
    charge = line[4]

    debug_log.debug(
        f'Metabolite "{identifier}" was identified as a manually. '
        + "curated metabolite."
    )
    return build_metabolite(
        identifier=identifier,
        name=name,
        compartment=compartment,
        charge=float(charge),
        formula=formula,
    )


def metabolite_from_data(
    data: cmod_retrieval.Data,
    compartment: str,
    model: cobra_core.Model = cobra_core.Model(),
) -> cobra_core.Metabolite:
    """
    Parses the metabolite information from a Data object and creates a
    Metabolite in given compartment. This function will try to figure out if
    the metabolite is already in the model

    Returns:
        Metabolite
    """
    # Try to obtain available translation
    xref = cmod_utils.find_intersection(
        dictlist=model.metabolites,
        query=data.attributes["xref"],
        revert=True,
    )

    if xref:
        data.identifier = xref

    # No need for replacement dictionary for metabolites
    metabolite = data.parse(model, compartment)
    if not isinstance(metabolite, cobra_core.Metabolite):
        raise AttributeError("Given object is not a Metabolite")

    debug_log.debug(f"Metabolite '{metabolite.id}' created from Data object")
    return metabolite


def convert_string_metabolite(
    line: str,
    model: cobra_core.Model,
    replacement: dict[str, str],
    directory: Path,
    database: Optional[str],
    model_id: Optional[str],
) -> cobra_core.Metabolite:
    """
    Transform a string into a Metabolite and returns it. This object can be
    either manually-curated or from a database. If the compound if found under
    a different name (cross-reference), this will be returned instead.

    Custom metabolite syntax:

    :code:`formatted_identifier, name, compartment, chemical_formula,
    molecular_charge`

    Metabolite from database:

    :code:`metabolite_identifier, compartment`

    Args:
        line (str): string with information of metabolite
        model (Model): Model to check for cross-references
        replacement (dict[str, str]): Dictionary with either the new identifier
            and/or the identifier of an object inside the model
        directory (Path): Path to directory where data is located
        database (str, optional): Name of database. Check
            :obj:`cobramod.retrieval.available_databases` for a list of names.
            The argument can be None ONLY for manually curated metabolites
        model_id (str, optional): Exclusive for BIGG. Name of model

    Returns:
        Metabolite
    """
    if line.count(",") > 1:
        new_metabolite = metabolite_from_string(
            line_string=line, replacement=replacement, model=model
        )
    else:
        # Prep string
        segment = (part.strip().rstrip() for part in line.split(sep=","))
        identifier = next(segment)
        compartment = next(segment)

        if replacement:
            old_metabolite = find_replacements(
                identifier=identifier,
                replace_dict=replacement,
                model=model,
                obj_type="metabolites",
            )
            if isinstance(old_metabolite, cobra_core.Metabolite):
                return old_metabolite

        identifier = replacement.get(identifier, identifier)

        # Retrieve from database
        data = cmod_retrieval.get_data(
            identifier, directory, database, model_id
        )
        new_metabolite = data.parse(model, compartment, replacement)

    if not isinstance(new_metabolite, cobra_core.Metabolite):
        raise TypeError("Given object is not a valid COBRApy Metabolite")
    return new_metabolite


def get_file_metabolites(
    model: cobra_core.Model,
    filename: Path,
    replacement: dict[str, str],
    directory: Path,
    database: Optional[str],
    model_id: Optional[str],
) -> list[cobra_core.Metabolite]:
    """
    Return a list with Metabolites from a file. If a metabolite is found in
    model under a different name (cross-reference), it will be included in the
    list instead.

    Custom metabolite syntax:

    :code:`formatted_identifier, name, compartment, chemical_formula,
    molecular_charge`

    Metabolite from database:

    :code:`metabolite_identifier, compartment`

    Args:
        model (Model): Model to get metabolites if available
        filename (Path): location of the file with metabolites
        replacement (dict[str, str]): Dictionary with either the new identifier
            and/or the identifier of object inside the model.
        directory (Path): Path to directory where data is located.
        database (str, optional): Name of database. Check
            :obj:`cobramod.retrieval.available_databases` for a list of names.
            The argument can be None ONLY for manually curated metabolites
        model_id (str, optional): Bigg-specific argument. Name of the model

    Raises:
        FileNotFoundError: If given file is not found

    Returns:
        list[Metabolite]
    """
    if not filename.exists():
        raise FileNotFoundError(
            f'Given file in "{str(filename)}" does not exists. '
            + "Please create the given file or verify that the path is correct"
        )
    # For each line, build and add metabolite. If a Metabolite is no properly
    # created, either raise an Error or use a default.
    with open(filename, "r") as f:
        new_metabolites: list[cobra_core.Metabolite] = list()
        for line in cmod_utils.read_lines(f=f):
            new_metabolites.append(
                convert_string_metabolite(
                    line, model, replacement, directory, database, model_id
                )
            )
    return new_metabolites


def get_reaction(
    data: cmod_retrieval.Data,
    compartment: str,
    replacement: dict[str, str],
    show_imbalance: bool,
    stop_imbalance: bool,
    model: cobra_core.Model,
) -> cobra_core.Reaction:
    """
    Creates a :class:`cobra:Reaction` from reaction Data. Location of the
    reaction can be set with the argument 'compartment'.

    This function searchs if the reaction and their metabolites are already in
    the model and will use them instead. Additionally, genes will be created
    and linked to the reaction.

    Args:
        data_dict (Data): Biochemical information of the reaction
        compartment (str): Locations of the reaction
        replacement (dict[str, str]): Original identifiers to be replaced
            Values are the new identifiers.
        model (Model): Model to search for duplicates and cross-references
        stop_imbalance (bool): If unbalanced reaction is found, stop process
        show_imbalance (bool): If unbalanced reaction is found, show output

    Returns:
        Reaction
    """
    if data.database == "META":
        msg = (
            f'Metabolic pathway information for reaction "{data.identifier}" '
            'comes from database "META". No genes will be created. Please use '
            "a specie-specific sub-database from BioCyc."
        )
        debug_log.warning(msg)

    xref = cmod_utils.find_intersection(
        dictlist=model.reactions,
        query=data.attributes["xref"],
        revert=True,
    )

    if xref:
        data.identifier = xref

    data.attributes["replacements"] = replacement
    reaction = data.parse(model, compartment, replacement)

    if not isinstance(reaction, cobra_core.Reaction):
        raise AttributeError("Given object is not a Metabolite")

    cmod_utils.check_imbalance(
        reaction=reaction,
        show_imbalance=show_imbalance,
        stop_imbalance=stop_imbalance,
    )

    return reaction


def reaction_from_string(
    line_string: str,
    directory: Path,
    database: Optional[str],
    stop_imbalance: bool,
    show_imbalance: bool,
    replacement: dict[str, str],
    model_id: Optional[str],
    model: cobra_core.Model = cobra_core.Model(),
) -> cobra_core.Reaction:
    """
    Returns a custom reaction from a string that includes the information
    of the reaction and its metabolites.

    If the metabolites are not in inside the model, they will be retrieved from
    the specified database. This Function will also search for cross-references
    in the model.

    Syntax:

    :code:`reaction_identifier, reaction_name | coefficient metabolite <->
    coefficient metabolite

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`4 OXYGEN-MOLECULE_c`**

    Else, include a model to retrieve identifiers from it.

    Reversibility will depend on the type of arrow. Options are:

    **`<--, -->, <=>, <->`**

    Args:
        line_string (str): string with information
        directory (Path): Path to directory where data is located
        database (str, optional): Name of database. Check
            :obj:`cobramod.retrieval.available_databases` for a list of names.
            The argument can be None ONLY for manually curated metabolites
        stop_imbalance (bool): If unbalanced reaction is found, stop process
        show_imbalance (bool): If unbalanced reaction is found, show output
        replacement (dict): Dictionary with either the new identifier and/or
            the identifier of an object inside the model
        model_id (optional, str): Bigg-specific argument. Name of the model
        model (Model, optional): A model to obtain metabolite objects from.
            Defaults to an empty Model

    Raises:
        NoDelimiter: If no delimiter "|" is found.

    Returns:
        Reaction
    """
    # Remove blank and obtains information
    try:
        info, reaction_str = [
            string.strip().rstrip() for string in line_string.split("|")
        ]
    except ValueError:
        raise cmod_error.NoDelimiter(line_string)

    # It is possible to pass only the identifier
    if "," in info:
        reaction_id, reaction_name = [
            string.strip().rstrip() for string in info.split(",")
        ]
    else:
        reaction_id = reaction_name = info

    # Create Base reaction and then fill it with its components.
    reaction = cobra_core.Reaction(id=reaction_id, name=reaction_name)

    cmod_retrieval.build_reaction_from_str(
        model,
        reaction,
        reaction_str,
        directory,
        database,
        model_id,
        replacement,
    )

    # Set new bounds and check imbalance
    # new_reaction.bounds = bounds
    cmod_utils.check_imbalance(
        reaction=reaction,
        show_imbalance=show_imbalance,
        stop_imbalance=stop_imbalance,
    )
    return reaction


def string_to_reaction(
    line: str,
    model: cobra_core.Model,
    directory: Path,
    database: Optional[str],
    stop_imbalance: bool,
    show_imbalance: bool,
    replacement: dict[str, str],
    model_id: Optional[str],
    genome: Optional[str],
) -> cobra_core.Reaction:
    """
    Returns a :class:`cobra.Reaction` from a string. It can be either
    manually-curated, or the identifier for a reaction in a database. The
    function will search for cross-references of the reaction and its
    metabolites if the database provides that information

    Syntax:

    :code:`reaction_identifier, compartment`

    For custom reactions:

    :code:`reaction_identifier, reaction_name | coefficient metabolite <->
    coefficient metabolite

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`4 OXYGEN-MOLECULE_c`**

    Else, include a model to retrieve identifiers from it.

    Reversibility will depend on the type of arrow that. Options are

    **`<--, -->, <=>, <->`**

    Args:
        line (str): String with custom reaction or identifier of reaction
        model (Model): Model to search for cross-references
        directory (Path): Path to directory, where data is located
        database (str, optional): Name of database. Check
            :obj:`cobramod.retrieval.available_databases` for a list of names.
            This argument can be empty ONLY for manually-curated reactions
        stop_imbalance (bool): If unbalanced reaction is found, stop process
        show_imbalance (bool): If unbalanced reaction is found, show output
        replacement (dict[str, str]): Dictionary with either the new identifier
            and/or the identifier of an object inside the model.
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model

    Returns:
        Reaction
    """
    if "|" in line:
        # Create custom reaction
        reaction = reaction_from_string(
            line_string=line,
            directory=directory,
            database=database,
            model=model,
            stop_imbalance=stop_imbalance,
            show_imbalance=show_imbalance,
            replacement=replacement,
            model_id=model_id,
        )

    else:
        # Must obtain from database
        identifier, compartment = (
            part.strip().rstrip() for part in line.split(",")
        )

        if database is None:
            raise TypeError(
                "Function 'string_to_reaction' cannot retrieve metabolic "
                "pathway information if no database if provided"
            )

        reaction = create_object(
            identifier,
            directory,
            database,
            compartment,
            replacement,
            show_imbalance,
            stop_imbalance,
            model,
            model_id,
            genome,
        )

        if not isinstance(reaction, cobra_core.Reaction):
            raise AttributeError("Given object is not a valid COBRApy Reaction")

        cmod_utils.check_imbalance(reaction, stop_imbalance, show_imbalance)

    return reaction


def get_file_reactions(
    model: cobra_core.Model,
    filename: Union[str, Path],
    directory: Path,
    stop_imbalance: bool,
    show_imbalance: bool,
    replacement: dict[str, str],
    database: Optional[str],
    model_id: Optional[str],
    genome: Optional[str],
) -> list[cobra_core.Reaction]:
    """
    Returns list with reactions from file. All reactions can be either created
    manually or retrieved from a database. For each reactions, its always
    checks for mass balance.

    Custom reactions:

    :code:`reaction_identifier, reaction_name | coefficient metabolite <->
    coefficient metabolite

    Identifiers of metabolites have to end with an underscore and a
    compartment:

    E.g **`4 OXYGEN-MOLECULE_c`**

    Else, include a model to retrieve identifiers from it. Reversibility will
    depend on the type of arrow that. Options are:

    **`<--, -->, <=>, <->`**

    From database:

    :code:`original_identifier, compartment`


    Args:
        model (Model): model to check for cross-references
        filename (Path, str): location of the file with reaction information
        directory (Path): Path to directory where data is located
        stop_imbalance (bool): If unbalanced reaction is found, stop process
        show_imbalance (bool): If unbalanced reaction is found, show output
        replacement (dict[str, str]): Dictionary with either the new identifier
            and/or the identifier of an object inside the model
        database (str, optional): Name of database. Check
            :obj:`cobramod.retrieval.available_databases` for a list of names.
            This argument can be empty ONLY for manually-curated reactions
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie

    Raises:
        FileNotFoundError: If file does not exists

    Returns:
        list[Reaction]
    """
    if isinstance(filename, str):
        filename = Path(filename).absolute()

    if not filename.exists():
        raise FileNotFoundError(
            f'File "{filename.name}" does not exist. '
            + "Check if the given path is correct."
        )
    with open(filename, "r") as f:
        new_reactions = list()
        for line in cmod_utils.read_lines(f=f):
            new_reactions.append(
                string_to_reaction(
                    line=line,
                    model=model,
                    directory=directory,
                    database=database,
                    stop_imbalance=stop_imbalance,
                    show_imbalance=show_imbalance,
                    replacement=replacement,
                    model_id=model_id,
                    genome=genome,
                )
            )
    return new_reactions


def create_object(
    identifier: str,
    directory: Union[Path, str],
    database: Optional[str],
    compartment: str,
    replacement: dict[str, str] = {},
    show_imbalance: bool = True,
    stop_imbalance: bool = False,
    model: cobra_core.Model = cobra_core.Model(),
    model_id: Optional[str] = None,
    genome: Optional[str] = None,
) -> cobra_core.Object:
    """
    Creates and returns a COBRApy object based on given identifier and
    database. Identifier names will be formatted.

    .. hint:: Hyphens will become underscores.

    Args:
        identifier (str): Original identifier of the database.
        directory (Path): Path to directory where data is stored.
        database (Optional[str]): Name of database. Check
            :obj:`cobramod.retrieval.available_databases` for a list of names.
        compartment (str): Location of the object. In case of reaction, all
            metabolites will be included in the same location.
        replacement (dict[str, str]): Dictionary with either the new identifier
            and/or the identifier of an object inside the model.

    Arguments for reactions:
        stop_imbalance (bool): If an unbalanced reaction is found, stop the
            process. Defaults to False.
        show_imbalance (bool): If an unbalanced reaction is found, print
            output. Defaults to True.
        model (Model): Model in which cross-references for reaction
            and metabolites might be found and used. If nothing is specified,
            it uses an empty model by default

    Special arguments for databases:
        model_id (str, optional): Exclusive for BIGG. Retrieve object from
            specified model. Internal functions use the default: "universal"
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            species involved. Genes will be obtained for this species.
            List available at https://www.genome.jp/kegg/catalog/org_list.html

    Returns:
        cobra.core.Object: A Reaction or Metabolite object
    """
    query = replacement.get(identifier, identifier)

    if not isinstance(directory, Path):
        directory = Path(directory).absolute()

    try:
        data = cmod_retrieval.get_data(
            directory=directory,
            identifier=query,
            database=database,
            # debug_level=10,
            model_id=model_id,
            genome=genome,
        )
    except requests.HTTPError:
        # If not found, then rename
        data = cmod_retrieval.get_data(
            directory=directory,
            identifier=identifier,
            database=database,
            # debug_level=10,
            model_id=model_id,
            genome=genome,
        )
        data.identifier = replacement.get(identifier, identifier)

    if data.mode == "Reaction":
        obj = get_reaction(
            data,
            compartment,
            replacement,
            show_imbalance,
            stop_imbalance,
            model,
        )
    else:
        obj = data.parse(model, compartment, replacement)

    return obj


def add_metabolites(
    model: cobra_core.Model,
    obj: Union[Path, str, list[cobra_core.Object], list[str]],
    directory: Union[Path, str] = Path.cwd().joinpath("data"),
    database: Optional[str] = None,
    **kwargs,
):
    """
    Adds given object into the model. The options are:

     - Path or str: A file with the metabolites. E. g:
        Path.cwd().joinpath("file_with_names.txt") or "./file_with_names.txt"
     - List[Metabolites]: A list with multiple Metabolite objects.

        Custom metabolite syntax:

        :code:`formatted_identifier, name, compartment, chemical_formula,
        molecular_charge`

        Metabolite from database:

        :code:`metabolite_identifier, compartment`

     - List[str]: A list with multiple str with the mentioned syntax.

     CobraMod will try to download the biochemical information and create
     instead custom metabolites if not found in the database. In case of
     only custom metabolites, it is recommended to use None for the argument
     database

    Args:
        model (Model): Model to be expanded and searched for metabolites.
        obj: A Path; a list with either strings or Metabolite objects,
            or a single string. See syntax above.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names. Defaults
            to None (This is useful for custom metabolites).
        directory (Path): Path to directory where the data is located. If
            nothing is specified, then it defaults to the directory "data" in
            the current working directory.

    Keyword Arguments:
        replacement (dict): Dictionary with either the new identifier and/or
            the identifier of an object inside the model.
        include_metanetx_specific_ec: Determines whether MetaNetX specific
            EC numbers should be taken over. These are generally not found in
            other databases. Furthermore, this could result in non-existing
            Brenda IDs being created. The default value is False.

    Raises:
        WrongSyntax (from str): If the syntax is not followed correctly as
            mentioned above.
        ValueError: If Keyword Arguments are missing.
        FileNotFoundError (from Path): if file does not exists
    """
    # Defaults
    include_metanetx_specific_ec = kwargs.pop(
        "include_metanetx_specific_ec", False
    )
    if not isinstance(directory, Path):
        directory = Path(directory).absolute()

    # Kwargs
    replacement: dict[str, str] = kwargs.pop("replacement", {})
    model_id: str = kwargs.pop("model_id", {})

    if isinstance(obj, str):
        obj = Path(obj).absolute()

    metabolites: list[cobra_core.Metabolite]
    if isinstance(obj, Path):
        metabolites = get_file_metabolites(
            model, obj, replacement, directory, database, model_id
        )

    elif isinstance(obj, list):
        metabolites = []

        for item in obj:
            if isinstance(item, cobra_core.Metabolite):
                metabolites.append(item)

            elif isinstance(item, str):
                metabolites.append(
                    convert_string_metabolite(
                        item, model, replacement, directory, database, model_id
                    )
                )

    else:
        raise cmod_error.WrongDataError(
            "Given object is not the type mentioned in the docstrings. "
            + "Please use a string, a list or a pathlib 'Path' object."
        )

    for metabolite in metabolites:
        crossreferences.add_crossreferences(
            object=metabolite,
            directory=directory,
            include_metanetx_specific_ec=include_metanetx_specific_ec,
        )

    cmod_utils.confirm_metabolite(model, metabolites)


def find_replacements(
    identifier: str,
    obj_type: str,
    replace_dict: dict[str, str],
    model: cobra_core.Model,
) -> Optional[cobra_core.Object]:
    """
    Returns either a COBRApy object if found in the model or a string with the
    new identifier, which can be used to create a new object or retrieve data.
    The replace_dict has two options for keys. One can be the unmodified
    identifier given by a database, or a formated identifier. E.g:
    "WATER_c": "C00001_c" or "ACETALD-DEHYDROG-RXN" : "R00228_c"

    Args:
        identifier (str): Original identifier to replace or rename.
        obj_type (str): Type of object to search for. Only options are
            "reactions" or "metabolites".
        replace_dict (dict): Dictionary with either the new identifier and or
            the identifier of object inside the model.
        model (Model):

    Returns:
        A COBRApy Object if replacement is found in the model. Otherwise
        nothing is returned
    """
    cobra_dictlist: dict[str, cobra_core.DictList] = {
        "reactions": model.reactions,
        "metabolites": model.metabolites,
    }
    # Check if KeyErrors might appear, otherwise continue
    try:
        new_identifier = replace_dict[identifier]
        # Return either COBRApy object or new string to search up in database
        replacement = cobra_dictlist[obj_type].get_by_id(new_identifier)
        if not isinstance(replacement, cobra_core.Object):
            raise TypeError(
                f"The object for {identifier} is not a valid COBRApy object!"
            )

        debug_log.debug(
            f"{obj_type[:-1]} '{identifier}' was found in the model and will "
            " replace {new_identifier}."
        )
        return replacement
    except KeyError:
        return None


def add_reactions(
    model: cobra_core.Model,
    obj: Union[str, Path, list[cobra_core.Object], list[str]],
    directory: Union[Path, str],
    database: Optional[str] = None,
    stop_imbalance: bool = False,
    show_imbalance: bool = True,
    **kwargs,
):
    """Adds given object into the model. The options are:

     - Path: A file with components. E. g. :
        Path.cwd().joinpath("file_with_names.txt")
     - List[Reactions]: A list with regular Reactions
     - List[str]: Either the identifier with its corresponding compartment or a
     string with all components. This applies for the Path option. E.g. :

        :code:`reaction_identifier, compartment`

        For custom reactions

        :code:`reaction_identifier, reaction_name | coefficient metabolite <->
        coefficient metabolite

        Identifiers of metabolites have to end with an underscore and a
        compartment:

            E.g. **`4 OXYGEN-MOLECULE_c`**

     - List[str]: A list with multiple str with the mentioned syntax.

     CobraMod will try to use the same metabolites and reactions inside
     of the model. If not found, then CobraMod will try to download its
     biochemical information or create custom objects.

    Args:
        model (Model): Model to expand and search for reactions.
        obj: A Path; a list with either strings or Reaction objects,
            or a single string. See syntax above.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names. Defaults
            to None (Useful for custom reactions).
        directory (Path): Path to directory where data is located.
        stop_imbalance (bool): If an unbalanced reaction is found,
            stop the process. Defaults to False.
        show_imbalance (bool): If an unbalanced reaction is found, print
            output. Defaults to True.

    Keyword Arguments:
        replacement (dict): Dictionary with either the new identifier and/or
            the identifier of an object inside the model.
        genome (str): Exclusive for KEGG. Abbreviation for the species
            involved. Genes will be obtained for this species.
        consider_sub_elements: Specifies whether additional cross-references
            should also be added to the subelements. For example, you can
            specify whether only the reaction or also its metabolites
            should be expanded. Defaults to True
        include_metanetx_specific_ec: Determines whether MetaNetX specific
            EC numbers should be taken over. These are generally not found in
            other databases. Furthermore, this could result in non-existing
            Brenda IDs being created. The default value is False.

    Raises:
        WrongSyntax (from str): If the syntax is not followed correctly as
            mentioned above.
        ValueError: If Keyword Arguments are missing.
        FileNotFoundError (from Path): If the file does not exists.

    """
    # Defaults
    replacement: dict[str, str] = kwargs.pop("replacement", dict())
    genome: str = kwargs.pop("genome", None)
    model_id: str = kwargs.pop("model_id", None)
    consider_sub_elements = kwargs.pop("consider_sub_elements", True)
    include_metanetx_specific_ec = kwargs.pop(
        "include_metanetx_specific_ec", False
    )

    if isinstance(directory, str):
        directory = Path(directory).absolute()

    if isinstance(obj, str):
        obj = Path(obj).absolute()

    # In case of a Path
    reactions: list[cobra_core.Reaction]
    if isinstance(obj, Path):
        reactions = get_file_reactions(
            model=model,
            filename=obj,
            directory=directory,
            database=database,
            replacement=replacement,
            stop_imbalance=stop_imbalance,
            show_imbalance=show_imbalance,
            genome=genome,
            model_id=model_id,
        )

    elif isinstance(obj, list):
        reactions = []

        for item in obj:
            if isinstance(item, str):
                reactions.append(
                    string_to_reaction(
                        line=item,
                        model=model,
                        directory=directory,
                        database=database,
                        replacement=replacement,
                        stop_imbalance=stop_imbalance,
                        show_imbalance=show_imbalance,
                        genome=genome,
                        model_id=model_id,
                    )
                )
            if isinstance(item, cobra_core.Reaction):
                reactions.append(item)

    # Raise error if wrong
    else:
        raise cmod_error.WrongDataError(
            "Given object is not the type mentioned in the docstrings. "
            + "Please use a string, a list or a pathlib 'Path' object."
        )

    for reaction in reactions:
        # add cross references only if reactions is non custom
        if not database:
            continue

        crossreferences.add_crossreferences(
            object=reaction,
            directory=directory,
            consider_sub_elements=consider_sub_elements,
            include_metanetx_specific_ec=include_metanetx_specific_ec,
        )

    cmod_utils.add_reactions_to_model(model=model, reactions=reactions)
