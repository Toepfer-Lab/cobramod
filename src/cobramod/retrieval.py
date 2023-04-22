"""
"""
from __future__ import annotations

import json
import xml.etree.ElementTree as et
from pathlib import Path
from typing import Any, Literal, Optional, Union
from warnings import warn

import cobra.core as cobra_core
import requests

import cobramod.parsing as cmod_parsing
import cobramod.core as cmod_core
import cobramod.utils as cmod_utils
from cobramod.debug import debug_log

BIOCYC = "https://websvc.biocyc.org/getxml?id="
PMN = "https://pmn.plantcyc.org/getxml?"
KEGG = "https://rest.kegg.jp/get/"
BIGG = "https://bigg.ucsd.edu/api/v2/models/universal/"

EXTENSIONS = [".xml", ".txt", ".json"]

databases = {"biocyc": BIOCYC, "plantcyc": PMN, "kegg": KEGG, "bigg": BIGG}


class Data:
    identifier: str
    mode: str
    database: str
    attributes: dict[str, Any]
    path: Path

    def __init__(
        self,
        identifier: str,
        attributes: dict[str, str],
        mode: Literal["Pathway", "Reaction", "Metabolite"],
        database: str,
        path: Path,
    ):
        self.mode = mode
        self.identifier = identifier
        self.attributes = attributes
        self.database = database
        self.path = path

    def __repr__(self):
        return f"Data [{self.mode}: {self.identifier}]"

    @classmethod
    def from_json(
        cls, entry: str, data: dict[str, Any], database: str, path: Path
    ):
        is_compound = data.get("formulae", None)

        if is_compound:
            mode = "Metabolite"
            attributes = cmod_parsing.bigg.parse_metabolite_attributes(data)
        else:
            mode = "Reaction"
            attributes = cmod_parsing.bigg.parse_reaction_attributes(
                data, entry
            )

        return cls(entry, attributes, mode, database, path)

    @classmethod
    def from_text(cls, entry: str, text: str, database: str, path: Path):
        data: dict[str, list[str]] = cmod_parsing.kegg.data_from_string(text)
        mode = data["ENTRY"][0].split()[-1]

        if mode == "Compound":
            mode = "Metabolite"
            attributes = cmod_parsing.kegg.parse_metabolite_attributes(
                data, entry
            )

        elif mode == "Reaction":
            attributes = cmod_parsing.kegg.parse_reaction_attributes(
                data, entry
            )
        else:
            raise AttributeError(f"Cannot parse type '{mode}'")

        return cls(entry, attributes, mode, database, path)

    @classmethod
    def from_xml(cls, entry: str, root: et.Element, database: str, path: Path):
        if root.findall("Compound") or root.findall("Protein"):
            mode = "Metabolite"
            attributes = cmod_parsing.biocyc.parse_metabolite_attributes(
                root, entry
            )

        elif root.findall("Reaction"):
            mode = "Reaction"
            attributes = cmod_parsing.biocyc.parse_reaction_attributes(
                root, entry
            )
        else:
            raise AttributeError(
                "Cannot infere the object type from given et.Element"
            )
        return cls(entry, attributes, mode, database, path)

    def parse(
        self, model: cobra_core.Model, compartment: str
    ) -> cobra_core.Object:
        identifier = f"{self.identifier.replace('-','_')}_{compartment}"

        if self.mode == "Metabolite":
            try:
                return model.metabolites.get_by_id(identifier)

            except KeyError:
                return cmod_utils.build_metabolite(
                    identifier,
                    self.attributes["formula"],
                    self.attributes["name"],
                    self.attributes["charge"],
                    compartment,
                )
        elif self.mode == "Reaction":
            replacement = cmod_utils.find_intersection(
                dictlist=model.reactions,
                query=self.attributes["xref"],
                revert=True,
            )
            if replacement:
                self.identifier = replacement

            try:
                return model.reactions.get_by_id(identifier)

            except KeyError:
                reaction = cobra_core.Reaction(
                    identifier, self.attributes["name"]
                )
                model.add_reactions([reaction])

                reaction_str = ""
                part: str
                for part in self.attributes["equation"].split(" "):
                    # Removal of prefixes
                    if cmod_utils.is_compound(part):
                        part = part[2:]
                        part += f"_{compartment}"

                    part = self.attributes["replacements"].get(part, part)

                    reaction_str += f"{part} "

                if self.attributes["transport"]:
                    reaction_str = cmod_utils.convert_to_transport(
                        reaction_str, compartment
                    )
                    msg = (
                        f'Reaction "{reaction.id}" has the same metabolite on'
                        " both sides of the equation (e.g transport reaction). "
                        "COBRApy ignores these metabolites. To avoid this, by "
                        "default, CobraMod will assign the reactant side to "
                        "the extracellular compartment. Please curate the "
                        "reaction if necessary."
                    )
                    debug_log.warning(msg=msg)
                    warn(message=msg, category=UserWarning)

                directory = self.path.parents[1]

                if self.database == "BIGG":
                    directory = self.path.parents[2]

                build_reaction_from_str(
                    model,
                    reaction,
                    reaction_str,
                    directory,
                    self.database,
                    # NOTE: add model_id
                    "",
                )
                return reaction

        else:
            raise AttributeError("Cannot parse given data. Contact maintainers")


def get_response(
    query: str, model_id: Optional[str] = None
) -> tuple[str, requests.Response]:
    """
    Request different servers to see if there a valid response and returns the
    database name with the response

    model_id is a BIGG-specific argument.

    The content type must be xml, json or text
    """
    if not model_id:
        model_id = "universal"

    for database, web in databases.items():
        url = web + query

        try:
            if database == "bigg":
                response, _ = cmod_parsing.bigg.find_url(model_id, query)
                # This extra attribute if later use to save the file in their
                # corresponding directories
                response.extra = Path(model_id)

            else:
                response = requests.get(url)
                response.raise_for_status()

            # Not a valid content-type to process
            header = response.headers.get("Content-Type", "").lower()

            if "html" in header:
                continue

            return database, response

        except (requests.HTTPError, requests.Timeout, requests.ConnectionError):
            continue

    raise AttributeError(
        f"No website was able to return information from the query: {query}"
    )


def file_to_Data_class(identifier: str, filename: Path) -> Data:
    f = open(filename, "r")
    text = f.read()
    f.close()

    suffix = filename.suffix
    parent = filename.parent.name
    bigg_parent = filename.parents[1].name

    if suffix == ".json":
        return Data.from_json(
            identifier, json.loads(text), bigg_parent, filename
        )

    elif suffix == ".xml":
        return Data.from_xml(identifier, et.fromstring(text), parent, filename)

    elif suffix == ".txt":
        return Data.from_text(identifier, text, parent, filename)

    else:
        raise AttributeError("Cannot parse given content type")


def write(name: str, directory: Path, response: requests.Response) -> Path:
    """
    Finds the content type and writes the data into disk and returns the Path
    of the new file
    """

    header = response.headers.get("Content-Type", "").lower()

    if "json" in header:
        prefix = ".json"

    elif "xml" in header:
        prefix = ".xml"

    elif "text" in header:
        prefix = ".txt"

    else:
        raise AttributeError("Cannot parse given content type")

    filename = directory.joinpath(name + prefix)

    with open(file=filename, mode="w+") as file:
        file.write(response.text)

    return filename


def get_files(directory: Path, identifier: str):
    for item in directory.iterdir():
        for pattern in EXTENSIONS:
            if not item.match(identifier + pattern):
                continue
            yield item


def get_data(
    identifier: str,
    directory: Union[str, Path],
    database: Optional[str],
    model_id: Optional[str] = None,
    # NOTE: add docstring
):
    if isinstance(directory, str):
        directory = Path(directory).absolute()

    # Try first locally and the query databases
    if not database:
        try:
            filename = next(directory.rglob(identifier + "*"))

        except StopIteration:
            database, response = get_response(identifier)

            if database == "bigg":
                extra: Path = getattr(response, "extra")

                if not extra:
                    raise AttributeError(
                        "The 'extra' attribute was not found in the response"
                    )
                directory = directory.joinpath(extra)
            filename = write(identifier, directory, response)
    else:
        try:
            if model_id:
                filename = next(
                    get_files(
                        directory.joinpath(database, model_id), identifier
                    )
                )
            else:
                filename = next(
                    get_files(directory.joinpath(database), identifier)
                )

        except StopIteration:
            database, response = get_response(identifier, model_id)

            if database == "bigg":
                extra: Path = getattr(response, "extra")

                if not extra:
                    raise AttributeError(
                        "The 'extra' attribute was not found in the response"
                    )
                directory = directory.joinpath("BIGG", extra)

            filename = write(identifier, directory, response)

    return file_to_Data_class(identifier, filename)


def build_reaction_from_str(
    model: cobra_core.Model,
    reaction: cobra_core.Reaction,
    reaction_str: str,
    directory: Path,
    database: str,
    model_id: str,
):
    position = cmod_utils.get_arrow_position(reaction_str)
    arrow = reaction_str[position : position + 3]
    reaction.bounds = cmod_utils.ARROWS[arrow]

    compounds: dict[str, float] = {}
    # Reactants
    for pair in reaction_str[:position].rstrip().strip().split("+"):
        try:
            coefficient, metabolite = pair.rstrip().strip().split(" ")
        except ValueError:
            coefficient = 1
            metabolite = pair.rstrip().strip()

        compounds[metabolite] = -1.0 * float(coefficient)

    # products
    for pair in reaction_str[position + 3 :].rstrip().strip().split("+"):
        try:
            coefficient, metabolite = pair.rstrip().strip().split(" ")
        except ValueError:
            coefficient = 1
            metabolite = pair.rstrip().strip()

        compounds[metabolite] = float(coefficient)

    for identifier, coefficient in compounds.items():
        compartment = identifier[-1]
        try:
            data = get_data(identifier[:-2], directory, database, model_id)
            metabolite = cmod_core.creation.get_metabolite(
                data, compartment, model
            )
        except AttributeError:
            # NOTE: add log. This is only for custom
            metabolite = cobra_core.Metabolite(
                identifier, name=identifier, compartment=compartment
            )

        # NOTE: add logs
        if not isinstance(metabolite, cobra_core.Metabolite):
            raise TypeError("Given object is not a valid COBRApy Metabolite")

        reaction.add_metabolites({metabolite: coefficient})
