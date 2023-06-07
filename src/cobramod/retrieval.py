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

import cobramod.utils as cmod_utils
from cobramod.core import creation, genes
from cobramod.debug import debug_log
from cobramod.parsing import bigg, biocyc, kegg, plantcyc, solcyc

BIOCYC = "https://websvc.biocyc.org/getxml?id="
PMN = "https://pmn.plantcyc.org/getxml?id="
KEGG = "https://rest.kegg.jp/get/"
SOLCYC = "https://solcyc.sgn.cornell.edu/getxml?id="

EXTENSIONS = [".xml", ".txt", ".json"]

# NOTE: bigg uses an extra function
databases = {
    "BIOCYC": BIOCYC,
    "PMN": PMN,
    "KEGG": KEGG,
    "BIGG": "",
    "SOL": SOLCYC,
}


class Data:
    identifier: str
    mode: str
    database: str
    attributes: dict[str, Any]
    path: Path
    model_id: str
    genome: str

    def __init__(
        self,
        identifier: str,
        attributes: dict[str, str],
        mode: Literal["Pathway", "Reaction", "Metabolite"],
        database: str,
        path: Path,
        model_id: str,
    ):
        self.mode = mode
        self.identifier = identifier
        self.attributes = attributes
        self.database = database
        self.path = path
        self.model_id = model_id
        self.genome = ""

    def __repr__(self):
        return f"Data [{self.mode}: {self.identifier}]"

    @classmethod
    def from_json(
        cls, entry: str, data: dict[str, Any], database: str, path: Path
    ):
        is_compound = data.get("formulae", data.get("formula", None))
        model_id = path.parent.name

        # Non universal models includes suffices
        if model_id != "universal":
            entry = entry[: entry.find("_")]

        if is_compound:
            mode = "Metabolite"
            attributes = bigg.parse_metabolite_attributes(data)
        else:
            mode = "Reaction"
            attributes = bigg.parse_reaction_attributes(data, entry)

        return cls(entry, attributes, mode, database, path, model_id)

    @classmethod
    def from_text(
        cls, entry: str, text: str, database: str, path: Path, genome: str
    ):
        data: dict[str, list[str]] = kegg.data_from_string(text)
        mode = data["ENTRY"][0].split()[-1]
        gene_path = path.parent.joinpath("GENES")

        if mode == "Compound":
            mode = "Metabolite"
            attributes = kegg.parse_metabolite_attributes(data, entry)

        elif mode == "Reaction":
            attributes = kegg.parse_reaction_attributes(
                data, entry, genome, gene_path
            )

        elif mode == "Module":
            mode = "Pathway"
            attributes = kegg.parse_pathway_attributes(data, entry)
        else:
            raise AttributeError(f"Cannot parse type '{mode}'")

        obj = cls(entry, attributes, mode, database, path, "")
        obj.genome = genome
        return obj

    @classmethod
    def from_xml(cls, entry: str, root: et.Element, database: str, path: Path):
        gene_path = path.parent.joinpath("GENES")

        if root.findall("Compound") or root.findall("Protein"):
            mode = "Metabolite"
            attributes = biocyc.parse_metabolite_attributes(root, entry)

        elif root.findall("Reaction"):
            mode = "Reaction"
            attributes = biocyc.parse_reaction_attributes(
                root, entry, gene_path
            )

        elif root.findall("Pathway"):
            mode = "Pathway"
            attributes = biocyc.parse_pathway_attributes(root, entry)
        else:
            raise AttributeError(
                "Cannot infere the object type from given et.Element"
            )
        return cls(entry, attributes, mode, database, path, "")

    def parse(
        self,
        model: cobra_core.Model,
        compartment: str,
        replacement: dict[str, str] = {},
    ) -> cobra_core.Object:
        identifier = replacement.get(self.identifier, self.identifier)
        identifier = f"{identifier.replace('-','_')}_{compartment}"

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
            xref = cmod_utils.find_intersection(
                dictlist=model.reactions,
                query=self.attributes["xref"],
                revert=True,
            )
            if xref:
                self.identifier = xref

            try:
                reaction = model.reactions.get_by_id(identifier)
                if not isinstance(reaction, cobra_core.Reaction):
                    raise TypeError("Not a valid COBRApy Reaction")
                genes.genes_to_reaction(reaction, self.attributes["genes"])
                return reaction

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

                    part = replacement.get(part, part)

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
                    self.model_id,
                    replacement,
                )

                genes.genes_to_reaction(reaction, self.attributes["genes"])

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

    subdatabase = ""

    biocyc_credentials = False

    try:
        # e.g sol:META:water
        database, subdatabase, identifier = query.split(":")
        url = f"{databases[database]}{subdatabase}:{identifier}"

    except ValueError:
        database, identifier = query.split(":")

        if database == "KEGG":
            url = databases["KEGG"]
        elif database == "BIGG":
            url = databases["BIGG"]
        else:
            url = f"{databases['BIOCYC']}{database}:{identifier}"
            biocyc_credentials = True

    if database == "BIGG":
        response, _ = bigg.find_url(model_id, query)
        response.raise_for_status()

        # This extra attribute if later use to save the file in their
        # corresponding directories
        response.extra = Path(model_id)  # type: ignore

        return database, response

    s = requests.Session()
    if biocyc_credentials:
        # FIXME: find a better way
        user, pwd = cmod_utils.get_credentials(
            Path.cwd().joinpath("credentials.txt")
        )
        s.post(
            "https://websvc.biocyc.org/credentials/login/",
            data={"email": user, "password": pwd},
        )

    response = s.get(url)
    s.close()
    response.raise_for_status()

    # Not a valid content-type to process
    header = response.headers.get("Content-Type", "").lower()

    if response.text[:300].find("Subscription") != -1 and "html" in header:
        raise requests.HTTPError(
            f"Cannot retrieve data from '{url}'. It is possible that a "
            "subscription is required. Try using 'META'"
        )

    return database, response


def file_to_Data_class(
    identifier: str, filename: Path, genome: Optional[str]
) -> Data:
    f = open(filename, "r")
    text = f.read()
    f.close()

    suffix = filename.suffix
    parent = filename.parent.name
    extra_dir = filename.parents[1].name

    if suffix == ".json":
        return Data.from_json(identifier, json.loads(text), extra_dir, filename)

    elif suffix == ".xml":
        if extra_dir == "PMN" or extra_dir == "SOL":
            parent = f"{extra_dir}:{parent}"
        return Data.from_xml(identifier, et.fromstring(text), parent, filename)

    elif suffix == ".txt":
        if not genome:
            genome = ""
        return Data.from_text(identifier, text, parent, filename, genome)

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
    genome: Optional[str] = None,
    # NOTE: add docstring
):
    if isinstance(directory, str):
        directory = Path(directory).absolute()

    # Biocyc db families
    family = ""
    if database is not None and database.find(":") != -1:
        family, database = (i.rstrip().strip() for i in database.split(":"))
        family = family.upper()

        directory = directory.joinpath(family)

    # Try first locally and the query databases
    if not database:
        try:
            filename = next(directory.rglob(identifier + "*"))

        except StopIteration:
            response_database, response = get_response(identifier)

            if response_database == "bigg":
                extra: Path = getattr(response, "extra")

                if not extra:
                    raise AttributeError(
                        "The 'extra' attribute was not found in the response"
                    )
                directory = directory.joinpath(extra)

            filename = write(identifier, directory, response)
    else:
        # Create dir if needed
        if not directory.joinpath(database).exists():
            directory.joinpath(database).mkdir()

        try:
            if model_id:
                if not directory.joinpath(database, model_id).exists():
                    directory.joinpath(database, model_id).mkdir()

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
            if family:
                query = f"{family}:{database}:{identifier}"

            else:
                query = f"{database}:{identifier}"

            response_database, response = get_response(query, model_id)

            if genome:
                kegg.retrieve_kegg_genes(directory, genome)

            if database != "KEGG" and database != "BIGG":
                if not family:
                    biocyc.retrieve_gene_information(
                        directory, identifier, database
                    )

                elif family == "PMN":
                    plantcyc.retrieve_gene_information(
                        directory, identifier, database
                    )
                elif family == "SOL":
                    solcyc.retrieve_gene_information(
                        directory, identifier, database
                    )

            if response_database == "bigg":
                extra: Path = getattr(response, "extra")

                if not extra:
                    raise AttributeError(
                        "The 'extra' attribute was not found in the response"
                    )
                directory = directory.joinpath("BIGG", extra)

                if not directory.exists():
                    directory.mkdir()

            else:
                directory = directory.joinpath(database)

            filename = write(identifier, directory, response)

    data = file_to_Data_class(identifier, filename, genome)
    return data


def build_reaction_from_str(
    model: cobra_core.Model,
    reaction: cobra_core.Reaction,
    reaction_str: str,
    directory: Path,
    database: Optional[str],
    model_id: Optional[str],
    replacement: dict[str, str],
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

    # It must be a biocyc family
    if database is not None and database.find(":") != -1:
        directory = directory.parent

    for identifier, coefficient in compounds.items():
        compartment = identifier[-1]

        identifier = replacement.get(identifier[:-2], identifier[:-2])
        identifier = f"{identifier}_{compartment}"

        try:
            data = get_data(identifier[:-2], directory, database, model_id)
            metabolite = creation.metabolite_from_data(data, compartment, model)

        except (requests.HTTPError, ValueError):
            # NOTE: add log. This is only for custom
            metabolite = cobra_core.Metabolite(
                identifier, name=identifier, compartment=compartment
            )

        # NOTE: add logs
        if not isinstance(metabolite, cobra_core.Metabolite):
            raise TypeError("Given object is not a valid COBRApy Metabolite")

        reaction.add_metabolites({metabolite: coefficient})
