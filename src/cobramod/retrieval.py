"""
Retrieval of data from multiple databases

This module includes the class 'Data', that is responsable for the internal
data representation for this package. The class Data includes methods for
parsing itself into a Reaction, Metabolite or Pathway.

Additionally, this package includes high level functions such as `get_data`,
`file_to_Data_class` or `get_response` that handles the data retrieval or
convertion
"""

from __future__ import annotations

import json
import urllib.parse
import warnings
import xml.etree.ElementTree as et
from pathlib import Path
from typing import Any, Literal, Optional, Union

import cobra.core as cobra_core
import requests

import cobramod.utils as cmod_utils
from cobramod.core import creation, genes
from cobramod.debug import debug_log
from cobramod.parsing import bigg, biocyc, kegg, plantcyc, solcyc
from cobramod.parsing import db_version as cmod_db

db_configuration = cmod_db.DataVersionConfigurator()

BIOCYC = "https://websvc.biocyc.org/getxml?id="
PMN = "https://pmn.plantcyc.org/getxml?id="
KEGG = "https://rest.kegg.jp/get/"
# NOTE: deprecation for v2
SOLCYC = "https://solcyc.sgn.cornell.edu/getxml?id="

EXTENSIONS = [".xml", ".txt", ".json"]


class Databases:
    """
    Simple object that shows the information about the database
    """

    def __init__(self):
        self.msg = (
            "CobraMod supports BioCyc, the Plant Metabolic Network (PMN), "
            "KEGG and BiGG Models repository. "
            "BioCyc includes around 18.000 sub-databases. A complete list "
            "for BioCyc can be found at: "
            "'https://biocyc.org/biocyc-pgdb-list.shtml'.\n"
            "The database-specific identifiers can be found in the URL of the "
            "respective data."
        )

    def __repr__(self):
        """
        Returns a string about the definition of "identifier" and where
        to obtain them.
        """
        self.msg2 = (
            "CobraMod uses abbreviations to represent the databases or "
            "sub-databases:\n"
            "Abbreviation   Database Name\n"
            "META           Biocyc, subdatabase MetaCyc\n"
            "ARA            Biocyc, subdatabase AraCyc\n"
            "\n"
            "KEGG           Kyoto encyclopedia of Genes and Genomes\n"
            "BIGG           BiGG Models\n"
            "\n"
            "PMN:META       Plantcyc, subdatabase META\n"
            "PMN:ARA        Plantcyc, subdatabase ARA\n"
            "This applies for all subdatabases from SolCyc, BioCyc and Plantcyc\n"
        )

        return self.msg + self.msg2

    def _repr_html_(self):
        """
        Returns a HTML string about the definition of "identifier" and where
        to obtain them.
        """
        return """
<p>CobraMod supports BioCyc, the Plant Metabolic Network (PMN), KEGG and BiGG Models repository. BioCyc includes around 18.000 sub-databases. A complete list for BioCyc can be found at: 'https://biocyc.org/biocyc-pgdb-list.shtml'.<br
    />The database-specific identifiers can be found in the URL of the respective data.CobraMod uses abbreviations to represent the databases or sub-databases:
</p>
<table style="width: 50%; border-collapse: collapse;" border="1">
    <tbody>
        <tr>
            <td style="width:
        50%;"><strong>Abbreviation</strong></td>
            <td style="width:
            50%;"><strong>Database name</strong></td>
        </tr>
        <tr>
            <td style="width: 50%;">META</td>
            <td style="width: 50%;">Biocyc, subdatabase MetaCyc</td>
        </tr>
        <tr>
            <td style="width:
                50%;">ARA</td>
            <td style="width: 50%;">Biocyc, subdatabase AraCyc
            </td>
        </tr>
        <tr>
            <td style="width: 50%;">KEGG</td>
            <td style="width: 50%;">Kyoto encyclopedia of Genes and Genomes
            </td>
        </tr>
        <tr>
            <td style="width: 50%;">BIGG</td>
            <td style="width: 50%;">BiGG Models</td>
        </tr>
        <tr>
            <td style="width: 50%;">PMN:META</td>
            <td style="width:
                    50%;">PlantCyc, subdatabase META</td>
        </tr>
        <tr>
            <td style="width: 50%;">PMN:ARA</td>
            <td style="width:
                        50%;">PlantCyc, subdatabase ARA</td>
        </tr>
    </tbody>
</table>
<p>This applies for all subdatabases from BioCyc and Plantcyc
</p>
        """


available_databases = Databases()

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
    version: str
    """
    CobraMod's Data object. This class is responsible for converting
    information into COBRApy objects
    """

    def __init__(
        self,
        identifier: str,
        attributes: dict[str, str],
        mode: Literal["Pathway", "Reaction", "Metabolite"],
        database: str,
        path: Path,
        model_id: str,
        version: str,
    ):
        self.mode = mode
        self.identifier = identifier
        self.attributes = attributes
        self.database = database
        self.path = path
        self.model_id = model_id
        self.genome = ""
        self.version = version

    def __repr__(self):
        return f"Data [{self.mode}: {self.identifier}]"

    @classmethod
    def from_json(
        cls, entry: str, data: dict[str, Any], database: str, path: Path
    ):
        """
        Creates an instance from a JSON file. Database: BIGG
        """
        is_compound = data.get("formulae", data.get("formula", None))
        model_id = path.parent.name

        if path.parents[1].joinpath("database_version").exists():
            with open(path.parents[1].joinpath("database_version"), "r") as f:
                version = json.load(f).get("bigg_models_version")
        else:
            info_response = requests.get(
                "http://bigg.ucsd.edu/api/v2/database_version"
            )
            info_response.raise_for_status()

            with (
                path.parents[1]
                .joinpath("database_version")
                .open(mode="w+") as f
            ):
                f.write(info_response.text)
            version = info_response.json().get("bigg_models_version")

        if version is None:
            raise Exception("Version of the BIGG JSON cannot be empty!")

        # Non universal models includes suffices
        if model_id != "universal":
            entry = entry[: entry.find("_")]

        mode: Literal["Pathway", "Reaction", "Metabolite"]
        if is_compound:
            mode = "Metabolite"
            attributes = bigg.parse_metabolite_attributes(data)
        else:
            mode = "Reaction"
            attributes = bigg.parse_reaction_attributes(data, entry)

        return cls(entry, attributes, mode, database, path, model_id, version)

    @classmethod
    def from_text(
        cls, entry: str, text: str, database: str, path: Path, genome: str
    ):
        """
        Creates an instance from plain text. Database: KEGG
        """
        data: dict[str, list[str]] = kegg.data_from_string(text)
        entry_mode = data["ENTRY"][0].split()[-1]
        gene_path = path.parent.joinpath("GENES")

        if path.parent.joinpath("database_version").exists():
            with open(path.parent.joinpath("database_version"), "r") as f:
                version = cmod_utils.kegg_info_to_version(f.read())
        else:
            info_response = requests.get("http://rest.kegg.jp/info/kegg")
            info_response.raise_for_status()

            with path.parent.joinpath("database_version").open(mode="w+") as f:
                f.write(info_response.text)
            version = cmod_utils.kegg_info_to_version(info_response.text)

        if version is None:
            raise Exception("Version of the BIGG JSON cannot be empty!")

        mode: Literal["Pathway", "Reaction", "Metabolite"]
        if entry_mode == "Compound":
            mode = "Metabolite"
            attributes = kegg.parse_metabolite_attributes(data, entry)

        elif entry_mode == "Reaction":
            mode = "Reaction"
            attributes = kegg.parse_reaction_attributes(
                data, entry, genome, gene_path
            )

        elif entry_mode == "Module":
            mode = "Pathway"
            attributes = kegg.parse_pathway_attributes(data, entry)
        else:
            raise AttributeError(f"Cannot parse type '{entry_mode}'")

        obj = cls(entry, attributes, mode, database, path, "", version)
        obj.genome = genome
        return obj

    @classmethod
    def from_xml(cls, entry: str, root: et.Element, database: str, path: Path):
        """
        Creates an instance from an XML file. Database: Biocyc families
        """
        version_element = root.find("metadata/PGDB")
        if version_element is not None:
            version: Optional[str] = version_element.attrib.get("version")
            assert version is not None
        else:
            raise Exception("Version for the XML cannot be empty!")

        gene_path = path.parent.joinpath("GENES")

        mode: Literal["Pathway", "Reaction", "Metabolite"]
        if (
            root.findall("Compound")
            or root.findall("Protein")
            or root.findall("RNA")
        ):
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
        return cls(entry, attributes, mode, database, path, "", version)

    def parse(
        self,
        model: cobra_core.Model,
        compartment: str,
        replacement: dict[str, str] = {},
    ) -> cobra_core.Object:
        """
        Parses information and creates a cobrapy Object. This can be a
        Metabolite or Reaction.

        Args:
            model (cobra.core.Model): Model to search for information
            compartment (str): Location of the object
            replacement (dict[str, str]): In case of Reactions, the dictionary
                represents the metabolites to replace. The key is the
                identifier of the object to replaced and the value is the new
                identifier to use instead. Defaults to an empty dictionary

        Returns:
            cobra.core.Object This can be a Reaction or Metabolites

        """
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
        encoded_id = urllib.parse.quote(identifier, safe="")

        if database == "KEGG":
            url = f"{databases['KEGG']}{encoded_id}"
        else:
            url = f"{databases['BIOCYC']}{database}:{encoded_id}"
            biocyc_credentials = True

    if database == "BIGG":
        query = query.removeprefix("BIGG:")
        response, _ = bigg.find_url(model_id, query)
        response.raise_for_status()

        # This extra attribute is later used to save the file in their
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
    """
    Creates a Data object from given file.

    Args:
        identifier (str): Identifier of the object. It should be included in
            the file itself
        filename (Path): Location of the file
        genome (Optional[str]): BIGG-specific argument. Genome involved

    Returns:
        Data
    """
    f = open(filename, "r")
    text = f.read()
    f.close()

    suffix = filename.suffix
    parent = filename.parent.name
    extra_dir = filename.parents[1].name

    if suffix == ".json":
        return Data.from_json(identifier, json.loads(text), extra_dir, filename)

    elif suffix == ".xml":
        if extra_dir == "PMN":
            parent = f"{extra_dir}:{parent}"
        # NOTE: Remove for v2.0.0
        elif extra_dir == "SOL":
            parent = f"{extra_dir}:{parent}"
            warnings.warn(
                "Database Solcyc is being deprecated for next version",
                DeprecationWarning,
            )
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
    """
    Returns a Generator with the possible files given a specific identifier
    """
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
) -> Data:
    """
    Retrieves the Data for given identifier. This function either retrieves
    from the server of the databases or locally.

    Args:
        identifier (str): Name of the object to retrieve
        dictionary (str or Path): Location of the files to retrieve or store
        database (Optional[str]): Name of the database. Check
            cobramod.retrieval.available_databases for more information
        model_id (Optional[str]): BIGG-specific argument. Name of the model to
            retrieve information
        genome: (Optional[str]): Name of the genome to retrieve

    Returns:
        Data
    """
    if isinstance(directory, str):
        directory = Path(directory).absolute()

    main_dir = directory

    # Biocyc db families
    family = ""
    if database is not None and database.find(":") != -1:
        family, database = (i.rstrip().strip() for i in database.split(":"))
        family = family.upper()

        directory = directory.joinpath(family)

    if not directory.exists():
        directory.mkdir()

    extra: Path
    # Try first locally and the query databases
    if not database:
        try:
            filename = next(directory.rglob(identifier + "*"))
            response_database = filename.parent.name.upper()

        except StopIteration:
            response_database, response = get_response(identifier)

            if response_database == "bigg":
                extra = getattr(response, "extra")

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
                response_database = "BIGG"

            else:
                filename = next(
                    get_files(directory.joinpath(database), identifier)
                )
                response_database = filename.parent.name.upper()

        except StopIteration:
            if family:
                query = f"{family}:{database}:{identifier}"

            else:
                query = f"{database}:{identifier}"

            response_database, response = get_response(query, model_id)

            if database == "KEGG":
                kegg.retrieve_kegg_genes(directory, identifier)

            if database != "KEGG" and database != "BIGG":
                if not family:
                    biocyc.retrieve_gene_information(
                        directory, identifier, database
                    )

                elif family == "PMN":
                    plantcyc.retrieve_gene_information(
                        directory, identifier, database
                    )
                # NOTE: Deprecation for 2.0.0
                elif family == "SOL":
                    warnings.warn(
                        "Database Solcyc is being deprecated for next version",
                        DeprecationWarning,
                    )
                    solcyc.retrieve_gene_information(
                        directory, identifier, database
                    )

            if response_database == "BIGG":
                extra = getattr(response, "extra")

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

    if family:
        response_database = f"{family}:{database}"

    db_configuration.check_database_version(
        main_dir, response_database, data.version
    )
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
    """
    Populates the reaction with Metabolites from given reaction string.
    Metabolites are retrieved from the servers or locally

    Args:
        model (cobra.core.Model): Model to search for metabolites
        reaction (cobra.core.Reaction): Empty reaction to populate
        reaction_str (str): valid reaction string
        directory (Path): Path for the location of the directory where the
            information is stored
        database (str, optional): Name of the database. Check
            cobramod.retrieval.available_databases for more information
        model_id (str, optional): BIGG-specific argument. Name of the model to
            retrieve information
        genome: (Optional[str]): Name of the genome to retrieve
        replacement (dict[str, str]): The dictionary represents the metabolites
            to replace. The key is the identifier of the object to replaced and
            the value is the new identifier to use instead
    """
    position = cmod_utils.get_arrow_position(reaction_str)
    arrow = reaction_str[position : position + 3]
    reaction.bounds = cmod_utils.ARROWS[arrow]

    compounds: dict[str, float] = {}
    # Reactants
    for pair in reaction_str[:position].rstrip().strip().split("+"):
        coefficient: Union[str, float]
        try:
            coefficient, metabolite_str = pair.rstrip().strip().split(" ")
        except ValueError:
            coefficient = 1
            metabolite_str = pair.rstrip().strip()

            # Empty part
            if metabolite_str == "":
                break

        compounds[metabolite_str] = -1.0 * float(coefficient)

    # products
    for pair in reaction_str[position + 3 :].rstrip().strip().split("+"):
        try:
            coefficient, metabolite_str = pair.rstrip().strip().split(" ")
        except ValueError:
            coefficient = 1
            metabolite_str = pair.rstrip().strip()

            # Empty part
            if metabolite_str == "":
                break

        compounds[metabolite_str] = float(coefficient)

    # It must be a biocyc family
    # if database is not None and database.find(":") != -1:
    #     directory = directory.parent

    for identifier, coefficient in compounds.items():
        compartment = identifier[-1]

        identifier = replacement.get(identifier[:-2], identifier[:-2])
        identifier = f"{identifier}_{compartment}"

        metabolite: cobra_core.Metabolite
        try:
            data = get_data(identifier[:-2], directory, database, model_id)
            metabolite = creation.metabolite_from_data(data, compartment, model)

        except (requests.HTTPError, ValueError):
            metabolite = cobra_core.Metabolite(
                identifier, name=identifier, compartment=compartment
            )
            debug_log.debug(f"Curated Metabolite '{metabolite.id}' created")

        # NOTE: add logs
        if not isinstance(metabolite, cobra_core.Metabolite):
            raise TypeError("Given object is not a valid COBRApy Metabolite")

        reaction.add_metabolites({metabolite: coefficient})
        debug_log.debug(
            f"Metabolite '{metabolite.id}' added to reaction '{reaction.id}'"
        )
