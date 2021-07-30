#!/usr/bin/env python3
"""Summary

This module is responsible for the short summary in the command
line and for other summaries in different formats.
"""
import warnings
from pathlib import Path
from typing import Dict

import pandas
from cobra import Model

from cobramod.debug import debug_log


class DataModel:
    """
    A class that can create a snapshot of a model and generates the summary
    based on those snapshots. Contains methods to identify differences and
    save them in multiple formats.
    """

    def __init__(self, lists: Dict[str, list]):
        self.metabolites = lists.get("metabolites", [])
        self.demands = lists.get("demands", [])
        self.exchanges = lists.get("exchanges", [])
        self.genes = lists.get("genes", [])
        self.groups = lists.get("groups", [])
        self.sinks = lists.get("sinks", [])

        # reaction includes only values not already in exchanges or sinks

        self.reactions = [
            reaction
            for reaction in lists.get("reactions") or []
            if reaction not in set(self.sinks)
            and reaction not in set(self.exchanges)
        ]

    @classmethod
    def from_model(cls, model: Model):
        """
        Method to create a DataModel object based on a model object.

        Args:
            model (Model): Model based on which a DataModel object
                is to be created.
        """

        characteristics = {
            "reactions": model.reactions,
            "metabolites": model.metabolites,
            "demands": model.demands,
            "exchanges": model.exchanges,
            "genes": model.genes,
            "groups": model.groups,
            "sinks": model.sinks,
        }

        data = {}

        for identifier, values in characteristics.items():
            try:
                data[identifier] = values.list_attr("id")
            except AttributeError:
                data[identifier] = []

        return cls(data)

    def diff(self, other):
        """
        Creates a new DataModel object consisting of
        the differences between the original and the
        passed DataModel object.

        Args:
            other (DataModel): DataModel to be compared with this one.
        """
        left = self - other
        right = other - self
        data = {
            "reactions": left.reactions + right.reactions,
            "metabolites": left.metabolites + right.metabolites,
            "demands": left.demands + right.demands,
            "exchanges": left.exchanges + right.exchanges,
            "genes": left.genes + right.genes,
            "groups": left.groups + right.groups,
            "sinks": left.sinks + right.sinks,
        }

        return DataModel(data)

    def __sub__(self, other):
        reactions = set(other.reactions)
        metabolites = set(other.metabolites)
        demands = set(other.demands)
        exchanges = set(other.exchanges)
        genes = set(other.genes)
        groups = set(other.groups)
        sinks = set(other.sinks)

        reactions = [x for x in self.reactions if x not in reactions]
        metabolites = [x for x in self.metabolites if x not in metabolites]
        demands = [x for x in self.demands if x not in demands]
        exchanges = [x for x in self.exchanges if x not in exchanges]
        genes = [x for x in self.genes if x not in genes]
        groups = [x for x in self.groups if x not in groups]
        sinks = [x for x in self.sinks if x not in sinks]

        data = {
            "reactions": reactions,
            "metabolites": metabolites,
            "demands": demands,
            "exchanges": exchanges,
            "genes": genes,
            "groups": groups,
            "sinks": sinks,
        }

        return DataModel(lists=data)

    def __str__(self):
        output = [
            "Reactions:\n",
            str(self.reactions),
            "\n",
            "Metabolites:\n",
            str(self.metabolites),
            "\n",
            "Exchange:\n",
            str(self.exchanges),
            "\n",
            "Demand:\n",
            str(self.demands),
            "\n",
            "Sinks:\n",
            str(self.sinks),
            "\n",
            "Genes:\n",
            str(self.genes),
            "\n",
            "Groups:\n",
            str(self.groups),
            "\n",
        ]

        return "".join(output)

    def _to_dataframe(
        self, model: Model = None, additions=None, deletions=None
    ):
        """
        Creates a pandas DataFrame based on a DataModel object.
        You can pass another DataModel object to include it in
        the DataFrame to have both DataModels in one DataFrame.
        Intended for internal use only.

        Args:
            model (Model): Model for the extraction of model id and name.
            additions (DataModel): DataModel that contains the new entities
                in the model.
            deletions (DataModel): DataModel that contains the removed entities
                in the model.
        """
        dtype = pandas.StringDtype()

        dictionary = {
            "Reactions": pandas.Series(self.reactions, dtype=dtype),
            "Exchange": pandas.Series(self.exchanges, dtype=dtype),
            "Demand": pandas.Series(self.demands, dtype=dtype),
            "Sinks": pandas.Series(self.sinks, dtype=dtype),
            "Metabolites": pandas.Series(self.metabolites, dtype=dtype),
            "Genes": pandas.Series(self.genes, dtype=dtype),
            "Groups": pandas.Series(self.groups, dtype=dtype),
        }

        if model is not None:
            model_identifier = {
                "Model identifier": pandas.Series(str(model.id), dtype=dtype),
                "Model name": pandas.Series(str(model.name), dtype=dtype),
            }

            # Preserving the order so that the original model
            # comes first and then the modifications
            model_identifier.update(dictionary)
            dictionary = model_identifier

        if additions is not None:
            addition = {
                "New in Reactions": pandas.Series(
                    additions.reactions, dtype=dtype
                ),
                "New in Exchange": pandas.Series(
                    additions.exchanges, dtype=dtype
                ),
                "New in Demand": pandas.Series(additions.demands, dtype=dtype),
                "New in Sinks": pandas.Series(additions.sinks, dtype=dtype),
                "New in Metabolites": pandas.Series(
                    additions.metabolites, dtype=dtype
                ),
                "New in Genes": pandas.Series(additions.genes, dtype=dtype),
                "New in Groups": pandas.Series(additions.groups, dtype=dtype),
            }
            dictionary.update(addition)

        if deletions is not None:
            deleted = {
                "Removed in Reactions": pandas.Series(
                    additions.reactions, dtype=dtype
                ),
                "Removed in Exchange": pandas.Series(
                    additions.exchanges, dtype=dtype
                ),
                "Removed in Demand": pandas.Series(
                    additions.demands, dtype=dtype
                ),
                "Removed in Sinks": pandas.Series(
                    additions.sinks, dtype=dtype
                ),
                "Removed in Metabolites": pandas.Series(
                    additions.metabolites, dtype=dtype
                ),
                "Removed in Genes": pandas.Series(
                    additions.genes, dtype=dtype
                ),
                "Removed in Groups": pandas.Series(
                    additions.groups, dtype=dtype
                ),
            }
            dictionary.update(deleted)

        return pandas.DataFrame(dictionary)

    def to_excl(
        self, path, model: Model = None, additions=None, deletions=None
    ):
        """
        Method to save a DataModel as an Excel file.
        Can also be used to save the changes between
        two points in time, in an Excel file format.
        For other formats see :func:`cobramod.summary.DataModel.to_csv`
        or :func:`cobramod.summary.DataModel.to_txt`.

        Args:
            path (Path): Location where the file is to be saved.
            model (Model): Model for the extraction of model id and name.
            additions (DataModel): DataModel that contains the new entities in
                the model.
            deletions (DataModel): DataModel that contains the remote entities
                in the model.
        """

        save = self._to_dataframe(model, additions, deletions)
        save.to_excel(path, index=False)

    def to_csv(
        self, path, model: Model = None, additions=None, deletions=None
    ):
        """
        Method to save a DataModel as a CSV file. Can also be used to save
        the changes between two points in time, as a CSV.
        For other formats see :func:`cobramod.summary.DataModel.to_excl` or
        :func:`cobramod.summary.DataModel.to_txt`.

        Args:
            path (Path): Location where the file is to be saved.
        model (Model): Model for the extraction of model id and name.
        additions (DataModel): DataModel that contains the new entities in
            the model.
        deletions (DataModel): DataModel that contains the remote entities
            in the model.
        """

        save = self._to_dataframe(model, additions, deletions)
        save.to_csv(path, index=False)

    def to_txt(
        self, path, model: Model = None, additions=None, deletions=None
    ):
        """
        Method to save a DataModel as a txt file. Can also be used to save
        the changes between two points in time, as a txt.
        For other formats see :func:`cobramod.summary.DataModel.to_excl` or
        :func:`cobramod.summary.DataModel.to_csv`.

        Args:
            path (Path): Location where the file is to be saved.
            model (Model): Model for the extraction of model id and name.
            additions (DataModel): DataModel that contains the new entities
            in the model.
            deletions (DataModel): DataModel that contains the remote entities
            in the model.
        """
        output = []

        if model is not None:
            output.extend(
                [
                    "Summary:",
                    f"Model identifier: {model.id}",
                    "Model name:",
                    str(model.name),
                ]
            )

        output.append(str(self))

        if additions is not None:
            output.append("New:")
            output.append(str(additions))

        if deletions is not None:
            output.append("Removed:")
            output.append(str(deletions))

        with open(file=str(path), mode="w") as file:
            file.writelines(line + "\n" for line in output)


def summary(model: Model, original: DataModel, filename: Path = None):
    """
    Produces the short summary and another one in the defined format.

    Args:
        model (Model): model with recent changes.
        original (DataModel): Object with data from the previous model. Use
            method :func:`cobramod.summary.DataModel`.
        filename (Path): Location where the summary should be stored.
            The file format is determined by the suffix of the filename.
            Thus, '.txt', '.csv' or '.xlsx' can be used.
    """

    if isinstance(filename, str):
        filename = Path(filename)

    new_values = DataModel.from_model(model)

    deletions = original - new_values
    additions = new_values - original

    # check if there are any changes otherwise notify the user
    num_changes = 0
    for value in vars(additions).values():
        num_changes += len(value)

    for value in vars(deletions).values():
        num_changes += len(value)

    if num_changes == 0:
        msg = "No changes in the model were detected!"
        debug_log.warning(msg)
        warnings.warn(msg, UserWarning)

    # output a short summary
    print(
        "{:13} {:^7} | {:^7} {:10}".format(
            "Number of", "new", "removed", "entities in"
        )
        + "\n"
        + "*"
        + ("=" * 21)
        + "|"
        + ("=" * 19)
        + "*"
        + "\n"
        "{:13} {:^7} | {:^7} {:10}".format(
            "Reactions", len(additions.reactions), len(deletions.reactions), ""
        )
        + "\n"
        "{:13} {:^7} | {:^7} {:10}".format(
            "Metabolites",
            len(additions.metabolites),
            len(deletions.metabolites),
            "",
        )
        + "\n"
        "{:13} {:^7} | {:^7} {:10}".format(
            "Exchange", len(additions.exchanges), len(deletions.exchanges), ""
        )
        + "\n"
        "{:13} {:^7} | {:^7} {:10}".format(
            "Demand", len(additions.demands), len(deletions.demands), ""
        )
        + "\n"
        "{:13} {:^7} | {:^7} {:10}".format(
            "Sinks", len(additions.sinks), len(deletions.sinks), ""
        )
        + "\n"
        "{:13} {:^7} | {:^7} {:10}".format(
            "Genes", len(additions.genes), len(deletions.genes), ""
        )
        + "\n"
        "{:13} {:^7} | {:^7} {:10}".format(
            "Groups", len(additions.groups), len(deletions.groups), ""
        )
        + "\n"
    )

    # if desired, create another summary in the requested format
    if filename is None:
        return

    options = {
        ".xlsx": new_values.to_excl,
        ".csv": new_values.to_csv,
        ".txt": new_values.to_txt,
    }

    try:
        file_format = filename.suffix
        options[file_format](filename, model, additions, deletions)
    except KeyError:
        debug_log.warning(
            "Unknown format therefore no summary was created."
            "Use '.xlsx', ',csv' or ',txt'."
        )
