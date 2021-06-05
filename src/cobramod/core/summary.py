#!/usr/bin/env python3
"""Summary

This module is responsible for the short summary in the command
line and for other summaries in different formats.
"""
import logging
import warnings
from pathlib import Path
from typing import Dict

import pandas
from cobra import Model

from cobramod.debug import debug_log


class DataModel:
    """
    A class that can create a snapshot of a model and generate the summary
    based on those snapshots. Contains methods to identify differences
    and save them in multiple formats.
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
        data = {
            "reactions": model.reactions.list_attr("id"),
            "metabolites": model.metabolites.list_attr("id"),
            "demands": model.demands.list_attr("id"),
            "exchanges": model.exchanges.list_attr("id"),
            "genes": model.genes.list_attr("id"),
            "groups": model.groups.list_attr("id"),
            "sinks": model.sinks.list_attr("id"),
        }
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

    def _to_dataframe(self, model: Model = None, original=None):
        """
        Creates a pandas DataFrame based on a DataModel object.
        You can pass another DataModel object to include it in
        the DataFrame to have both DataModels in one DataFrame.
        Intended for internal use only.

        Args:
            model (Model): model for the extraction of model id and name.
            original (DataModel): Object with data from previous model. Use
                method :func:`cobramod.summary.DataModel`.
        """
        dtype = pandas.StringDtype()

        dictionary = {
            "Changed reactions": pandas.Series(self.reactions, dtype=dtype),
            "Changed exchange": pandas.Series(self.exchanges, dtype=dtype),
            "Changed demand": pandas.Series(self.demands, dtype=dtype),
            "Changed sinks": pandas.Series(self.sinks, dtype=dtype),
            "Changed metabolites": pandas.Series(
                self.metabolites, dtype=dtype
            ),
            "Changed genes": pandas.Series(self.genes, dtype=dtype),
            "Changed groups": pandas.Series(self.groups, dtype=dtype),
        }

        if model is not None:
            original = {
                "Model identifier": pandas.Series(str(model.id), dtype=dtype),
                "Model name": pandas.Series(str(model.name), dtype=dtype),
                "Reactions": pandas.Series(original.reactions, dtype=dtype),
                "Exchange": pandas.Series(original.exchanges, dtype=dtype),
                "Demand": pandas.Series(original.demands, dtype=dtype),
                "Sinks": pandas.Series(original.sinks, dtype=dtype),
                "Metabolites": pandas.Series(
                    original.metabolites, dtype=dtype
                ),
                "Genes": pandas.Series(original.genes, dtype=dtype),
                "Groups": pandas.Series(original.groups, dtype=dtype),
            }

            # Preserving the order so that the original model
            # comes first and then the modifications
            original.update(dictionary)
            dictionary = original

        return pandas.DataFrame(dictionary)

    def to_excl(self, path, model: Model = None, original=None):
        """
        Method to save a DataModel as an Excel file.
        Can also be used to save the changes between
        two points in time, in an Excel file format.
        For other formats see :func:`cobramod.summary.DataModel.to_csv`
        or :func:`cobramod.summary.DataModel.to_txt`.

        Args:
            path (Path): Location where the file is to be saved.
            model (Model): model for the extraction of model id and name.
            original (DataModel): DataModel of the original model
        """

        save = self._to_dataframe(model, original)
        save.to_excel(path, index=False)

    def to_csv(self, path, model: Model = None, original=None):
        """
        Method to save a DataModel as a csv file. Can also be used to save
        the changes between two points in time, as a csv.
        For other formats see :func:`cobramod.summary.DataModel.to_excl` or
        :func:`cobramod.summary.DataModel.to_txt`.

        Args:
            path (Path): Location where the file is to be saved.
            model (Model): model for the extraction of model id and name.
            original (DataModel): DataModel of the original model
        """

        save = self._to_dataframe(model, original)
        save.to_csv(path, index=False)

    def to_txt(self, path, model: Model = None, original=None):
        """
        Method to save a DataModel as a txt file. Can also be used to save
        the changes between two points in time, as a txt.
        For other formats see :func:`cobramod.summary.DataModel.to_excl` or
        :func:`cobramod.summary.DataModel.to_csv`.

        Args:
            path (Path): Location where the file is to be saved.
            model (Model): model for the extraction of model id and name.
            original (DataModel): DataModel of the original model
        """
        output = []

        if model is not None:
            output.extend(
                [
                    "Summary:",
                    f"Model identifier: {model.id}",
                    "Model name:",
                    str(model.name),
                    str(original),
                    "Changes:",
                ]
            )

        output.append(str(self))

        with open(file=str(path), mode="w") as file:
            file.writelines(line + "\n" for line in output)


def summary(
    model: Model, original: DataModel, filename: Path = None,
):
    """
        Produces the short summary and another one in the defined format.

        Args:
            model (Model): model with recent changes.
            original (DataModel): Object with data from previous model. Use
                method :func:`cobramod.summary.DataModel`.
            filename (Path): Location where the summary should be stored. The
                file format is determined by the suffix of the filename.
                Thus, '.txt', '.csv' or '.xlsx' can be used.
    """

    if isinstance(filename, str):
        filename = Path(filename)

    new_values = DataModel.from_model(model)

    diff = original.diff(new_values)

    # check if there are any changes otherwise notify the user
    num_changes = 0
    for value in vars(diff).values():
        num_changes += len(value)

    if num_changes == 0:
        msg = "No changes in the model were detected!"
        debug_log.warning(msg)
        warnings.warn(msg, UserWarning)

    # output a short summary
    print(
        "Changes:\n"
        "Reactions \t" + str(len(diff.reactions)) + "\n"
        "Metabolites\t" + str(len(diff.metabolites)) + "\n"
        "Exchange \t" + str(len(diff.exchanges)) + "\n"
        "Demand\t\t" + str(len(diff.demands)) + "\n"
        "Sinks\t\t" + str(len(diff.sinks)) + "\n"
        "Genes\t\t" + str(len(diff.genes)) + "\n"
        "Groups\t\t" + str(len(diff.groups)) + "\n"
        "(Includes additions & deletions)\n"
    )

    # if desired, create another summary in the requested format
    if filename is None:
        return

    options = {
        ".xlsx": diff.to_excl,
        ".csv": diff.to_csv,
        ".txt": diff.to_txt,
    }

    try:
        file_format = filename.suffix
        options[file_format](filename, model, new_values)
    except KeyError:
        logging.warning(
            "Unknown format therefore no summary was created."
            "Use '.xlsx', ',csv' or ',txt'."
        )
