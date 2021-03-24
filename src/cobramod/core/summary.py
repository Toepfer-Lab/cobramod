import os
from pathlib import Path
from typing import Dict

import pandas
from cobra import Model


class DataModel:
    """
    A class that can create a snapshot of a model and generate the summary
    based on those snapshots. Contains methods to identify differences
    and save them in multiple formats.
    """

    def __init__(self, lists: Dict[str, list]):
        self.metabolites = lists.get("metabolites")
        self.demands = lists.get("demands")
        self.exchanges = lists.get("exchanges")
        self.genes = lists.get("genes")
        self.groups = lists.get("groups")
        self.sinks = lists.get("sinks")

        # reaction includes only values not already in exchanges or sinks

        self.reactions = [reaction for reaction in lists.get("reactions") or [] if
                          reaction not in set(self.sinks) and reaction not in set(self.exchanges)]

    @classmethod
    def from_model(cls, model: Model):
        """
        Method to create a DataModel object based on a model object.

        Args:
            model (Model): Model based on which a DataModel object is to be created.
        """
        data = {'reactions': model.reactions.list_attr("id"),
                'metabolites': model.metabolites.list_attr("id"),
                'demands': model.demands.list_attr("id"),
                'exchanges': model.exchanges.list_attr("id"),
                'genes': model.genes.list_attr("id"),
                'groups': model.groups.list_attr("id"),
                'sinks': model.sinks.list_attr("id")}
        return cls(data)

    def diff(self, other):
        """
        Creates a new DataModel object consisting of the differences between the
        original and the passed DataModel object.

        Args:
            other (DataModel): DataModel to be compared with this one.
        """
        left = self - other
        right = other - self
        data = {'reactions': left.reactions + right.reactions,
                'metabolites': left.metabolites + right.metabolites,
                'demands': left.demands + right.demands,
                'exchanges': left.exchanges + right.exchanges,
                'genes': left.genes + right.genes,
                'groups': left.groups + right.groups,
                'sinks': left.sinks + right.sinks}

        return DataModel(data)

    def __sub__(self, other):
        data = {'reactions': [x for x in self.reactions if x not in set(other.reactions)],
                'metabolites': [x for x in self.metabolites if x not in set(other.metabolites)],
                'demands': [x for x in self.demands if x not in set(other.demands)],
                'exchanges': [x for x in self.exchanges if x not in set(other.exchanges)],
                'genes': [x for x in self.genes if x not in set(other.genes)],
                'groups': [x for x in self.groups if x not in set(other.groups)],
                'sinks': [x for x in self.sinks if x not in set(other.sinks)]}

        return DataModel(lists=data)

    def __str__(self):
        output = ['Reactions:\n',
                  str(self.reactions), "\n",
                  "Metabolites:\n",
                  str(self.metabolites), "\n",
                  "Exchange:\n",
                  str(self.exchanges), "\n",
                  "Demand:\n",
                  str(self.demands), "\n",
                  "Sinks:\n",
                  str(self.sinks), "\n",
                  "Genes:\n",
                  str(self.genes), "\n",
                  "Groups:\n",
                  str(self.groups), "\n"]

        return ''.join(output)

    def _to_dataframe(self, model: Model = None, original=None):
        """
        Creates a pandas DataFrame based on a DataModel object. You can pass another DataModel
        object to include it in the DataFrame to have both DataModels in one DataFrame.
        Intended for internal use only.

        Args:
            model (Model): model for the extraction of model id and name.
            original (DataModel): Object with data from previous model. Use
                method :func:`cobramod.summary.DataModel`.
        """

        dictionary = {
            "Changed reactions": pandas.Series(self.reactions),
            "Changed exchange": pandas.Series(self.exchanges),
            "Changed demand": pandas.Series(self.demands),
            "Changed sinks": pandas.Series(self.sinks),
            "Changed metabolites": pandas.Series(self.metabolites),
            "Changed genes": pandas.Series(self.genes),
            "Changed groups": pandas.Series(self.groups),
        }

        if model is not None:
            original = {
                "Model identifier": pandas.Series(str(model.id)),
                "Model name": pandas.Series(str(model.name)),
                "Reactions": pandas.Series(original.reactions),
                "Exchange": pandas.Series(original.exchanges),
                "Demand": pandas.Series(original.demands),
                "Sinks": pandas.Series(original.sinks),
                "Metabolites": pandas.Series(original.metabolites),
                "Genes": pandas.Series(original.genes),
                "Groups": pandas.Series(original.groups),
            }

            # Preserving the order so that the original model comes first and then the modifications
            original.update(dictionary)
            dictionary = original

        return pandas.DataFrame(dictionary)

    def to_excl(self, path, model: Model = None, original=None):
        """
        Method to save a DataModel as an Excel file. Can also be used to save the
        changes between two points in time, in an Excel file format.
        For other formats see :func:`cobramod.summary.DataModel.to_csv` or
        :func:`cobramod.summary.DataModel.to_txt`.

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
                    str(self),
                    "Changes:",
                ]
            )

        output.append(str(self.diff(original)))

        with open(file=str(path), mode="w+") as f:
            f.writelines(line + "\n" for line in output)


def summary(model: Model,
            original,
            file_format: str = None,
            filename: Path = None
            ):
    """
        Produces the short summary and another one in the defined format.

        Args:
            model (Model): model with recent changes.
            original (DataModel): Object with data from previous model. Use
                method :func:`cobramod.summary.DataModel`.
            file_format (str): The format in which the further summary should be generated.
                If no additional summary is desired, this should be None.
            filename (Path): Location where the summary should be stored.
    """
    if filename is None:
        filename = Path(os.getcwd()) / "summary"

    new_values = DataModel.from_model(model)

    diff = original.diff(new_values)
    print("Changes:\n"
          "Reactions \t" + str(len(diff.reactions)) + "\n"
          "Metabolites\t" + str(len(diff.metabolites)) + "\n"
          "Exchange \t" + str(len(diff.exchanges)) + "\n"
          "Demand\t\t" + str(len(diff.demands)) + "\n"
          "Sinks\t\t" + str(len(diff.sinks)) + "\n"
          "Genes\t\t" + str(len(diff.genes)) + "\n"
          "Groups\t\t" + str(len(diff.groups)) + "\n")

    if file_format is None:
        return
    elif file_format is "excel":
        diff.to_excl(filename.with_suffix(".xlsx"), model, new_values)
    elif file_format is "csv":
        diff.to_csv(filename.with_suffix(".csv"), model, new_values)
    elif file_format is "txt":
        diff.to_txt(filename.with_suffix(".txt"), model, new_values)
    else:
        print("No known format. Use 'excel', 'csv' or 'txt'")
