"""Modules genes

This module is responsible for the creation of COBRApy Genes from the parsed
data of the specific parsers.
"""

from typing import Any

import cobra.core as cobra_core

from cobramod.debug import debug_log


def genes_to_reaction(reaction: cobra_core.Reaction, data_dict: dict[str, Any]):
    """
    Adds the corresponding gene-reaction rule to the Reaction and correct
    its name according to the passed dictionary
    """
    rule: str = data_dict["rule"]

    # This line will create the corresponding genes of the reaction
    reaction.gene_reaction_rule = rule

    # Update name
    for gene in reaction.genes:
        gene.name = data_dict["genes"][gene.id]

        if not gene.name:
            gene.name = gene.id

        debug_log.debug(
            f'Gene "{gene.id}" was added with gene name "{gene.name}" to '
            f'reaction "{reaction.id}".'
        )
