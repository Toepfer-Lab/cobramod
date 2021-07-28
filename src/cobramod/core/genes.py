#!/usr/bin/env python3
"""Modules genes

This module is responsible for the creation of COBRApy Genes from the parsed
data of the specific parsers.
"""
from cobra import Reaction

from cobramod.debug import debug_log


def _genes_to_reaction(reaction: Reaction, data_dict: dict):
    """
    Adds the corresponding gene-reaction rule to the Reaction and correct
    its name according to the passed dictionary
    """
    rule = data_dict["GENES"]["rule"]
    # This line will create the corresponding genes of the reaction
    reaction.gene_reaction_rule = rule
    # Update name
    for gene in reaction.genes:
        gene.name = data_dict["GENES"]["genes"][gene.id]
        if not gene.name:
            gene.name = gene.id
        debug_log.info(
            f'Gene "{gene.id}" was added with gene name "{gene.name}" to '
            f'reaction "{reaction.id}".'
        )
