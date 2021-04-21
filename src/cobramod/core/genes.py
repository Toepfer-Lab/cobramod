#!/usr/bin/env python3
"""Modules genes

This module is responsable for the creation of COBRApy Genes from the parsed
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
        debug_log.info(
            f'For reaction "{reaction.id}", gene "{gene.id}" was created and '
            f'its name changed to "{gene.name}".'
        )
