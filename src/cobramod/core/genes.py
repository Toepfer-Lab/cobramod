#!/usr/bin/env python3
"""Modules genes

This module is responsable for the creation of COBRApy Genes from the parsed
data of the specific parsers.
"""
from cobra import Reaction


def _genes_to_reaction(reaction: Reaction, data_dict: dict):
    rule = data_dict["GENES"]["rule"] or ""
    reaction.gene_reaction_rule = rule
    for gene in reaction.genes:
        gene.name = data_dict["GENES"][gene.id]
        pass
