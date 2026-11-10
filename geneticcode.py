# Functions for working with the genetic code.
# Codons are written with DNA letters (T instead of U).

import numpy as np
import pandas as pd

BASES = "TCAG"

# The standard genetic code (NCBI table 1). The letters are in the same order
# as the codons made below: TTT, TTC, TTA, TTG, TCT, ... GGG.  '*' means stop.
AMINO_ACID_LETTERS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"


def make_genetic_code():
    code = {}
    i = 0
    for first in BASES:
        for second in BASES:
            for third in BASES:
                codon = first + second + third
                code[codon] = AMINO_ACID_LETTERS[i]
                i += 1
    return code


GENETIC_CODE = make_genetic_code()

# the 20 amino acids (without stop), in alphabetical order
AMINO_ACIDS = sorted(set(AMINO_ACID_LETTERS) - {"*"})


def load_properties():
    return pd.read_csv("data/amino_acids.csv", index_col="aa")


def all_mutations():
    """Every way of changing one letter in a codon.
    Returns a list of (old_codon, new_codon, position)."""
    mutations = []
    for codon in GENETIC_CODE:
        for position in range(3):
            for base in BASES:
                if base != codon[position]:
                    new_codon = codon[:position] + base + codon[position + 1:]
                    mutations.append((codon, new_codon, position + 1))
    return mutations


MUTATIONS = all_mutations()  # 64 codons x 9 = 576, only needs making once


def error_cost(code, values):
    """Average squared change in an amino acid property over all single-letter
    mutations. Mutations to or from a stop codon are skipped.

    code   : dict codon -> amino acid letter
    values : dict amino acid letter -> number (e.g. polar requirement)
    """
    before = []
    after = []
    for old_codon, new_codon, position in MUTATIONS:
        old_aa = code[old_codon]
        new_aa = code[new_codon]
        if old_aa == "*" or new_aa == "*":
            continue
        before.append(values[old_aa])
        after.append(values[new_aa])

    before = np.array(before)
    after = np.array(after)
    return np.mean((before - after) ** 2)


def random_code():
    """A fake genetic code: same groups of codons as the real one,
    but the 20 amino acids are shuffled between the groups."""
    shuffled = np.random.permutation(AMINO_ACIDS)
    swap = dict(zip(AMINO_ACIDS, shuffled))
    swap["*"] = "*"  # stop codons stay stop codons

    fake = {}
    for codon in GENETIC_CODE:
        fake[codon] = swap[GENETIC_CODE[codon]]
    return fake


def damage_per_amino_acid(values):
    """For each amino acid: the average squared change caused by a mutation
    in one of its codons (using the real code). Low = well protected."""
    changes = {}  # amino acid -> list of squared changes
    for old_codon, new_codon, position in MUTATIONS:
        old_aa = GENETIC_CODE[old_codon]
        new_aa = GENETIC_CODE[new_codon]
        if old_aa == "*" or new_aa == "*":
            continue
        if old_aa not in changes:
            changes[old_aa] = []
        changes[old_aa].append((values[old_aa] - values[new_aa]) ** 2)

    damage = {}
    for aa in changes:
        damage[aa] = np.mean(changes[aa])
    return damage
