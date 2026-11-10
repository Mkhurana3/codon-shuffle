# Run with:  pytest   (from the main project folder)

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geneticcode import GENETIC_CODE, MUTATIONS, error_cost, load_properties, random_code
from alphafold import read_plddt


def test_genetic_code():
    assert len(GENETIC_CODE) == 64
    assert GENETIC_CODE["ATG"] == "M"   # start codon = methionine
    assert GENETIC_CODE["TAA"] == "*"   # stop


def test_number_of_mutations():
    assert len(MUTATIONS) == 64 * 9


def test_same_answer_as_the_paper():
    # Haig & Hurst (1991) got 5.194 for the real code
    values = load_properties()["polar_requirement"].to_dict()
    assert round(error_cost(GENETIC_CODE, values), 3) == 5.194


def test_random_code_keeps_stop_codons():
    fake = random_code()
    for codon in ["TAA", "TAG", "TGA"]:
        assert fake[codon] == "*"


def test_read_plddt():
    # two lines in the same column layout as an AlphaFold PDB file
    line1 = "ATOM      2  CA  MET A   1       0.000   0.000   0.000  1.00 45.10"
    line2 = "ATOM      9  CA  LYS A   2       0.000   0.000   0.000  1.00 91.37"
    assert read_plddt(line1 + "\n" + line2) == [("M", 45.10), ("K", 91.37)]
