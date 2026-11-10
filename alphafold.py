# Download AlphaFold predictions and read the confidence score (pLDDT)
# of every amino acid.
#
# In AlphaFold's PDB files, pLDDT is stored in the column that is normally
# used for the "B-factor" (characters 61-66 of each ATOM line).

import json
import urllib.request

import pandas as pd

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def download(url):
    with urllib.request.urlopen(url) as response:
        return response.read().decode()


def read_plddt(pdb_text):
    """List of (amino acid, pLDDT), one per residue.
    Each residue has one 'CA' (alpha carbon) atom, so we only read those lines."""
    residues = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            aa = THREE_TO_ONE.get(line[17:20], "X")
            plddt = float(line[60:66])
            residues.append((aa, plddt))
    return residues


def get_protein(uniprot_id):
    """Ask the AlphaFold database for the protein's file, download it and read it."""
    info = json.loads(download("https://alphafold.ebi.ac.uk/api/prediction/" + uniprot_id))
    pdb_text = download(info[0]["pdbUrl"])
    return read_plddt(pdb_text)


def get_all_proteins(proteins):
    """proteins: dict uniprot_id -> name. Returns one table with every residue."""
    rows = []
    for uniprot_id, name in proteins.items():
        try:
            residues = get_protein(uniprot_id)
        except Exception as error:
            print("could not get", name, "-", error)
            continue
        print(name, len(residues), "residues")
        for aa, plddt in residues:
            rows.append({"protein": name, "aa": aa, "plddt": plddt})
    return pd.DataFrame(rows)
