import os
from pathlib import Path

from data_readers import TLC_TO_OLC

import pandas as pd

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

READ_PATHS = Path("pdb_files")


def modify_b_factor(z_table: pd.Series):

    pdb_names = list(filter(lambda x: x.endswith(".pdb"), os.listdir(READ_PATHS)))
    pdb_io = PDBIO()

    for name in pdb_names:

        pdb_id = name.replace(".pdb", "")
        protein: Model = PDBParser(QUIET=True).get_structure(pdb_id, str(READ_PATHS / name))[0]

        current_z_scores = z_table.filter(like=pdb_id)

        atom: Atom
        for atom in protein["S"].get_atoms():
            atom.bfactor = 0.0

        resi: Residue
        for resi in protein["L"].get_residues():

            resi_olc = TLC_TO_OLC[resi.resname]
            resi_num = resi.id[1]

            resi_abs_z_scores = current_z_scores.filter(like=f":{resi_olc}L{resi_num}").abs()

            for atom in resi.get_atoms():
                atom.bfactor = resi_abs_z_scores.mean() if len(resi_abs_z_scores) != 0 else 0.

        new_name = f"{pdb_id}_bfacmod.pdb"
        pdb_io.set_structure(protein)
        pdb_io.save(str(READ_PATHS / "outputs" / new_name))
