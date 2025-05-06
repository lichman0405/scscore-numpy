"""
Utility functions for SCScorer.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

FP_LEN = 1024
FP_RAD = 2

def smiles_to_fingerprint(smiles):
    """
    Convert a SMILES string to a molecular fingerprint.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros((FP_LEN,), dtype=bool)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(
        mol, FP_RAD, nBits=FP_LEN, useChirality=True), dtype=bool)