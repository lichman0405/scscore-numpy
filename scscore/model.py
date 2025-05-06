"""
SCScorer: A standalone NumPy-based model for scoring molecular synthetic complexity.

This module provides functionality to load a pre-trained model and perform scoring on given SMILES strings.
"""

import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle

# Constants
FP_LEN = 2048  # Updated fingerprint length to match model weight dimensions
FP_RAD = 2  # Fingerprint radius
SCORE_SCALE = 5.0  # Score scaling factor

class SCScorer:
    def __init__(self, score_scale=SCORE_SCALE):
        self.vars = []
        self.score_scale = score_scale
        self._restored = False

    def restore(self):
        """
        Load model weights from a pre-defined path within the project.
        """
        weight_path = os.path.join(os.path.dirname(__file__), 'models', 'model.ckpt-10654.as_numpy.pickle')
        self._load_vars(weight_path)

        def mol_to_fp(mol):
            if mol is None:
                return np.zeros((FP_LEN,), dtype=bool)
            return np.array(AllChem.GetMorganFingerprintAsBitVect(
                mol, FP_RAD, nBits=FP_LEN, useChirality=True), dtype=bool)

        self.mol_to_fp = mol_to_fp
        self._restored = True
        print(f"Restored variables from {weight_path}")

    def _load_vars(self, weight_path):
        """
        Load model variables (weights and biases) from a pickle file.

        Args:
            weight_path (str): Path to the model weights file.
        """
        with open(weight_path, 'rb') as f:
            # Use latin1 encoding to handle Python 2.x saved pickles
            self.vars = pickle.load(f, encoding='latin1')
        for i, var in enumerate(self.vars):
            print(f"Var {i}: shape = {var.shape}")  # Debug: Print shapes of variables

    def smi_to_fp(self, smi):
        """
        Convert a SMILES string to a fingerprint.

        Args:
            smi (str): SMILES string.

        Returns:
            np.ndarray: Fingerprint array.
        """
        if not smi:
            return np.zeros((FP_LEN,), dtype=bool)
        mol = Chem.MolFromSmiles(smi)
        fp = self.mol_to_fp(mol)
        return fp

    def apply(self, x):
        """
        Apply the model to an input fingerprint.

        Args:
            x (np.ndarray): Input fingerprint.

        Returns:
            float: Predicted score.
        """
        if not self._restored:
            raise ValueError("Model weights must be restored before applying.")

        for i in range(0, len(self.vars), 2):
            W = self.vars[i]
            b = self.vars[i + 1]
            if x.shape[0] != W.shape[0]:  # Check if input dimension matches weight dimension
                raise ValueError(f"Input dimension {x.shape[0]} does not match weight dimension {W.shape[0]}")
            x = np.matmul(x, W) + b
            if i < len(self.vars) - 2:
                x = np.maximum(x, 0)  # ReLU activation
        return 1 + (self.score_scale - 1) * (1 / (1 + np.exp(-x)))

    def get_score_from_smi(self, smi):
        """
        Predict the synthetic complexity score for a given SMILES string.

        Args:
            smi (str): SMILES string.

        Returns:
            float: Synthetic complexity score.
        """
        fp = self.smi_to_fp(smi)
        score = self.apply(fp)
        return float(score.item())  # Convert NumPy array to Python float