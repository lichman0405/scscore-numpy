"""
Test Cases for SCScorer

This script provides a set of predefined SMILES strings to test the SCScorer model.
Each test case includes a SMILES string and its expected behavior.

Usage:
Run this script to validate the SCScorer functionality with various test cases.
"""

from scscore.model import SCScorer

def run_test_cases():
    # Initialize SCScorer and restore model
    scorer = SCScorer()
    scorer.restore()

    # Define test cases
    test_cases = [
        {"smiles": "CCO", "description": "Simple ethanol molecule"},
        {"smiles": "C1=CC=CC=C1", "description": "Benzene ring"},
        {"smiles": "CC(C)C(=O)O", "description": "Simple carboxylic acid"},
        {"smiles": "C[C@@H](O)[C@H](O)CO", "description": "Glyceraldehyde molecule"},
        {"smiles": "", "description": "Empty SMILES string"},
        {"smiles": "invalid_smiles", "description": "Invalid SMILES string"},
    ]

    # Run each test case
    for i, test in enumerate(test_cases):
        smiles = test["smiles"]
        description = test["description"]

        try:
            score = scorer.get_score_from_smi(smiles)
            print(f"Test Case {i + 1}: {description}")
            print(f"  SMILES: {smiles}")
            print(f"  Synthetic Complexity Score: {score}\n")
        except Exception as e:
            print(f"Test Case {i + 1}: {description}")
            print(f"  SMILES: {smiles}")
            print(f"  Error: {e}\n")

if __name__ == "__main__":
    run_test_cases()