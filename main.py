"""
Main entry point for SCScorer.
"""

from scscore.model import SCScorer

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="SCScorer: Predict synthetic complexity score for a given SMILES.")
    parser.add_argument("smiles", type=str, help="Input SMILES string.")
    args = parser.parse_args()

    scorer = SCScorer()
    scorer.restore()
    score = scorer.get_score_from_smi(args.smiles)
    print(f"SCScore for '{args.smiles}': {score}")