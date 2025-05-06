from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from scscore.model import SCScorer

# Initialize FastAPI app
app = FastAPI()

# Initialize SCScorer model
scorer = SCScorer()
scorer.restore()

# Define request body schema
class SMILESRequest(BaseModel):
    smiles: str

@app.post("/calculate_score")
async def calculate_score(smiles_request: SMILESRequest):
    """
    Calculate the synthetic complexity score for a given SMILES string.
    """
    smiles = smiles_request.smiles
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES string cannot be empty.")

    try:
        # Calculate the score
        score = scorer.get_score_from_smi(smiles)
        return {
            "smiles": smiles,
            "score": score
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error calculating score: {str(e)}")