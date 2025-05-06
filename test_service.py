import requests

# Define the URL of the FastAPI service
url = "http://127.0.0.1:8000/calculate_score"

# Test cases
test_cases = [
    {"smiles": "CCO", "description": "Simple ethanol molecule"},
    {"smiles": "C1=CC=CC=C1", "description": "Benzene ring"},
    {"smiles": "CC(C)C(=O)O", "description": "Simple carboxylic acid"},
    {"smiles": "C[C@@H](O)[C@H](O)CO", "description": "Glyceraldehyde molecule"},
    {"smiles": "CCCCCCCCCCCCCCCCCCCCCC", "description": "Long alkane chain"},
    {"smiles": "C1CCCC1", "description": "Cyclopentane"},
    {"smiles": "C1=CC=CC=C1C(CO)O", "description": "Phenyl alcohol"},
    {"smiles": "CC(C)(C)C(=O)O", "description": "Tert-butyl carboxylic acid"},
    {"smiles": "", "description": "Empty SMILES string"},
    {"smiles": "invalid_smiles", "description": "Invalid SMILES string"},
    {"smiles": "CCN(CC)CC", "description": "Simple amine"},
    {"smiles": "C1=CC=C2C=CC=CC2=C1", "description": "Naphthalene"},
    {"smiles": "C1=CC=C(C=C1)C2=CC=CC=C2", "description": "Biphenyl"},
    {"smiles": "C1CCCCC1C(=O)O", "description": "Cyclohexane carboxylic acid"},
    {"smiles": "C(C(C(C(C(=O)O)O)O)O)O", "description": "Simple sugar (glucose-like)"},
]

# Run each test case
for i, test in enumerate(test_cases):
    smiles = test["smiles"]
    description = test["description"]

    print(f"Running Test Case {i + 1}: {description}")
    print(f"  Input SMILES: {smiles}")
    
    try:
        # Send POST request
        response = requests.post(url, json={"smiles": smiles})
        
        # Print response
        if response.status_code == 200:
            print(f"  Response: {response.json()}\n")
        else:
            print(f"  Error: {response.json()['detail']}\n")
    except Exception as e:
        print(f"  Exception occurred: {e}\n")