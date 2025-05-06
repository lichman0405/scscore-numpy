# SCScore-Numpy

`SCScore-Numpy` is a molecular synthetic complexity scoring model implemented using NumPy. This project provides a standalone Python-based model to calculate synthetic complexity scores for molecules represented by SMILES strings. Additionally, it offers a **FastAPI** service for convenient HTTP-based scoring.

[中文文档](README-cn.md)

---

## Features
- **Synthetic Complexity Scoring**: Calculate the synthetic complexity score for a molecule using a pre-trained model.
- **SMILES Input**: Accept standard SMILES strings as input for scoring.
- **FastAPI Service**: A RESTful API for scoring molecules via HTTP.
- **Testing Scripts**: Includes unit tests and service test cases to ensure reliability.

---

## Quick Start

### Prerequisites
1. **Clone the repository**:
   ```bash
   git clone https://github.com/lichman0405/scscore-numpy.git
   cd scscore-numpy
   ```

2. **Set up the environment**:
   Create a Python virtual environment and install the required dependencies:
   ```bash
   conda create -n scscore-numpy python=3.8 -y
   conda activate scscore-numpy
   pip install -r requirements.txt
   ```

3. **Download Model Weights**:
   Place the model file `model.ckpt-10654.as_numpy.pickle` into the `scscore/models/` directory.

---

## Usage

### Command Line
Run the following script to compute the score for a single SMILES string:
```bash
python main.py
```

Ensure `main.py` specifies the correct SMILES string for scoring.

---

### FastAPI Service
Run the following command to start the FastAPI service:
```bash
uvicorn fastapi_service:app --reload --host 0.0.0.0 --port 8000
```

After starting, you can interact with the API using the following endpoints:
- **Swagger UI**: `http://127.0.0.1:8000/docs`
- **ReDoc Documentation**: `http://127.0.0.1:8000/redoc`

#### Example API Request
**URL**: `POST /calculate_score`

**Request Body**:
```json
{
  "smiles": "C1=CC=C2C=CC=CC2=C1"
}
```

**Response**:
```json
{
  "smiles": "C1=CC=C2C=CC=CC2=C1",
  "score": 1.5987794399261475
}
```

#### PowerShell Test
Use the following `curl` command to test the API:
```powershell
curl -X POST http://127.0.0.1:8000/calculate_score -H "Content-Type: application/json" -d '{"smiles": "C1=CC=C2C=CC=CC2=C1"}'
```

---

## Testing

### Unit Tests
Run the following command to execute unit tests:
```bash
python -m unittest discover tests
```

### Service Tests
The `test_service.py` file contains test cases for the FastAPI service. Run the following command to execute:
```bash
python test_service.py
```

---

## Project Structure
```plaintext
scscore-numpy/
│
├── scscore/                     # Core code module
│   ├── __init__.py              # Module initialization
│   ├── model.py                 # SCScore NumPy inference model
│   ├── utils.py                 # Utility functions
│   └── models/                  # Directory for model files
│       └── model.ckpt-10654.as_numpy.pickle # Model file
│
├── tests/                       # Directory for test scripts
│   ├── __init__.py              # Makes it a module
│   └── test_model.py            # Tests for inference functionality
│
├── main.py                      # Entry point script
├── fastapi_service.py           # FastAPI service implementation
├── test_service.py              # Script for testing the service
├── requirements.txt             # Project dependencies
├── README-en.md                 # English documentation
└── README-cn.md                 # Chinese documentation
```

---

## Dependencies
The project relies on the following Python libraries:
- **FastAPI**
- **Uvicorn**
- **NumPy**
- **RDKit**

Install all dependencies via:
```bash
pip install -r requirements.txt
```

---

## FAQ

### 1. How are invalid SMILES strings handled?
The FastAPI service will return a 400 error code with a message such as:
```json
{
  "detail": "SMILES string cannot be empty."
}
```

### 2. What should I do if the model file is missing?
Ensure that the file `model.ckpt-10654.as_numpy.pickle` exists in the `scscore/models/` directory.

---

## Contributing
We welcome contributions! Feel free to submit issues or pull requests. Before contributing, please ensure all tests pass.

---

## License
This project is licensed under the MIT License.