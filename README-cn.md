# SCScore-Numpy

SCScore-Numpy 是一个基于 NumPy 的分子合成复杂性评分模型。本项目实现了一个独立的 Python 模型，用于加载预训练的权重并计算给定分子（SMILES 表示）的合成复杂性分数。此外，项目还提供了一个基于 **FastAPI** 的 Web 服务，支持通过 HTTP 接口计算分数。

---

## 功能概览
1. **分子合成复杂性评分**：使用预训练模型计算分子的合成复杂性分数。
2. **支持 SMILES 输入**：接收标准的 SMILES 字符串进行评分。
3. **FastAPI 服务**：通过 RESTful API 提供分数计算功能。
4. **测试脚本**：包含单元测试和服务测试用例。

---

## 快速开始

### 环境准备
1. 克隆项目：
   ```bash
   git clone https://github.com/your-repo/scscore-numpy.git
   cd scscore-numpy
   ```

2. 创建虚拟环境并安装依赖：
   ```bash
   conda create -n scscore-numpy python=3.8 -y
   conda activate scscore-numpy
   pip install -r requirements.txt
   ```

3. 下载模型权重：
   请确保模型文件 `model.ckpt-10654.as_numpy.pickle` 已放置在 `scscore/models/` 目录下。

---

## 使用方法

### 1. 命令行运行
运行以下脚本，直接调用模型对单个 SMILES 字符串进行评分：
```bash
python main.py
```

确保 `main.py` 中已指定正确的 SMILES 字符串。

---

### 2. FastAPI 服务
#### 启动服务
运行以下命令启动 FastAPI 服务：
```bash
uvicorn fastapi_service:app --reload --host 0.0.0.0 --port 8000
```

服务启动后，可通过以下 URL 访问：
- Swagger 文档：`http://127.0.0.1:8000/docs`
- ReDoc 文档：`http://127.0.0.1:8000/redoc`

#### API 示例
**Endpoint**: `POST /calculate_score`

**请求示例**：
```json
{
  "smiles": "C1=CC=C2C=CC=CC2=C1"
}
```

**响应示例**：
```json
{
  "smiles": "C1=CC=C2C=CC=CC2=C1",
  "score": 1.5987794399261475
}
```

**PowerShell 测试命令**：
```powershell
curl -X POST http://127.0.0.1:8000/calculate_score -H "Content-Type: application/json" --data '{"smiles": "C1=CC=C2C=CC=CC2=C1"}'
```

---

## 测试

### 1. 单元测试
运行以下命令执行单元测试：
```bash
python -m unittest discover tests
```

### 2. 服务测试
文件 `test_service.py` 包含对 FastAPI 服务的测试用例，运行以下命令执行：
```bash
python test_service.py
```

#### 示例输出
```plaintext
Running Test Case 1: Simple ethanol molecule
  Input SMILES: CCO
  Response: {'smiles': 'CCO', 'score': 2.345678}

Running Test Case 2: Benzene ring
  Input SMILES: C1=CC=CC=C1
  Response: {'smiles': 'C1=CC=CC=C1', 'score': 3.567890}

Running Test Case 3: Simple carboxylic acid
  Input SMILES: CC(C)C(=O)O
  Response: {'smiles': 'CC(C)C(=O)O', 'score': 2.987654}
```

---

## 项目结构
```plaintext
scscore-numpy/
│
├── scscore/                     # 核心代码模块
│   ├── __init__.py              # 模块初始化
│   ├── model.py                 # SCScore NumPy 推理模型
│   ├── utils.py                 # 工具函数（可选）
│   └── models/                  # 模型文件存放目录
│       └── model.ckpt-10654.as_numpy.pickle # 模型文件
│
├── tests/                       # 测试代码目录
│   ├── __init__.py              # 使其成为模块
│   └── test_model.py            # 测试推理功能
│
├── main.py                      # 程序入口
├── fastapi_service.py           # FastAPI 服务实现
├── test_service.py              # 测试服务的脚本
├── requirements.txt             # 项目依赖
└── README.md                    # 项目说明文档
```

---

## 依赖
以下是项目主要依赖的 Python 包：
```plaintext
fastapi
uvicorn
numpy
rdkit
```

安装所有依赖：
```bash
pip install -r requirements.txt
```

---

## 常见问题
### 1. 如何处理无效的 SMILES 字符串？
FastAPI 服务会返回 400 错误码，并提示如下：
```json
{
  "detail": "SMILES string cannot be empty."
}
```

### 2. 模型文件丢失怎么办？
请确保 `model.ckpt-10654.as_numpy.pickle` 文件存在于 `scscore/models/` 目录下。

---

## 贡献
欢迎提交 Issue 或 Pull Request 来贡献代码！在提交之前，请确保通过所有单元测试。

---

## 许可证
本项目基于 MIT 协议开源。