"""
Unit tests for SCScorer.
"""

import unittest
import os
import sys

# 添加项目根目录到模块搜索路径
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from scscore.model import SCScorer

class TestSCScorer(unittest.TestCase):
    def setUp(self):
        self.scorer = SCScorer()
        self.scorer.restore()

    def test_score(self):
        smiles = "CCO"
        score = self.scorer.get_score_from_smi(smiles)
        self.assertIsInstance(score, float)
        self.assertGreater(score, 0)

if __name__ == "__main__":
    unittest.main()