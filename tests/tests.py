# Author: Katherina Cortes
# Date: June 8, 2022
# Purpose: test code

import os, sys

cur_path=os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, cur_path+"/..")

import unittest
from src import plotGenerator, database, alignment, dataExploration


class DatabaseCreation(unittest.TestCase):
    def test(self):
        self.assert_()

class AlignmentTest(unittest.TestCase):
    def testGlobalA(self):
        s1 = 'cat'
        s2 = 'ccatb'
        alignment.globalAlign(s1,s2)
        self.assertEqual()

    def testGlobalRev(self):
        s1 = 'cat'
        s2 = 'ccatb'

        score1 = alignment.globalAlign(s1, s2)
        score2 = alignment.globalAlign(s2, s1)
        self.assertEqual(score1, score2)


if __name__ == '__main__':
    unittest.main()
