# Author: Katherina Cortes
# Date: June 8, 2022
# Purpose: test code

import os, sys

cur_path=os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, cur_path+"/..")

import unittest
from src import plotGenerator, database, alignment, dataExploration, kmerBinning



class KmerBinningTests(unittest.TestCase):
    # test kmers from virus genomes are correctly made
    def testKmerCreation(self):
        data = {'1': 'thecats', '2': 'thecatty'}
        kmers = database.getkmers(data, 5)
        self.assertEqual(kmers, {'theca': ['1', '2'], 'hecat': ['1', '2'], 'ecats': ['1'], 'ecatt': ['2'], 'catty': ['2']})


class AlignmentTest(unittest.TestCase):
    # test correct alignment
    def testGlobalA(self):
        s1 = 'cat'
        s2 = 'ccatb'
        score = alignment.globalAlign(s1,s2)
        self.assertEqual(score, 9)

    # test it doenst matter which sequence first or second
    def testGlobalRev(self):
        s1 = 'cat'
        s2 = 'ccatb'

        score1 = alignment.globalAlign(s1, s2)
        score2 = alignment.globalAlign(s2, s1)
        self.assertEqual(score1, score2)


if __name__ == '__main__':
    unittest.main()
