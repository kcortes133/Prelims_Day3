import os, sys

cur_path=os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, cur_path+"/..")

import unittest
from src import *


class tester(unittest.TestCase):
    def test(self):
        self.assert_()



if __name__ == '__main__':
    unittest.main()
