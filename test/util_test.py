#!/usr/bin/env python3
import sys
import unittest
import xmlrunner
import assemble.util as util

class UtilTest(unittest.TestCase):  # pylint: disable-msg=R0904

    def test_split_list(self):
        inlist = [1, 2, 3, 4, 5, 6]
        self.assertEquals([[1], [2], [3], [4], [5], [6]],
                          util.split_list(inlist, 6))
        self.assertEquals([[1, 2, 3], [4, 5, 6]],
                          util.split_list(inlist, 2))
        self.assertEquals([[1, 2], [3, 4], [5, 6]],
                          util.split_list(inlist, 3))
        self.assertEquals([[1], [2, 3], [4], [5, 6]],
                          util.split_list(inlist, 4))

if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(UtilTest)]

    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))
