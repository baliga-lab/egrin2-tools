#!/usr/bin/env python3
import sys
import unittest
import xmlrunner
import mongomock
import pandas
import numpy as np

import query.egrin2_query as e2q


class QueryTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test suite for the egrin2_query module"""

    def test_remove_list_duplicates(self):
        self.assertEquals([1, 2, 3], e2q.remove_list_duplicates([1, 2, 3]))
        self.assertEquals([1, 2, 3], e2q.remove_list_duplicates([1, 2, 3, 1, 2, 3]))

    def test_rsd(self):
        """don't care for empty list"""
        #self.assertTrue(np.isnan(e2q.rsd([])))
        self.assertEquals(0, e2q.rsd([1]))
        self.assertEquals(1/3, e2q.rsd([1, 2]))

    def test_find_match(self):
        df = pandas.DataFrame([[42, 621, 53], [14, 244, 621]], index=['row 1', 'row 2'], columns=['col 1', 'col 2', 'col 3'])
        self.assertEquals(42, e2q.find_match(42, df, 'col 1'))
        self.assertEquals(53, e2q.find_match(42, df, 'col 3'))
        self.assertEquals(244, e2q.find_match(14, df, 'col 2'))
        self.assertEquals(53, e2q.find_match(621, df, 'col 3'))
        self.assertEquals(None, e2q.find_match(4711, df, 'col 3'))


if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(QueryTest)]

    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))

