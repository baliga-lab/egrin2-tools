#!/usr/bin/env python3
import sys
import unittest
import xmlrunner
import numpy as np
import pandas as pd

import assemble.resample as resample

class ResampleTest(unittest.TestCase):  # pylint: disable-msg=R0904

    def test_rsd(self):
        a = np.array([1, 2, 3])
        self.assertAlmostEquals(0.408248290464, resample.rsd(a))

    def test_make_col_id_rsd_groups(self):
        df = pd.DataFrame(np.array([[1, 2.0, 2.12231], [2, 2.3, 2.31231], [3, 0.123123, 0.1212221],
                                    [4, -1.2313213, -1.212131231], [5, 4.123123123, 4.012313123],
                                    [1, -0.1231232, -0.1121121], [2, 3.1231231, 3.012313213],
                                    [3, -0.2231232, -0.2121121], [4, -1.1231231, -1.012313213],
                                    [5, -0.2231232, -0.2121121], [1, -1.1231231, -1.012313213],
                                    [2, -0.2231232, -0.2121121], [3, -1.1231231, -1.012313213]]),
                          index=['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11',
                                 'g12', 'g13'],
                          columns=['col_id', 'raw_expression', 'standardized_expression'])
        nrows = 2
        nresamples = 1
        gb = resample.make_col_id_rsd_groups(df, nrows, nresamples)
        col_ids = sorted(map(int, gb.groups.keys()))
        # there is currently not really a good way to test the results, but testing that we
        # have all groups is a start
        self.assertEquals(col_ids, [1, 2, 3, 4, 5])


if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(ResampleTest)]

    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))
