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


if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(ResampleTest)]

    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))
