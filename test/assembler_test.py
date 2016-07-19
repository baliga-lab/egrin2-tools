#!/usr/bin/env python3
import sys
import unittest
import xmlrunner
import sqlite3
import os

import assemble.assemble_sqlite as asl

TEST_DATABASE = 'assemble_sqlite_test.db'

class SqliteDBTest(unittest.TestCase):

    def setUp(self):
        if os.path.exists(TEST_DATABASE):
            os.remove(TEST_DATABASE)
        self.conn = sqlite3.connect(TEST_DATABASE)
        self.db = asl.SqliteDB(self.conn)

    def test(self):
        pass

if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(SqliteDBTest)]

    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))
