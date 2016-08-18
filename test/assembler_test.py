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
        asl.create_tables(self.conn)
        self.db = asl.SqliteDB(self.conn)

    def tearDown(self):
        if os.path.exists(TEST_DATABASE):
            os.remove(TEST_DATABASE)

    def test_insert_cols(self):
        ids = asl.db_insert_cols(self.conn, ['col1', 'col2', 'col3'])
        self.assertEquals(len(ids), 3)
        cur = self.conn.cursor()
        cur.execute('select count(*) from columns')
        self.assertEquals(3, cur.fetchone()[0])

    def test_insert_rows(self):
        ids = asl.db_insert_rows(self.conn, ['row1', 'row2', 'row3'])
        self.assertEquals(len(ids), 3)
        cur = self.conn.cursor()
        cur.execute('select count(*) from rows')
        self.assertEquals(3, cur.fetchone()[0])


if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(SqliteDBTest)]

    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))
