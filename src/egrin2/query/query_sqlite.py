import sqlite3

def genes_with_name(conn, name):
    cursor = conn.cursor()
    cursor.execute('select rowid from rows where name=?', (name, ))
    results = cursor.fetchall()
    if len(results) > 0:
        return results
    else:
        return None
