import sqlite3

def genes_with_name(db, name):
    cursor = db.cursor()
