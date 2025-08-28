import sys
from src import database


msa, tax, db = sys.argv[1:]
database.export_sintax(msa, tax, db)
