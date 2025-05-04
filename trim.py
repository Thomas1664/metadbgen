from src import database
import sys


in_meta, fasta, out_meta = sys.argv[1:]
database.import_trim(in_meta, fasta, out_meta)
