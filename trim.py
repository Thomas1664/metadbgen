from src import database, primer
import sys


in_meta, fasta, out_meta, cutadapt_info, cutadapt_out = sys.argv[1:]
database.import_trim(in_meta, fasta, out_meta)
primer.parse_cutadapt_info(cutadapt_info, cutadapt_out)
