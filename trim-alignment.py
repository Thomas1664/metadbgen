import sys
from src import primer


msa, cutadapt_info, meta, fasta_out = sys.argv[1:]
primer.trim_msa(cutadapt_info, meta, msa, fasta_out)
