import csv
import sys
from glob import glob
from src import genbank_parser

gb_dir, out_file = sys.argv[1:]

files = glob(f'{gb_dir}/cox2_*.gb')

with open(out_file, 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    genbank_parser.dataclass_to_csv_header(csv_writer, genbank_parser.ParsedRecord)
    for file in files:
        genbank_parser.parse_genbank_file(file, csv_writer)

