from src import taxonomy
import sys

if __name__ == '__main__':
    input_xml, out_tax, filtered_tax = sys.argv[1:]
    taxonomy.parse_taxonomy(input_xml, out_tax)
    taxonomy.filter_taxonomy(out_tax, filtered_tax)
