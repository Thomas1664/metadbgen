import sys

if len(sys.argv) != 2:
    print("Usage: python filter_cutadapt_info.py <input_file>")
    sys.exit(1)

input_path = sys.argv[1]

lines = ""

# This script filters out any sequences which were not cut by cutadapt.
# We can't use Pandas here because for sequences that were not cut
# have fewer columns which breaks Pandas.
with open(input_path, "r") as file:
    for line in file.readlines():
        # we only care about the first 2 columns
        split = line.split('\t', 2)
        # Avoid out-of-bounds access
        # When this happens, the file is most likely invalid.
        if len(split) < 3:
            break
        # Skip rows where none of the primers was cut
        if split[1] != '-1':
            lines += f'{line}\n'
# Write only rows for sequences where at least one of the primers was cut
with open(input_path, "w") as file:
    file.write(lines)
