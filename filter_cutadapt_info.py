import contextlib
import sys


input_path, output_path, *rest = sys.argv[1:]
untrimmed_path = rest[0] if rest else None

# This script filters out any sequences which were not cut by cutadapt.
# We can't use Pandas here because sequences that were not cut
# have fewer columns which breaks Pandas.
with open(input_path, "r") as in_file, open(output_path, "w") as out_file, contextlib.ExitStack() as stack:
    untrimmed_file = stack.enter_context(open(untrimmed_path, "w")) if untrimmed_path else None
    for line in in_file:
        # we only care about the first 2 columns
        split = line.split('\t', 2)
        # Avoid out-of-bounds access
        # When this happens, the file is most likely invalid.
        if len(split) < 3:
            break
        # Skip rows where none of the primers was cut
        if split[1] != '-1':
            # Write only rows for sequences where at least one of the primers was cut
            out_file.write(f'{line}\n')
        elif untrimmed_file is not None:
            untrimmed_file.write(line)
