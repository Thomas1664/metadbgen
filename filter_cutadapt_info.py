import pandas as pd
import sys

if len(sys.argv) != 2:
    print("Usage: python filter_cutadapt_info.py <input_file>")
    sys.exit(1)

input_path = sys.argv[1]

df = pd.read_csv(input_path, sep="\t", header=None)
df = df[df.iloc[:, 1] != -1]
df.to_csv(input_path, sep="\t", header=False, index=False)
