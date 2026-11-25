import pandas as pd
import csv
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Path to input phenotype matrix")
parser.add_argument("pheno_list", help="Path to file with phenotype names")
parser.add_argument("-b", "--batch-size", type=int, default=25, help="Number of phenotypes per batch")
parser.add_argument("-p", "--prefix", type=str, default="phenotypes", help="output prefix")

args = parser.parse_args()



input_file = args.input_file
batch = args.batch_size
prefix = args.prefix

with open(args.pheno_list, encoding="utf-8") as f:
    wanted_phenotypes = [line.rstrip("\r\n") for line in f if line.strip()]



columns = pd.read_csv(input_file, sep='\t', nrows=0).columns.tolist()
assert columns[:2] == ['FID', 'IID'], "Expect first 2 columns to be FID, IID"

missing = set(wanted_phenotypes) - set(columns)
if missing:
    raise ValueError(f"Columns not found: {sorted(missing)}")

indices = [columns.index(name) for name in wanted_phenotypes]

outputs = []
for part, start in enumerate(range(0, len(wanted_phenotypes), batch), start=1):
    group = ['FID', 'IID'] + wanted_phenotypes[start:start + batch]
    group_indices = [0, 1] + indices[start:start + batch]

    out = open(f"{prefix}_{part:03}.tsv", "w", newline="")
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerow(group)
    outputs.append((out, writer, group_indices))

try:
    with open(input_file, newline="") as src:
        reader = csv.reader(src, delimiter="\t")
        next(reader)
        for row in reader:
            for out, writer, group_indices in outputs:
                writer.writerow([row[i] for i in group_indices])

finally:
    for out, _, _ in outputs:
        out.close()
