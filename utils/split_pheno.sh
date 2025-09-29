#!/bin/bash

# Usage: ./split_phenotype_file.sh phenotype.tsv P

set -e

INPUT_FILE="$1"
P="$2"

if [[ -z "$INPUT_FILE" || -z "$P" ]]; then
    echo "Usage: $0 <phenotype_file.tsv> <number_of_phenotypes_per_split>"
    exit 1
fi

# Read the header
HEADER=$(head -n 1 "$INPUT_FILE")
# Extract phenotype names (excluding FID and IID)
PHENOS=($(echo "$HEADER" | cut -f3-))
NUM_PHENOS=${#PHENOS[@]}

# Compute number of chunks
NUM_CHUNKS=$(( (NUM_PHENOS + P - 1) / P ))

for ((i=0; i<NUM_CHUNKS; i++)); do
    START=$((i * P + 3))  # +3 to skip FID and IID
    END=$((START + P - 1))
    OUTFILE="phenotype_split_$i.tsv"
    cut -f1,2,"$START"-"$END" "$INPUT_FILE" > "$OUTFILE"
done
