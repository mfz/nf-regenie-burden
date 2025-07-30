#!/bin/bash

input="$1"  # Replace with your actual file name
first="$2"
last="$3"

# Read the header and split it into an array
IFS=$'\t' read -r -a columns < <(head -n 1 "$input")

# Loop over phenotype columns (columns 3 to end, index 2+)
for ((i=$(($first + 1)); i<=$(($last + 1)); i++)); do
            pheno="${columns[$i]}"
                col_num=$((i + 1))  # 1-based indexing for cut
                    cut -f1,2,$col_num "$input" > "${pheno}.tsv"
            done
