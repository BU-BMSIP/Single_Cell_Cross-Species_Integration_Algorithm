#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <protein_fasta_file>"
    exit 1
fi

fasta_file="$1"

# Count the number of lines starting with ">"
seq_count=$(grep -c '^>' "$fasta_file")

echo "Number of protein sequences: $seq_count"
