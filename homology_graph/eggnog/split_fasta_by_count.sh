#!/bin/bash

# Check input
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_fasta_file>"
    exit 1
fi

input="$1"
max_seqs=100000
count=0
file_index=1
out_file=$(printf "split_%03d.fasta" "$file_index")

# Create or empty the first output file
> "$out_file"

while read -r line; do
    if [[ $line == ">"* ]]; then
        ((count++))
        if ((count > max_seqs)); then
            ((file_index++))
            out_file=$(printf "split_%03d.fasta" "$file_index")
            > "$out_file"
            count=1
        fi
    fi
    echo "$line" >> "$out_file"
done < "$input"

echo "Done. Generated $file_index file(s)."
