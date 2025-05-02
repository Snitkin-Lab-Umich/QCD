#!/bin/bash

# Arguments passed to the script
input_spades_assembly=$1
out_dir=$2

# Evaluate if there are any contigs greater than 1000
contig_count=$(awk '/^>/ { if (seq && length(seq) > 1000) { count++ } seq="" } !/^>/ { seq=seq $0 } END { if (seq && length(seq) > 1000) { count++ } print count }' "$input_spades_assembly")

# Check if contig_count is empty or zero
if [ -z "$contig_count" ] || [ "$contig_count" -eq 0 ]; then
    echo "Contig sizes are all ≤ 1000. Running QUAST with --min-contig 100."
    quast.py "$input_spades_assembly" -o "$out_dir" --min-contig 100 
else
    echo "Contig sizes ≥ 1000. Running QUAST with --contig-thresholds 0,1000,5000,10000,25000,50000."
    quast.py "$input_spades_assembly" -o "$out_dir" --contig-thresholds 0,1000,5000,10000,25000,50000
fi

