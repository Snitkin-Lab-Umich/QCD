#!/bin/bash

# Arguments passed to the script
#bioawk_params=$1
input_spades_assembly=$1
output_spades_l1000_assembly=$2
out_dir=$3
prefix=$4


#set -x  # Enable script debugging

# Evaluate bioawk_params as a command
contig_count=$(bioawk -c fastx '{if(length($seq) > 1000) {print ">"$name; print $seq}}' $input_spades_assembly | grep -c '>')
#echo $bioawk_params
echo "Contig count: $contig_count"
#contig_count=$(bioawk -c fastx '{if(length($seq) > 1000) {print ">"$name; print $seq}}' results/2024-06-21_Project_Marimba_11091-ES_QCD/11091-ES-12_S137/spades/contigs.fasta | grep -c '>')

#this works
#./bioawk.sh results/2024-06-21_Project_Marimba_11091-ES_QCD/11091-ES-12_S137/spades/contigs.fasta results/2024-06-21_Project_Marimba_11091-ES_QCD/11091-ES-12_S137/spades/11091-ES-12_S137_contigs_l1000.fasta results/2024-06-21_Project_Marimba_11091-ES_QCD/11091-ES-12_S137/spades/ 11091-ES-12_S137

if [ "$contig_count" -eq 0 ]; then
    echo "Copying file $input_spades_assembly to $output_spades_l1000_assembly"
    cp "$input_spades_assembly" "$output_spades_l1000_assembly" &&
    grep '>' "$output_spades_l1000_assembly" > "$out_dir/spades_assembly_header_info.txt" &&
    sed -i "s/>NODE_/>${prefix}_/g" "$output_spades_l1000_assembly" &&
    sed -i 's/_length_[0-9]\+_cov_[0-9.]*//g' "$output_spades_l1000_assembly"
else
    echo "Filtering file"
    bioawk -c fastx '{if(length($seq) > 1000) {print ">"$name; print $seq}}' "$input_spades_assembly" > "$output_spades_l1000_assembly" &&
    grep '>' "$output_spades_l1000_assembly" > "$out_dir/spades_assembly_header_info.txt" &&
    sed -i "s/>NODE_/>${prefix}_/g" "$output_spades_l1000_assembly" &&
    sed -i 's/_length_[0-9]\+_cov_[0-9.]*//g' "$output_spades_l1000_assembly"
fi
