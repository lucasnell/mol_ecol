#!/usr/bin/env bash

# Runs grinder simulator on fasta file(s)


tot_reads=100
out_name=gr_out
out_dir=~/Desktop/grinder_out
references=""
for f in *.fasta
do
    references+="-reference_file $f "
done



grinder "$references" -total_reads $tot_reads -base_name $out_name \
    -output_dir $out_dir \
    -length_bias 0 -fastq_output 1 -qual_levels 30 10 -unidirectional 1

# -length_bias 0  --> Larger fragments contribute the same amount as small ones
# -fastq_output 1 --> Generate FASTQ output
# -qual_levels 30 10 --> Specifies the good and bad quality scores, respectively
# -unidirectional 1 --> Proceed unidirectionally, from forward strand only
