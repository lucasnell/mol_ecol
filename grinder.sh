#!/usr/bin/env bash

# Runs grinder simulator on fasta file(s)

cd "~/Desktop/gd/fasta"
tot_reads=100
out_name=gr_out
out_dir=./_grinder_out
references="_ApeKI_s0102.fa"



# grinder "$references" -total_reads $tot_reads -base_name $out_name \
grinder -rf $references -total_reads $tot_reads -base_name $out_name \
    -output_dir $out_dir \
    -length_bias 0 -fastq_output 1 -qual_levels 30 10 -insert_dist 200
    # -length_bias 0 -fastq_output 1 -qual_levels 30 10 -unidirectional 1

# -length_bias 0  --> Larger fragments contribute the same amount as small ones
# -fastq_output 1 --> Generate FASTQ output
# -qual_levels 30 10 --> Specifies the good and bad quality scores, respectively
# -unidirectional 1 --> Proceed unidirectionally, from forward strand only
# -md <mutation_dist>... | -mutation_dist <mutation_dist>...
#         Introduce sequencing errors in the reads, under the form of
#         mutations (substitutions, insertions and deletions) at positions
#         that follow a specified distribution (with replacement): model
#         (uniform, linear, poly4), model parameters. For example, for a
#         uniform 0.1% error rate, use: uniform 0.1. To simulate Sanger
#         errors, use a linear model where the errror rate is 1% at the 5' end
#         of reads and 2% at the 3' end: linear 1 2. To model Illumina errors
#         using the 4th degree polynome 3e-3 + 3.3e-8 * i^4 (Korbel et al
#         2009), use: poly4 3e-3 3.3e-8. Use the <mutation_ratio> option to
#         alter how many of these mutations are substitutions or indels.
#         Default: uniform 0 0
# -mr <mutation_ratio>... | -mutation_ratio <mutation_ratio>...
#         Indicate the percentage of substitutions and the number of indels
#         (insertions and deletions). For example, use '80 20' (4
#         substitutions for each indel) for Sanger reads. Note that this
#         parameter has no effect unless you specify the <mutation_dist>
#         option. Default: 80 20
# -id <insert_dist>... | -insert_dist <insert_dist>...
#         Create paired-end or mate-pair reads spanning the given insert
#         length. Important: the insert is defined in the biological sense,
#         i.e. its length includes the length of both reads and of the stretch
#         of DNA between them: 0 : off, or: insert size distribution in bp, in
#         the same format as the read length distribution (a typical value is
#         2,500 bp for mate pairs) Two distinct reads are generated whether or
#         not the mate pair overlaps. Default: 0



