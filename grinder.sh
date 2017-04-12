#!/usr/bin/env bash

# Runs grinder simulator on fasta file(s)

cd "/Volumes/64gb/fasta"
# tot_reads=62500000  # Quarter lane
tot_reads=2500000  # 1/100 lane
out_name=gr_out
out_dir="../fastq/test"
reference="ApeKI.fa.gz"
abund=../abundances/ApeKI_one.txt  # Just one library


echo -e "Start:\n"`date +%H:%M:%S` >> ~/Desktop/grinder_times
gunzip -c $reference | \
    grinder -rf - -total_reads $tot_reads -base_name $out_name \
    -output_dir $out_dir \
    -length_bias 0 -fastq_output 1 -qual_levels 30 10 -insert_dist 200 \
    # -abundance_file $abund \
    -mutation_dist poly4 3e-3 3.3e-8  # (Korbel et al 2009)
    # -homopolymer_dist "N(n, 0.03494 + n * 0.06856)" \  # (Balzer et al. 2010)
echo -e "End:\n"`date +%H:%M:%S` >> ~/Desktop/grinder_times



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
# -hd <homopolymer_dist> | -homopolymer_dist <homopolymer_dist>
#         Introduce sequencing errors in the reads under the form of
#         homopolymeric stretches (e.g. AAA, CCCCC) using a specified model
#         where the homopolymer length follows a normal distribution N(mean,
#         standard deviation) that is function of the homopolymer length n:
# 
#           Margulies: N(n, 0.15 * n)              ,  Margulies et al. 2005.
#           Richter  : N(n, 0.15 * sqrt(n))        ,  Richter et al. 2008.
#           Balzer   : N(n, 0.03494 + n * 0.06856) ,  Balzer et al. 2010.
# 
#         Default: 0
# 
# -id <insert_dist>... | -insert_dist <insert_dist>...
#         Create paired-end or mate-pair reads spanning the given insert
#         length. Important: the insert is defined in the biological sense,
#         i.e. its length includes the length of both reads and of the stretch
#         of DNA between them: 0 : off, or: insert size distribution in bp, in
#         the same format as the read length distribution (a typical value is
#         2,500 bp for mate pairs) Two distinct reads are generated whether or
#         not the mate pair overlaps. Default: 0
# -af <abundance_file> | -abundance_file <abundance_file>
#         Specify the relative abundance of the reference sequences manually
#         in an input file. Each line of the file should contain a sequence
#         name and its relative abundance (%), e.g. 'seqABC 82.1' or 'seqABC
#         82.1 10.2' if you are specifying two different libraries.
# -nl <num_libraries> | -num_libraries <num_libraries>
#         Number of independent libraries to create. Specify how diverse and
#         similar they should be with <diversity>, <shared_perc> and
#         <permuted_perc>. Assign them different MID tags with
#         <multiplex_mids>. Default: 1
# -rs <random_seed> | -random_seed <random_seed>
#         Seed number to use for the pseudo-random number generator.

