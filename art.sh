#!/usr/bin/env bash

# Runs ART simulator on fasta file(s)

cd "~/Desktop/gd/fasta"
tot_reads=100
out_name=gr_out
out_dir=../grinder_out
digest="ApeKI"
references=""
for f in "${digest}"*.fa.gz
do
    references+="-reference_file $f "
done




# -f --fcov  the fold of read coverage to be simulated or number of reads/read pairs 
#            generated for each amplicon
# -i   --in       the filename of input DNA/RNA reference
# -l   --len      the length of reads to be simulated
# -o   --out      the prefix of output filename
# -sam --samout   indicate to generate SAM alignment file
# -ss  --seqSys   The name of Illumina sequencing system of the built-in profile used for simulation
# NOTE: all built-in sequencing system ID names are:
#     GA1 - GenomeAnalyzer I (36bp,44bp)
#     GA2 - GenomeAnalyzer II (50bp, 75bp) 
#     HS10 - HiSeq 1000 (100bp)
#     HS20 - HiSeq 2000 (100bp)             <- ***
#     HS25 - HiSeq 2500 (125bp, 150bp)      <- ***
#     HSXn - HiSeqX PCR free (150bp)
#     HSXt - HiSeqX TruSeq (150bp)
#     MinS - MiniSeq TruSeq (50bp)
#     MSv1 - MiSeq v1 (250bp)
#     MSv3 - MiSeq v3 (250bp)
#     NS50 - NextSeq500 v2 (75bp)



