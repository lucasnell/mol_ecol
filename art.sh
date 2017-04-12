

# Runs ART simulator on fasta file(s)

cd "/Volumes/64gb/fasta"
rcount=100
out_name=../fastq/art_test/art_out
prof=../art_profile/HiSeq2kL100R1.txt
reference="ApeKI.fa"


# -1 and -2 arguments make it use the same profile for both R1 and R2
# -nf 0 removes "N"-masking
# -ir2 and -dr2 make R2 insertion and deletiion rates equal those for R1
# -rs sets random number generator



art_illumina -i $reference -l 100 -c $rcount -m 200 -s 0 -p -o ${out_name} \
    -1 $prof -2 $prof -na \
    -ir2 0.00009 -dr2 0.00011 -nf 0 \
    -rs 1492027404
# The above took 1 hr, 15 min


# This takes a while, too
cd ${out_name/art_out/}
for f in *.fq
do
    gzip $f
done



# ==========================================
# ==========================================

# Some of ART's arguments

# ==========================================
# ==========================================



# -1   --qprof1   the first-read quality profile
# -2   --qprof2   the second-read quality profile
# -f --fcov       the fold of read coverage to be simulated or number of reads/read pairs 
#                     generated for each amplicon
# -c   --rcount   number of reads/read pairs to be generated per sequence/amplicon 
#                     (not be used together with -f/--fcov)
# -i   --in       the filename of input DNA/RNA reference
# -l   --len      the length of reads to be simulated
# -o   --out      the prefix of output filename
# -sam --samout   indicate to generate SAM alignment file
# -m   --mflen    the mean size of DNA/RNA fragments for paired-end simulations
# -s   --sdev     the standard deviation of DNA/RNA fragment size for paired-end simulations.
# -nf  --maskN    the cutoff frequency of 'N' in a window size of the read length for 
#                     masking genomic regions 
#                     NOTE: default: '-nf 1' to mask all regions with 'N'. Use '-nf 0' 
#                       to turn off masking
# -na  --noALN    do not output ALN alignment file
# -p   --paired   indicate a paired-end read simulation or to generate reads from both 
#                     ends of amplicons
# -ss  --seqSys   The name of Illumina sequencing system of the built-in profile used for 
#                 simulation
# -ir  --insRate  the first-read insertion rate (default: 0.00009)
# -ir2 --insRate2 the second-read insertion rate (default: 0.00015)
# -dr  --delRate  the first-read deletion rate (default:  0.00011)
# -dr2 --delRate2 the second-read deletion rate (default: 0.00023)
# -rs  --rndSeed  the seed for random number generator (default: system time in second)
#                     NOTE: using a fixed seed to generate two identical datasets from 
#                     different runs
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




