

# Align fastq files to aphid genome


cd "/Volumes/750gb"

cd ./ref

bwa index aphid_genome.fa
samtools faidx aphid_genome.fa

cd "/Volumes/750gb/fastq"

cores=3
name=apeki


# echo $(date +%H:%M:%S) > ~/Desktop/times.fa
# bwa mem -t $cores ../ref/aphid_genome.fa ${name}.fq | \
#     samtools view -bh -@ $(expr ${cores} - 1) - > ../bam/${name}.bam
# echo $(date +%H:%M:%S) >> ~/Desktop/times.fa
# # Took ~4 hrs with 100x coverage of all samples


for i in {1..10}
do
    echo -e "Started $i\t\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa
    bwa mem -t $cores -v 2 ../ref/aphid_genome.fa ${name}_${i}.fq | \
        samtools view -bh -@ $(expr ${cores} - 1) - > ../bam/${name}_${i}.bam
    echo -e "Finished aligning\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa
    cd ../bam
    samtools sort -o ${name}_${i}_sorted.bam -T ${name}_${i}_s -@ ${cores} ${name}_${i}.bam
    echo -e "Finished sorting\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa
    samtools view -q 20 -b ${name}_${i}_sorted.bam -@ $(expr ${cores} - 1) > ${name}_${i}_hi.bam
    echo -e "Finished filtering\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa
    mv ${name}_${i}_sorted.bam ${name}_${i}.bam
    samtools index ${name}_${i}_hi.bam
    echo -e "Finished indexing\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa
    echo -e "\n\n" >> ~/Desktop/times.fa
    cd ../fastq
done







# cd ../bam
# # Sort reads (this took ~1 hr with 100x coverage):
# samtools sort -o ${name}_sorted.bam -T ${name}_s -@ ${cores} ${name}.bam
# # Keep only high-quality reads:
# samtools view -q 20 -b ${name}_sorted.bam -@ $(expr ${cores} - 1) > ${name}_hi.bam
# # I only use sorted BAM files
# mv ${name}_sorted.bam ${name}.bam
# # Index bc sometimes this is useful
# samtools index ${name}_hi.bam


bams="${name}_hi.bam"
for i in {1..10}
do
    bams+=" ${name}_${i}_hi.bam"
done

# Now making mpileup file to input to popoolation2 (took 30 min for 100x coverage)
samtools mpileup -f ../ref/aphid_genome.fa \
    -d 10000 \
    -B ${bams} > ../allele_freq/${name}.mpileup
    # -B ${name}_hi.bam > ${name}_small.mpileup
    # -r GL349621.1:17500 \

# Creating synchronized files
# According to the popoolation2 documentation 
# (https://sourceforge.net/p/popoolation2/wiki/Tutorial/)...
#       "Synchronized files... basically contain the allele frequencies for every 
#       population at every base in the reference genome in a concise format. 
#       Note that the synchronized file format contains the allele frequencies after 
#       filtering for base quality."

# 
# echo -e "mpileup2sync.jar started\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa
# java -ea -Xmx4g -jar ../popoolation2/mpileup2sync.jar --input ${name}.mpileup \
#     --output ${name}_lo.sync --fastq-type illumina --min-qual 0 --threads ${cores}
# echo -e "mpileup2sync.jar ended\t\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa


# Takes ~3 hrs
echo -e "mpileup2sync.pl started\t\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa
perl ../popoolation2/mpileup2sync.pl --input ${name}.mpileup \
    --output ${name}.sync --fastq-type illumina --min-qual 1
echo -e "mpileup2sync.pl ended\t\t" $(date +%H:%M:%S) >> ~/Desktop/times.fa



