

# Align fastq files to aphid genome


cd "/Volumes/750gb"

cd ./ref

bwa index aphid_genome.fa

# cd "/Volumes/750gb/fastq"

cd ~/Desktop/aln/fq

cores=3
name=apeki


echo $(date +%H:%M:%S) > ~/Desktop/times.fa
bwa mem -t $cores ../ref/aphid_genome.fa ${name}.fq | \
    samtools view -bh -@ $(expr ${cores} - 1) - > ../bam/${name}.bam
echo $(date +%H:%M:%S) >> ~/Desktop/times.fa
# Took ~4 hrs with 100x coverage


# echo $(date +%H:%M:%S) >> ~/Desktop/times.fa
samtools sort -o ${name}_sorted.bam -T ${name}_s -@ ${cores} ${name}.bam
# echo $(date +%H:%M:%S) >> ~/Desktop/times.fa
# Took ~1 hr with 100x coverage

mv ${name}_sorted.bam ${name}.bam

# Keep only high-quality reads
samtools view -q 20 -b ${name}.bam -@ $(expr ${cores} - 1) > ${name}_hi.bam
terminal-notifier -message DONE
