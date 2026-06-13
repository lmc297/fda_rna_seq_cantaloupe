# Loop through each raw fastq.gz file
for f in ./raw_reads/*.fastq.gz
do
# Run FastQC
fastqc -t 1 $f
done
