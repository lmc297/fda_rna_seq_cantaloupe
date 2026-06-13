for f in ./raw_reads/*.fastq.gz
do
fastqc -t 1 $f
done
