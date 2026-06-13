#!/bin/bash

# Loop through each set of trimmed reads
for f in trimmed_reads_stranded/*_trimmed_1_val_1.fq.gz
do
# Map to the LT2 reference genome with bwa mem
# Output a sorted bam file using samtools
bwa mem -t 4 reference_genome_lt2/salmonella_GCF_000006945.2.fna $f ${f%_trimmed_1_val_1.fq.gz}_trimmed_2_val_2.fq.gz | samtools sort -@ 4 -O bam - > ${f%_trimmed_1_val_1.fq.gz}_lt2.bam
# Index the sorted bam with samtools
samtools index -@ 4 ${f%_trimmed_1_val_1.fq.gz}_lt2.bam
# Output the number of reads mapped to the LT2 genome
samtools view -F 0x4 ${f%_trimmed_1_val_1.fq.gz}_lt2.bam | cut -f 1 | sort | uniq | wc -l > ${f%_trimmed_1_val_1.fq.gz}_lt2.counts
done
