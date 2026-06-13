#!/bin/bash

# Loop through bam files
# These include tons of canteloupe reads that don't map to LT2
for f in bam_pass1/*.bam
do
# Index the bam file
samtools index -@ 1 $f
# Output only reads that map to LT2 (to be used in subsequent steps)
samtools view -F 0x4 -@ 1 -b -o ${f%.bam}_mapped.bam $f
# As a sanity check, output unmapped melon reads
samtools view -f 0x4 -@ 1 -b -o ${f%.bam}_unmapped_sanity.bam $f
done
