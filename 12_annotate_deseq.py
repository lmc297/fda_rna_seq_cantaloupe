#!/usr/bin/python

import sys, os, glob
from Bio import SeqIO

eggnog = {}
with open("reference_genome_lt2/eggnog_mapper/out.emapper.annotations", "r") as infile:
	for line in infile:
		if not line.startswith("##"):
			if "#query" in line:
				eheader = line.strip()
			else:
				seqid = line.split("\t")[0].strip()
				eggnog[seqid] = line.strip() 

d = {}
dlocus = {}
with open("reference_genome_lt2/salmonella_GCF_000006945.2.ffn", "r") as infile:
	for record in SeqIO.parse(infile, "fasta"):	
		seqdes = str(record.description).strip()
		seqid = str(record.id).strip()
		d[seqid] = seqdes
		locus = seqdes.split("locus_tag=")[1].strip()
		locus = locus.split("]")[0].strip()	
		dlocus[locus] = seqid

for f in glob.glob("./deseq_**.tsv"):
	fname = f.split("/")[-1].strip()
	with open(f, "r") as infile:
		for line in infile:
			if "baseMean" not in line:
				locus = line.split("\t")[0].split("_")[0].strip()
				if locus in dlocus.keys():
					seqid = dlocus[locus]
				else:
					seqid = "no_ffn"
				if seqid in d.keys():
					annot = d[seqid]
				else:
					annot = "no_ffn"
				if seqid in eggnog.keys():
					egg = eggnog[seqid]
				else:
					egg = "\t".join(["no_eggnog"] * len(eheader.split("\t"))).strip()
				final = line.strip() + "\t" + annot + "\t" + egg
				with open("annot_" + fname, "a") as outfile:
					print(final, file = outfile)
			else:
				final = line.strip() + "\t" + "ncbi_annotation" + "\t" + eheader
				with open("annot_" + fname, "a") as outfile:
					print(final, file = outfile)
