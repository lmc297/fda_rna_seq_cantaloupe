#!/usr/bin/env python

import sys, os, glob
import multiprocessing

from multiprocessing.pool import ThreadPool

def run_trimgalore(f1):
  # Get file containing reverse reads
  f2 = f1.replace("_1.fastq.gz", "_2.fastq.gz")
  # Run trim_galore
  cmd = "trim_galore --cores 8 --fastqc --paired --gzip --stranded_illumina -o trimmed_reads_stranded {0} {1}".format(f1, f2)
  os.system(cmd)

if __name__ == "__main__":
  # Get number of threads
  pool = multiprocessing.Pool(int(sys.argv[1]))
  # Collect .fastq.gz files
  tasks = glob.glob("trimmed_fastp/*_1.fastq.gz")
  # Run run_trimgalore
  pool.map(run_trimgalore, tasks)                                            
