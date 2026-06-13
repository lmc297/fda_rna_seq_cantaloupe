#!/usr/bin/env python

import sys, os, glob
import multiprocessing

from multiprocessing.pool import ThreadPool

# run FastQC on each trimmed file (*.fastq.gz)
def run_fastqc(f):
  cmd = "fastqc -t 1 {0}".format(f)
  os.system(cmd)

if __name__ == "__main__":
  # Specify number of threads
  pool = multiprocessing.Pool(int(sys.argv[1]))
  # Collect trimmed fastq.gz files
  tasks = glob.glob("trimmed_fastp/*.fastq.gz")
  # Run run_fastqc
  pool.map(run_fastqc, tasks)
