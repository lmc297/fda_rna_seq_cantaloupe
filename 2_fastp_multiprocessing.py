#!/usr/bin/env python

import sys, os, glob
import multiprocessing

from multiprocessing.pool import ThreadPool

def run_fastp(f1):
    # Get reverse reads file
    f2 = f1.replace("_1.fastq.gz", "_2.fastq.gz")
    # Get sample basename
    prefix = f1.split("/")[-1].strip()
    # Name output trimmed read files
    o1 = "trimmed_fastp/" + prefix.replace("_1.fastq.gz", "_trimmed_1.fastq.gz")
    o2 = o1.replace("_trimmed_1.fastq.gz", "_trimmed_2.fastq.gz")
    # Name log files
    prefix = "fastp_logs/" + prefix.replace("_1.fastq.gz", "")
    # Run fastp
    cmd = "fastp -i {0} -I {1} -o {2} -O {3} --detect_adapter_for_pe --length_required 36 --trim_poly_x --correction --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --thread 8 -j {4} -h {5}".format(f1, f2, o1, o2, prefix + ".json", prefix + ".html")
    os.system(cmd)

if __name__ == "__main__":
    # Get number of threads
    pool = multiprocessing.Pool(int(sys.argv[1]))
    # Make results directories
    os.mkdir("trimmed_fastp")
    os.mkdir("fastp_logs")
    # Collect raw reads
    tasks = glob.glob("raw_reads/*_1.fastq.gz")
    # Run run_fastp
    pool.map(run_fastp, tasks)
