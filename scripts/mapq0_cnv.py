import subprocess
import argparse
import tempfile
import os.path
import pysam
import pysamstats
import numpy as np
import sys

#checks if file exists
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file {} does not exist!".format(arg))
    else:
        return open(arg, 'r')  # return an open file handle

#takes in bam file and distance for merge function and outputs bed file showing regions with map qual 0
def merge(bam_file, distance, cutoff, output):
    with tempfile.NamedTemporaryFile(suffix='.bam') as tf:
        c0 = subprocess.run(["samtools view -q 1 -U {} -o /dev/null {}".format(tf.name, bam_file)], shell=True)
        c1 = subprocess.run(["bedtools bamtobed -i {} | bedtools merge -d {} -c 1 -o count | awk '$4 > {}'".format(tf.name, distance, cutoff)], shell=True, capture_output=True, text=True)
        for line in c1.stdout[0:-1].split('\n'):
            a = line.split('\t')
            print(a[0], a[1], a[2], 'mapq0', a[3], sep='\t', file=output)

#calculates the threshold for classifying a significant change in number of reads (+3sigma)
def threshold(reads):
    mean = np.mean(reads)
    std = np.std(reads)
    thresh_pos = mean + (3 * std)
    thresh_neg = mean - (3 * std)
    return thresh_pos, thresh_neg

#calcualtes the mean of the 1000 base region to compare to threshold value
def mean_of_array(reads,start, pos):
    mean = np.mean(reads[start:pos], dtype=int)
    return mean

#evaluates the number of reads in regions of 100 bases, if the mean of the region is greater than the threshold it is printed out
def cnv_calling(bam_file, output):
    mybam = pysam.AlignmentFile(bam_file)
    for contig in mybam.header.references:
        stats = pysamstats.load_coverage(mybam, chrom=contig)
        reads = stats.reads_pp
        position = stats.pos
        chrom = stats.chrom
        thresh_up, thresh_down = threshold(reads)
        start = 0
        big_start = 0
        last_end = 0
        state = "base"
        for i in range(1000, len(reads), 1000):
            mean = mean_of_array(reads, start, i)
            every.append(mean)
            if mean >= thresh_up:
                if last_end != start:
                    big_start = start
                last_end = i
                state = "on"
                prev_i = i
            else: 
                if state == "on":
                    print(chrom[big_start].decode("utf-8"), position[big_start], position[prev_i], 'cnv', sep='\t', file=output)
                    state = "off"
                else:
                    state = "off"
            start = i

def main():
    parser = argparse.ArgumentParser(description='Creating a bed file of map qual 0 regions and a txt file of CNV calls')
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input", dest="filename", required=True,
                    help="aligned and sorted bam file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-d", "--distance", dest="distance",
                    help="maximum distance between features allowed for features to be merged. Default is 0.", metavar="INT",
                    default=0, type= int)
    parser.add_argument("-c", "--cutoff", dest="cutoff",
                    help="minimum number of mapq0 reads present. Default is 50.", metavar="INT",
                    default=50, type= int)
    parser.add_argument(
            "-o", "--output", type=argparse.FileType('w'), dest="output", default=sys.stdout,
            help="output bed file (default: standard out)")
    args = parser.parse_args()

    merge(args.filename.name, args.distance, args.cutoff, args.output)
    cnv_calling(args.filename.name, args.output)


if __name__ == '__main__':
    main()