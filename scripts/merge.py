import subprocess
import argparse
import tempfile
import os.path

#checks if file exists
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file {} does not exist!".format(arg))
    else:
        return open(arg, 'r')  # return an open file handle

#takes in bam file and distance for merge function and outputs bed file showing regions with map qual 0
def merge(bam_file, distance):
	tf = tempfile.NamedTemporaryFile(delete=False, suffix='.bam')
	c0 = subprocess.run(["samtools view -q 1 -U {} -o /dev/null {}".format(tf.name, bam_file)], shell=True)
	c1 = subprocess.Popen(["bedtools bamtobed -i {} | bedtools merge -d {}".format(tf.name, distance)], shell=True, stdout=open("merged.bed", "w"))


def main():
	parser = argparse.ArgumentParser(description='Creating a bed file of map qual 0 regions')
	required = parser.add_argument_group('required arguments')
	required.add_argument("-i", dest="filename", required=True,
                    help="aligned and sorted bam file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
	parser.add_argument("-d", dest="integer",
                    help="maximum distance between features allowed for features to be merged. Default is 0.", metavar="INT",
                    default=0, type= int)
	args = parser.parse_args()

	bam_file = args.filename.name
	distance = args.integer

	merge(bam_file, distance)


if __name__ == '__main__':
    main()