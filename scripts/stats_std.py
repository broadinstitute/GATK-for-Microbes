import pysam
import pysamstats
import numpy as np
import subprocess

def threshold(reads):
	mean = np.mean(reads)
	#print(mean)
	std = np.std(reads)
	#print(std)
	thresh_pos = mean + (3 * std)
	thresh_neg = mean - (3 * std)
	return thresh_pos, thresh_neg

def mean_of_array(reads,start, pos):
	mean = np.mean(reads[start:pos], dtype=int)
	return mean

def main():
	mybam = pysam.AlignmentFile('UMB.aligned.sorted.bam')

	file = 'output_std_umb_contig_final.txt'

	with open(file, "w") as f:
		print("Index Chrom Pos Thresh Diff", file=f)
		for contig in mybam.header.references:  # Do the analysis per contig
			stats = pysamstats.load_coverage(mybam, chrom=contig)
			reads = stats.reads_pp
			position = stats.pos
			chrom = stats.chrom

			#print(len(reads))
			#print(len(position))
			#print(len(chrom))

			thresh_up, thresh_down = threshold(reads)
			#print(thresh_up)
			#print(thresh_down)

			start = 0
			big_start = 0
			last_end = 0
			state = "base"

			for i in range(1000, len(reads), 1000):
				mean = mean_of_array(reads, start, i)
				if mean >= thresh_up:
					if last_end != start:
						big_start = start
					last_end = i
					state = "on"
					prev_i = i
				else: 
					if state == "on":
						print(f'{chrom[big_start].decode("utf-8")}:{position[big_start]} {chrom[prev_i].decode("utf-8"):>20}:{position[prev_i]}')
						print(f'{chrom[big_start].decode("utf-8")}:{position[big_start]} {chrom[prev_i].decode("utf-8"):>20}:{position[prev_i]}', file=f)
						state = "off"
					else:
						state = "off"
				start = i

if __name__ == '__main__':
    main()


#a = pysamstats.load_coverage(mybam, chrom='NZ_CP060658.1')
#plt.plot(a.pos, a.reads_all)
#plt.show()