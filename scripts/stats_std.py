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
	mybam = pysam.AlignmentFile('MGH283_align_sort.bam')

	file = 'output_std_mgh_contig.txt'

	with open(file, "w") as f:
		print("Index Chrom Pos Thresh Diff", file=f)
		for contig in mybam.header.references:  # Do the analysis per contig
			stats = pysamstats.load_coverage(mybam, chrom=contig)

			reads = stats.reads_pp
			position = stats.pos
			chrom = stats.chrom

			print(len(reads))
			print(len(position))
			print(len(chrom))

			thresh_up, thresh_down = threshold(reads)
			#print(thresh_up)
			#print(thresh_down)

			start = 0
			big_start = 0
			last_end = 0

			for i in range(1000, len(reads), 1000):
				if i == ((len(reads)) - 1):
					current = reads[i]
					next_one = reads[0]
					next_pos = position[0]
				else:
					current = reads[i]
					next_one = reads[i+1]
					next_pos = position[i+1]

				mean = mean_of_array(reads, start, i)
				if mean >= thresh_up:
					if last_end == start:
						big_start_chrom = str(chrom[big_start]).lstrip('b\'').rstrip('\'')
						i_chrom = str(chrom[i]).lstrip('b\'').rstrip('\'')
						print(f'{big_start:<15} {i:<15} {mean:<10} +3sigma {big_start_chrom:>20}:{position[big_start]} {i_chrom:>20}:{position[i]}')
						print(f'{big_start:<15} {i:<15} {mean:<10} +3sigma {big_start_chrom:>20}:{position[big_start]} {i_chrom:>20}:{position[i]}', file=f)
					else:
						start_chrom = str(chrom[start]).lstrip('b\'').rstrip('\'')
						i_chrom = str(chrom[i]).lstrip('b\'').rstrip('\'')
						print(f'{start:<15} {i:<15} {mean:<10} +3sigma {start_chrom:>20}:{position[start]} {i_chrom:>20}:{position[i]}')
						print(f'{start:<15} {i:<15} {mean:<10} +3sigma {start_chrom:>20}:{position[start]} {i_chrom:>20}:{position[i]}', file=f)
						big_start = start
					last_end = i
				start = i

	c1 = subprocess.Popen(["tac {} | sort -u -k 1,1 | tac".format(file)], shell=True, stdout=open("output_std_mgh_contig_final.txt", "w"))

if __name__ == '__main__':
    main()


#a = pysamstats.load_coverage(mybam, chrom='NZ_CP060658.1')
#plt.plot(a.pos, a.reads_all)
#plt.show()