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

	# iterate over statistics, one record at a time
	#stats = pysamstats.load_coverage(mybam, chrom='NZ_CP023386.1', start= 2000000, end=3000000)
	stats = pysamstats.load_coverage(mybam)

	reads = stats.reads_all
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

	file = 'output_std.txt'

	with open(file, "w") as f:
		print("Index Chrom Pos Thresh Diff", file=f)
		for i in range(len(reads)):
			if i == ((len(reads)) - 1):
				current = reads[i]
				next_one = reads[0]
				next_pos = position[0]
			else:
				current = reads[i]
				next_one = reads[i+1]
				next_pos = position[i+1]

			if (i % 1000) == 0 and i != 0:
				#print(str(i) + ' : ' + str(position[i]))
				mean = mean_of_array(reads, start, i)
				if mean >= thresh_up:
					if last_end == start:
						print(str(big_start) + ' ' + str(i) + ' ' + str(mean) + ' +3sigma ' + str(chrom[big_start]).lstrip('b\'').rstrip('\'') + ':' + str(position[big_start]) + ' ' + str(chrom[i]).lstrip('b\'').rstrip('\'') + ':' + str(position[i]), file=f)
						print(str(big_start) + ' ' + str(i) + ' ' + str(mean) + ' +3sigma ' + str(chrom[big_start]).lstrip('b\'').rstrip('\'') + ':' + str(position[big_start]) + ' ' + str(chrom[i]).lstrip('b\'').rstrip('\'') + ':' + str(position[i]))
					else:
						print(str(start) + ' ' + str(i) + ' ' + str(mean) + ' +3sigma ' + str(chrom[start]).lstrip('b\'').rstrip('\'') + ':' + str(position[start]) + ' ' + str(chrom[i]).lstrip('b\'').rstrip('\'') + ':' + str(position[i]), file=f)
						print(str(start) + ' ' + str(i) + ' ' + str(mean) + ' +3sigma ' + str(chrom[start]).lstrip('b\'').rstrip('\'') + ':' + str(position[start]) + ' ' + str(chrom[i]).lstrip('b\'').rstrip('\'') + ':' + str(position[i]))
						big_start = start
					last_end = i
				else:
					print(str(start) + ' : ' + str(i) + ' ' + str(position[start]) + ' : ' + str(position[i]) + ' ' + str(mean))
				start = i

	c1 = subprocess.Popen(["tac {} | sort -u -k 1,1 | tac"].format(file), shell=True, stdout=open("output_std_final.txt", "w"))

if __name__ == '__main__':
    main()


#a = pysamstats.load_coverage(mybam, chrom='NZ_CP060658.1')
#plt.plot(a.pos, a.reads_all)
#plt.show()