import pysam
import pysamstats
import numpy as np
import matplotlib.pyplot as plt
import datetime

#x = np.arange(100)
#pos = 10
#indices = range(pos - 50, pos)
#print(x.take(indices, mode='wrap'))
#print(len(x.take(indices, mode='wrap')))

#posistheindexofthepos
def threshold(reads,pos):
	thresh = np.mean(reads[np.r_[-10000:pos]], dtype=int) / 2
	#coverage median will be new value
	#indices = range(pos-10000,pos)
	#thresh = np.mean(reads.take(indices, mode='wrap'))
	return thresh

def threshold_down(reads,pos):
	#thresh = np.mean(reads[np.r_[-10000:pos]], dtype=int)
	return thresh

def main():
	mybam = pysam.AlignmentFile('UMB.aligned.sorted.bam')

	# iterate over statistics, one record at a time
	#stats = pysamstats.load_coverage(mybam, chrom='NZ_CP023386.1', start=(2500000), end=3000000)
	stats = pysamstats.load_coverage(mybam)

	reads = stats.reads_all
	position = stats.pos
	chrom = stats.chrom

	print(len(reads))
	print(len(position))
	print(len(chrom))

	with open("output9.txt", "w") as f:
		print("Index Chrom Pos Thresh Diff", file=f)
		for i in range(len(reads)):
			if (i % 200) == 0:
				thresh = threshold(reads,i)
			if (i % 50000) == 0:
				print('currently: ' + str(i) + ' ' + str(chrom[i]).lstrip('b\'').rstrip('\'') + ' ' + str(position[i]) + ' ' + str(datetime.datetime.now().time()))
			if i == (len(reads)) - 1:
				#double check this is right
				current = reads[i]
				next_one = reads[0]
			else:
				current = reads[i]
				next_one = reads[i+1]
			if abs(next_one-current) > thresh:
				print(str(i) + ' ' + str(chrom[i]).lstrip('b\'').rstrip('\'') + ' ' + str(position[i]) + ' ' + str(thresh) + ' ' + str(next_one-current), file=f)
				print(str(i) + ' '+ str(chrom[i]).lstrip('b\'').rstrip('\'') + ' ' + str(position[i]) + ' ' + str(thresh) + ' ' + str(next_one-current))


if __name__ == '__main__':
    main()


#a = pysamstats.load_coverage(mybam, chrom='NZ_CP060658.1')
#plt.plot(a.pos, a.reads_all)
#plt.show()