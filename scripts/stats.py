import pysam
import pysamstats
import numpy as np
import matplotlib.pyplot as plt

#hello = [1,2,3,4,5,6]
#print(hello[2:6])

#posistheindexofthepos
def threshold(reads,pos):
	total = 0
	if pos < 10000:
		num = 10000 - pos
		ending = reads[-num:]
		for i in ending:
			total += i
		for j in range(pos):
				total += reads[j]
	elif pos > 10000:
		num = pos - 10000
		prev = reads[num:pos]
		for i in prev:
			total += i
	else:
		for i in range(pos):
			total +=reads[i]
	thresh = (total / 10000) / 2
	return thresh

def main():
	mybam = pysam.AlignmentFile('UMB.aligned.sorted.bam')

	chrom = []
	position = []
	reads = []

	# iterate over statistics, one record at a time
	for rec in pysamstats.stat_coverage(mybam, chrom='NZ_CP023386.1', start=(2848699-10500), end=2848710):
		#print(rec['chrom'], rec['pos'], rec['reads_all'])
		chrom.append(rec['chrom'])
		position.append(rec['pos'])
		reads.append(rec['reads_all'])

	print(len(chrom))
	print(len(position))
	print(len(reads))

	with open("output2.txt", "w") as f:
		print("Chrom Pos Thresh Diff", file=f)
		for i in range(len(reads)):
			thresh = threshold(reads,i)
			if i == (len(reads)) - 1:
				#double check this is right
				current = reads[i]
				next_one = reads[0]
			else:
				current = reads[i]
				next_one = reads[i+1]
			if abs(next_one-current) > thresh:
				#print(str(i) + ' ' + str(chrom[i]) + ' ' + str(position[i]) + ' ' + str(thresh) + ' ' + str(next_one-current), file=f)
				print(str(i) + ' '+ str(chrom[i]) + ' ' + str(position[i]) + ' ' + str(thresh) + ' ' + str(next_one-current))


if __name__ == '__main__':
    main()


#a = pysamstats.load_coverage(mybam, chrom='NZ_CP060658.1')
#plt.plot(a.pos, a.reads_all)
#plt.show()