
def lcs_backtrack(v, w):
	backtrack = [[None]*(len(w)+1) for i in range(len(v)+1)]
	s = [[None]*(len(w)+1) for i in range(len(v)+1)]

	for i in range(len(v) + 1):
		s[i][0] = 0
	for j in range(len(w) + 1):
		s[0][j] = 0
	for i in range(1, len(v) + 1):
		for j in range(1, len(w) + 1):
			match = 0
			if v[i-1] == w[j-1]:
				match = 1
			s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + match)
			if s[i][j] == s[i-1][j]:
				backtrack[i][j] = "↓"
			elif s[i][j] == s[i][j-1]:
				backtrack[i][j] = "→"
			elif s[i][j] == s[i-1][j-1] + match:
				backtrack[i][j] = "↘"

	#xprint(backtrack)
	return backtrack

#A longest common subsequence of s and t
#OutputLCS(Backtrack, v, |v|, |w|).
def output_lcs(backtrack, v, i, j):
	if i == 0 or j == 0:
		return ''
	if backtrack[i][j] == "↓":
		return output_lcs(backtrack, v, i - 1, j)
	elif backtrack[i][j] == "→":
		return output_lcs(backtrack, v, i, j - 1)
	else:
		return output_lcs(backtrack, v, i - 1, j - 1) + v[i]


def main():

	v = 'AACTTGG' #AACCTTGG
	w = 'ACACTGTGA' #ACACTGTGA

	#v = 'TACACCAAGCAAACGTGTCTCGTAAGCCAACTGGGGGACAATCGTGGGACTTTCCTCAGTGAAGGCGCAGTGTTGGCGCTCGCACTGTGCGTGATCGTGTACTGCGATTGGAGAGCCTCCTACCGGCAGTGTCTTGGGTTAAGTAGCTGGACTGAAATCGCGTAGCAATAGTGTGGACAATAGACAAAAGGAGCTGGGCTATATCCCGGATGATAGGGGTATTATTGTTTGCCCCCCTTTGGCAAAGCCCTTGCCTTTAATCCATCGCAGATTCGGAAGTATCTGACAGCAGATTAGCGAGGCCGCGAGTAGCGATAATTCACGCCTACATCCGTCCTGATGTGGCTGTCAAATTGAGTCGAGTGATCGGTCGTAGGGTAATGAGGAAGTGTCACAAGGGTGGGCAAAAACAATGCAAACTCTCGAATGCAGGCGGTCTGCCCTAGGTGGAAGAAACCGCGTAGGTATGGCTCCCGAACTGCCGTTGTTCGGGGTTCGATTAGGGGCTCGCGATCCGGCGCGCGCGAGGCACGAGATTGATATTAGCCGAATGAGTCGAACAATGTCTCAGGGATCGCCGTGGAAGGCTGGTTCATATGTCAAAGATTGAATTTGAGAACACCACCATGCCTGTGAGGTCAATCTCCAACTAATCAAGGCCTGGTTCAAAGAAACCTAATGTCATGCTGACTGCGCGTTCCATAGAGGGGTCGGCCATGCACAACGGGCCTATTCTAGGCATATCGATCCAACGGTGGATCCCGTTAAGTATAACTATCCCAGACGTAGTAAGATCGAGCGATTGACTCGAACGGCCCACCACTTATTTTCACGCCCCTGCCATCCATGCACCTCAACGTGAACGTGCCGTTAGGTGAGAACACAGGTTGCCCCGCCGCATTCACATAAACTGACGTAGACTTGGGCAGCAGACTTAAAACTCATGAATTAGTCAGAATCGCCCGCGATAATGTCTG'
	#w = 'TCCCCCCGTCTGACCCAGGATGAGAAGTGTGTTACGCTCCCGCAGTAAGCCAGTCGATATTATTCGCCTTTCGAAGGGACACGATTTAAATTACTATCCAGCTACGCCCTCACTACAGTGTTGCGCTTGCATTGATGACCATGCCCCACGGAGAGACAATGAGTGACTTTCTTAGTGCTTGACCTTGTTTATGTTCTAAAGCCGTATGTTAGCGCGCCTGTAAACCATGCGGAAGGCAGGGGGTCTTCGCGTACGCGAGAGGCTGTGGAAAAAAGATTTCCGGTGTTGTAGCGCGAACCCTGTAACCAACGGGAGATGCATTTCCGCGAGATATGAACCAATAAAATTGCGCAGGAGGGGTAGGCGGGAGCGAATGGGCGCCCCCACTGGAGGAGCGAAGCACTTAATTGATCCTCATGATCATAATCCGAAGCGGGGGAAATTGTAGGTCGCATGTCCGCCCGGAGGATAAGTAAACCTCCCACCTGGTTACCAAACCGGTCAGCTGTGGCGCAACGCAAAACCCTAGTATGTCGCAGAGGCAATATATCGCGCGGTTTTAATCGAGCCTTGTATGCTCGAGACGAGACGGTCCAGTTGGATTGGGTGCTGTCGGACATCCGGATATGGTTGCCTTAAGCACAATTTGTCTGCACTCAGTCTCTCAACAAGCGCGACTCTTACTTATTAATCCAGTTAGCTCGACTAATTGGAGCCACAGTAACAGTCCCTCTTAAGAAGGCAGGCTGGTATCCTGCGCCTCCGAAGCCCAGTAGCAGTAGGGGCGAATGATTACAATCCTGGTACGATTCACCAGGGCTCTTGTTCGGTTCATCTTCTGTGCCCCCGTCCCGTCAAGGGGGTCAGACGAAACGTATTGCCTAGTCCTTTGTTCGCTCTATAAGGTGGCCAGGGGGTTAACGGGAGCTCAAGAACCGAGGGGCCACTGATTTGTACACCTGTTGCATCCACGAATGTCAGTCTGAGACTGCTTG'
	
	backtrack = lcs_backtrack(v, w)
	lcs = output_lcs(backtrack, v, len(v), len(w))
	print(lcs)

	#expected out "AACTGG"

if __name__ == '__main__':
    main()