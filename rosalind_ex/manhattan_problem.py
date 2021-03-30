
# want to find the length of a longest path from source (0, 0) to sink (n, m) 
def manhattan_tourist(n, m, down, right):
	#s = [[0]*(m+1)]*(n+1)
	s = [[0]*(m+1) for i in range(n+1)]

	#values for first column
	for i in range(1, n+1):
		s[i][0] = s[i-1][0] + down[i-1][0]

	#values for first row
	for j in range(1, m+1):
		s[0][j] = s[0][j-1] + right[0][j-1]
	
	#values for inside of matrix
	for i in range(1, n+1):
		for j in range(1, m+1):
			s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])

	#print(s)

	return s[n][m]



def main():
	n = 4
	m = 4
	down = [[1, 0, 2, 4, 3],
			[4, 6, 5, 2, 1],
			[4, 4, 5, 2, 1],
			[5, 6, 8, 5, 3]]

	right = [[3, 2, 4, 0],
			[3, 2, 4, 2],
			[0, 7, 3, 3],
			[3, 3, 0, 2],
			[1, 3, 2, 2]]

	#print(down)
	#print(right)

	max_dist = manhattan_tourist(n, m, down, right)

	#answer should be 34
	print(max_dist)

if __name__ == '__main__':
    main()