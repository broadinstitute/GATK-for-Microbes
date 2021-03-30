
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



LongestPath(Graph, source, sink)
  for each node b in Graph
    sb ← −∞
  ssource ← 0
  topologically order Graph
  for each node b in Graph (following the topological order)
    sb ← maxall predecessors a of node b {sa + weight of edge from a to b}
  return ssink