import numpy as np

N = 10
LENGTH = 10

A = np.zeros([N, LENGTH],dtype=np.int)
f = open("minitest.txt", "r")

f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
for i in range(N):
	# for j in range(N):
	B = list(f.readline())
	for j in range(LENGTH):
		A[i,j] = int(B[j])

np.savetxt('minitest.out', A.transpose(), delimiter = '', fmt = "%d")
