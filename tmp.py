import numpy as np


comm_size = 3
size = 3

my_cols = {rank: size / comm_size for rank in range(comm_size)}

my_cols_ind = {rank:[rank+comm_size*i for i in range(my_cols[rank])] for rank in range(comm_size)}
for r in range(comm_size):
	if r < size%comm_size:
		my_cols[r] += 1
		my_cols_ind[r].append(r+(my_cols[r]-1)*comm_size)

for r in range(comm_size):
	print r
	for ind in my_cols_ind[r]:
		print "    ", ind

print '-'*10

for ind in range(size)[::-1]:
	print ind
	for r in range(comm_size):
		print "    ", r, [i for i in range(ind+1, size) if i in my_cols_ind[r]]
