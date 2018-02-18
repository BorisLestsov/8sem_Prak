import numpy as np

def norm(a0):
	return np.sqrt(a0.dot(a0))


with open("data/mat.txt", 'r') as f:
	f.readline()
	a = []
	for i in f:
		a.append([float(e) for e in i.split()])




mat = np.array(a)
print mat

for ind in range(3):

	I = np.zeros(shape=(3,3))
	for i in range(3):
		I[i,i] = 1

	a0 = mat[:, ind].copy()

	norm_a0 = norm(a0)
	x = a0.copy()

	x[ind] -= norm_a0

	x /= norm(x)

	U = I - 2*np.outer(x,x)

	for j in range(ind):
		U[:, j] = 0
		U[j, :] = 0
		U[j, j] = 1
	print U

	mat = U.dot(mat)

	print mat