import numpy as np
import matplotlib.pyplot as plt

# compile RDF file into fortran, and create a module
get_ipython().system('f2py3 -c -m statistics ../statistics.f90')

# import module
import statistics as stat

# import parameters
params = open('../data/params.data', 'r')
boxL, nhis, M, sigma, epsil, mass, dt = params.readline().split()
N = int(params.readline())
M = int(M)
boxL = float(boxL)
nhis = int(nhis)

rho = M / boxL**3

# import trajectory
r = np.loadtxt('../data/posicions.data').T
r = np.reshape(r.T, (N, M, 3))

# initialize g(r)
dr, gr = stat.statistics.declarate_radial_dist(boxL, nhis)

# accumulate statistics 
for i in range(N):
    gi = stat.statistics.accumulate_radial_dist(M, boxL, r[i,:,:].T, nhis, dr)
    gr += gi

gg = stat.statistics.compute_radial_dist(nhis, dr, N, rho, M, gr)
rhis = dr * np.arange(nhis)


fig = plt.figure(figsize=(3,3))
plt.plot(rhis[:-1], gg[:-1])

x1, x2, y1, y2 = plt.axis()
# plt.axis([x1, 4, y1, y2])
# plt.locator_params('x', nbins=6, tight='True')
plt.xlabel(r"r [$\sigma$]")
plt.ylabel("g(r)")
plt.legend(prop={'size':10})
fig.tight_layout() 
# plt.savefig("../plots/radial_dist.pdf")

plt.show()
