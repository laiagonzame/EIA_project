import numpy as np
import matplotlib.pyplot as plt

# compile RDF file into fortran, and create a module
get_ipython().system('f2py3 -c ../statistics.f90 -m statistics')

# import module
import statistics as stat

# TODO import parameters
boxL = 10
nhis = 200
N = 5
M = 400

rho = M / boxL**3

# test
r = np.array(np.random.randn(3,M,N), order='FORTRAN')
#r = np.array(np.random.randn(3,M,N), order='C')

# TODO import trajectory

dr, gr = stat.statistics.declarate_radial_dist(boxL, nhis)

# accumulate statistics 
for i in range(N):
    #gi = np.zeros(nhis)
    #stat.statistics.accumulate_radial_dist(M, boxL, r[:,:,i], nhis, dr, gi)

    gi = stat.statistics.accumulate_radial_dist(M, boxL, r[:,:,i], nhis, dr)
    gr += gi

gr = compute_radial_dist(nhis, dr, N, rho)
r = dr * np.arange(nhis)

# 
# fig = plt.figure(figsize=(3,3))
# plt.plot(r[:-1], gr[:-1])
# 
# x1, x2, y1, y2 = plt.axis()
# # plt.axis([x1, 4, y1, y2])
# # plt.locator_params('x', nbins=6, tight='True')
# plt.xlabel(r"r [$\sigma$]")
# plt.ylabel("g(r)")
# plt.legend(prop={'size':10})
# fig.tight_layout() 
# # plt.savefig("../plots/radial_dist.pdf")
# 
# plt.show()
