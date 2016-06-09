import numpy as np
import matplotlib.pyplot as plt

# compile MSD file into fortran, and create a module
get_ipython().system('f2py3 -c -m statistics ../statistics.f90')

# import module
import statistics as stat

# import parameters
params = open('../data/params.data', 'r')
boxL, nhis, M, sigma, epsil, mass, dt, kBTref, tterm, stepwrite = params.readline().split()
N = int(params.readline())
M = int(M)
boxL = float(boxL)
nhis = int(nhis)
dt = float(dt)
stepwrite = int(stepwrite)
tterm = float(tterm)

rho = M / boxL**3
t = tterm + dt * stepwrite * np.arange(N)

# import trajectory
r = np.loadtxt('../data/posicions_reals.data').T
r = np.reshape(r.T, (N, M, 3))
r0 = r[0,:,:]

msd = []
# get statistics 
for i in range(N):
    msd.append(stat.statistics.desp_cuad_medio(M, r0.T, r[i,:,:].T))



fig = plt.figure(figsize=(3,3))
plt.plot(t, msd)
x1, x2, y1, y2 = plt.axis()
# plt.axis([x1, 4, y1, y2])
# plt.locator_params('x', nbins=6, tight='True')
plt.xlabel(r"t")
plt.ylabel("MSD(t)")
#plt.legend(prop={'size':10})
fig.tight_layout() 
# plt.savefig("../plots/radial_dist.pdf")

plt.show()
