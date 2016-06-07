import numpy as np
import matplotlib.pyplot as plt

# IMPORTING DATA
# TODO check time and energy units
# TODO check where we are saving figures

# compile RDF file into fortran, and create a module
# get_ipython().system('f2py3 -c -m ener_temp ../temperature.f90')

# import module
import ener_temp as stat

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

nf = 3 * M - 3
rho = M / boxL**3
t = tterm + dt * stepwrite * np.arange(N)

# import trajectory
r = np.loadtxt('../data/posicions.data').T
v = np.loadtxt('../data/velocitats.data').T
r = np.reshape(r.T, (N, M, 3))
v = np.reshape(v.T, (N, M, 3))

temp = []
ecin = []
# get temperature
for i in range(N):
    temp.append(stat.temperature.compute_temperature(v[i,:,:].T, M))
    ecin.append(stat.temperature.kinetic_energy(nf,temp[i]))

plt.figure()
plt.plot(t, ecin, label='Ec')
# plt.plot(t, epot, label='Ep')
# plt.plot(t, ecin+epot, label='Et')
plt.xlabel('Time')
plt.ylabel('Energy / kB T')
plt.legend(loc=0)

plt.show()
