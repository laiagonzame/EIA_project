import os
import shutil
import numpy as np
import matplotlib.pyplot as plt


# define frozen angles
M = 100 * (np.arange(8) + 1)
M = M.astype(int)

time = []
for Mi in M:
    ifile = '../serie/data/prof' + str(Mi) + '.txt'
    data = open(ifile, 'r')
    for i in range(6):
        data.readline()
    line_n = data.readline()
    while line_n.strip():
        line = line_n
        line_n = data.readline()
    line = line.split()
    time.append(float(line[1]))

plt.loglog(M, time, 'o')
plt.plot(M, M**2/M[0]**2, '--')
plt.xlabel('Number of particles')
plt.ylabel('Execution time [s]')
plt.show()
