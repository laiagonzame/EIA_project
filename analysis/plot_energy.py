import numpy as np
import matplotlib.pyplot as plt

# IMPORTING DATA
# TODO check the outputs files and cd to the data folder
# TODO check time and energy units
# TODO check where we are saving figures
t, ecin, epot = np.loadtxt('data/energy.data').T

plt.figure()
plt.plot(t, ecin, label='Ec')
plt.plot(t, epot, label='Ep')
plt.plot(t, ecin+epot, label='Et')
plt.xlabel('Time')
plt.ylabel('Energy / kB T')
plt.legend(loc=0)

plt.show()
