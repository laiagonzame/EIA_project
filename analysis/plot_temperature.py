import numpy as np
import matplotlib.pyplot as plt

# IMPORTING DATA
# TODO check the outputs files and cd to the data folder
# TODO check time and temperature units
# TODO check where we are saving figures
t, temp, ecin = np.loadtxt('../data/energy-temp.data').T

plt.figure()
plt.plot(t, temp)
plt.xlabel('Time')
plt.ylabel('Temperature')
# plt.savefig('plots/temperature.pdf')
plt.show()
