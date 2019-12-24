import matplotlib.pyplot as plt
import numpy as np
import re

fig, ax = plt.subplots()

#Start V100 plot

dsize = []
dbandwidth = []

# Collect the data from the file, ignore empty lines
data = open('scaling.data', 'r')

for line in data:
    if re.match('^Data',line):
        line = line.strip('\n')
        words = line.split()
        dummy, dummy, dummy, dsize_in, dummy, dummy, dummy, dtime_in = line.split()
        dsize.append(float(dsize_in)*8.0/1024.0/1024.0)
        dbandwidth.append(float(dsize_in)*8.0/1024.0/1024.0/float(dtime_in))
        

plt.plot(dsize, dbandwidth, "o", linestyle='-')

#axes = plt.gca() # get current axes
#axes.set_xlim([0,1000])
#axes.set_ylim([1,7000])

ax.grid()

plt.xlabel('Data Size (MB)',fontsize=16)
plt.ylabel('Bandwidth (MB/sec)',fontsize=16)

fig.tight_layout()
plt.savefig("HDF5_scaling.pdf")
plt.savefig("HDF5_scaling.svg")
plt.savefig("HDF5_scaling.png", dpi=900)

plt.show()
