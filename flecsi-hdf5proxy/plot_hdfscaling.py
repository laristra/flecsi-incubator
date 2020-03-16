import matplotlib.pyplot as plt
import numpy as np
import re
import glob, os
import getpass

fig, ax = plt.subplots()

read32_dsize = []
read32_dbandwidth = []
write32_dsize = []
write32_dbandwidth = []
read64_dsize = []
read64_dbandwidth = []
write64_dsize = []
write64_dbandwidth = []

# Collect the data from the file, ignore empty lines

USER=getpass.getuser()
FILE_PATH="/lustre/ttscratch1/" + USER + "/hdf5proxy/slurm*"

for file in glob.glob(FILE_PATH):
   print file
   data = open(file, 'r')
   for line in data:
       line = line.strip('\n')
       if re.match('^srun',line):
           print line
           line_save = line
       if re.match('^Time',line):
           words_save = line_save.split()
           nodes_in = words_save[2]
           nrankspernode_in = words_save[4]
           print line
           words = line.split()
           dsize_in = words[9]
           dbandwidth_in = words[16]
           if re.match('32',words_save[4]):
              if re.match('Read',words[13]):
                 read32_dsize.append(float(dsize_in)*1024.0/1024.0)
                 read32_dbandwidth.append(float(dbandwidth_in))
                 print ""
              if re.match('Write',words[13]):
                 write32_dsize.append(float(dsize_in)*1024.0/1024.0)
                 write32_dbandwidth.append(float(dbandwidth_in))
           if re.match('64',words_save[4]):
              if re.match('Read',words[13]):
                 read64_dsize.append(float(dsize_in)*1024.0/1024.0)
                 read64_dbandwidth.append(float(dbandwidth_in))
                 print ""
              if re.match('Write',words[13]):
                 write64_dsize.append(float(dsize_in)*1024.0/1024.0)
                 write64_dbandwidth.append(float(dbandwidth_in))
        
read32_dsize, read32_dbandwidth = zip(*sorted(zip(read32_dsize, read32_dbandwidth)))
read64_dsize, read64_dbandwidth = zip(*sorted(zip(read64_dsize, read64_dbandwidth)))
write32_dsize, write32_dbandwidth = zip(*sorted(zip(write32_dsize, write32_dbandwidth)))
write64_dsize, write64_dbandwidth = zip(*sorted(zip(write64_dsize, write64_dbandwidth)))
plt.plot(read32_dsize, read32_dbandwidth, "o", linestyle='-', label="read32")
plt.plot(write32_dsize, write32_dbandwidth, "x", linestyle='--', label="write32")
plt.plot(read64_dsize, read64_dbandwidth, "o", linestyle='-', label="read64")
plt.plot(write64_dsize, write64_dbandwidth, "x", linestyle='--', label="write64")

#axes = plt.gca() # get current axes
#axes.set_xlim([0,1000])
#axes.set_ylim([1,7000])

ax.grid()

plt.xlabel('Data Size (MB)',fontsize=16)
plt.ylabel('Bandwidth (MB/sec)',fontsize=16)
plt.legend(loc="upper left")

fig.tight_layout()
plt.savefig("HDF5_scaling.pdf")
plt.savefig("HDF5_scaling.svg")
plt.savefig("HDF5_scaling.png", dpi=900)

plt.show()
