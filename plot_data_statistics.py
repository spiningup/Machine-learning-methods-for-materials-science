from read_json import *
import sys
from pylab import *
import numpy as np

mset = read_json(sys.argv[1])

volerror = []
for atoms in mset:
    if atoms.exptvol is not None:
        volerror.append((atoms.calcvol - atoms.exptvol) / atoms.exptvol * 100)
print np.mean(volerror), np.min(volerror), np.max(volerror)

fig = figure(figsize=(7.5,5.2))
hist(volerror, 50)
xlabel("Error in calculated volume [%]", fontsize=18)
ylabel("Number of ICSD entries", fontsize=18)
savefig('volume_error.pdf', bbox_inches='tight')                                  


#################################################
Eref = attribute_tolist(mset, attr="Eref")
fig = figure(figsize=(7.5,5.2))
hist(Eref, 50)
xlabel("Cohesive Energy [eV]", fontsize=18)
ylabel("Number of ICSD entries", fontsize=18)
savefig('hist_Ecoh.pdf', bbox_inches='tight')                                  





#show()
