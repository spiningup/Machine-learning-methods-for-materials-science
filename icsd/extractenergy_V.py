import commands
import os
import json
import numpy as np
from pylada import vasp

jobdir = '/scratch/vstevano/ICSD/O'
allformula = os.listdir('%s'%(jobdir))
print "number of formulas", len(allformula)

for formula in allformula:
    if "pylada" in formula or \
            "mag" in formula or \
            ".py" in formula or \
            ".pyc" in formula or \
            ".txt" in formula or \
            "tmp" in formula or \
            "pickle" in formula: 
        continue
    allcifs = os.listdir('%s/%s'%(jobdir, formula))
    for cif in allcifs:
        if ".cif" in cif: continue
        icsdno = cif[-6:]
        if not os.path.exists('%s/%s/icsd_%s'%(jobdir, formula, icsdno)):
            continue
        subdirs = os.listdir('%s/%s/icsd_%s'%(jobdir, formula, icsdno))
        E = []
        dirs = []
        for subdir in subdirs:
            a = vasp.Extract('%s/%s/icsd_%s/%s'%(jobdir,formula, icsdno, subdir))
            if a.success:
                E.append(a.total_energy.magnitude)
                dirs.append(subdir)
                #print "%s/%s"%(icsdno,subdir), a.total_energy
    #    if len(E) > 0: print icsdno, min(E)
        if len(E) > 1:
            minidx = np.argsort(E)[0]
            minE = min(E)
            minEdir = dirs[minidx]
            E.pop(minidx)
            dirs.pop(minidx)
            if minEdir == "non-magnetic":
                if np.mean(E) - minE < 1: # print the ones with energy difference smaller than 1 eV
                    print "%s/%20s,%6.2f, %3d, %6.2f, %6.2f"%(icsdno, minEdir, minE, len(E)+1, np.mean(E), np.std(E))
            else:
                nonmagidx = [j for j in range(len(dirs)) if dirs[j] == "non-magnetic"][0]
                nonmagE = E[nonmagidx]
                E.pop(nonmagidx)
                dirs.pop(nonmagidx)
                hsE = []
                lsE = []
                hsdir = []
                lsdir = []
                for dd, ee in zip(dirs, E):
                    if "hs" in dd:
                        hsE.append(ee)
                        hsdir.append(dd)
                    elif "ls" in dd:
                        lsE.append(ee)
                        lsdir.append(dd)

                if np.mean(E) - minE > 0.1 and np.std(E) > 0: 
                    if len(hsE) > 0:
                        print "%s/%20s, %7.2f, %3d, non-magnetic, %7.2f, %7.2f, %7.2f, %7.2f, %7.2f"%(icsdno, minEdir, minE, len(E)+2, nonmagE, np.mean(hsE), np.std(hsE), np.mean(lsE), np.std(lsE))
                    else:
                        print "%s/%20s,%7.2f, %3d, non-magnetic, %7.2f, %7.2f, %7.2f"%(icsdno, minEdir, minE, len(E)+2, nonmagE, np.mean(E), np.std(E))




