import commands
import os
import json
from pylada import vasp

jobdir = '/scratch/jyan/ML_natoms_10'

allcifs =  commands.getoutput("ls %s/*.cif"%(jobdir)).split('\n')
print "number of cif files", len(allcifs)

for cif in allcifs:
    icsdno = cif[-10:-4]
    if not os.path.exists('%s/%s'%(jobdir, icsdno)):
        continue
    subdirs = os.listdir('%s/%s'%(jobdir, icsdno))
    for subdir in subdirs:
        a = vasp.Extract('%s/%s/%s'%(jobdir,icsdno, subdir))
        if a.success:
            print "%s/%s"%(icsdno,subdir)
#            os.system("rm -rf %s/%s/WAVECAR")
#            os.system("rm -rf %s/%s/CHGCAR")
#            os.system("rm -rf %s/%s/CHG")
