import commands
import os
import json
from pylada import vasp

jobdir = os.getcwd()

allcifs =  commands.getoutput("ls %s/*.cif"%(jobdir)).split('\n')
print "number of cif files", len(allcifs)

timing = []
for cif in allcifs:
    icsdno = cif[-10:-4]
    if not os.path.exists('%s/%s'%(jobdir, icsdno)):
        continue
    subdirs = os.listdir('%s/%s'%(jobdir, icsdno))

    for subdir in subdirs:
        a = vasp.Extract('%s/%s/%s'%(jobdir,icsdno, subdir))
        if a.success:
            computing_time = os.path.getctime('%s/%s/%s/pbsout'%(jobdir, icsdno, subdir)) - os.path.getctime('%s/%s/%s/relax_cellshape/0/INCAR'%(jobdir, icsdno, subdir))
#            print '%s/%20s : %6.2f'%(icsdno, subdir, computing_time / 60.) # in minute
            timing.append((icsdno, subdir, computing_time/60.))

json.dump(timing, open("timing.json", "w"))
