import commands
import os
import json
from pylada import vasp

jobdir = '/scratch/jyan/ML_natoms_10'

allcifs =  commands.getoutput("ls %s/*.cif"%(jobdir)).split('\n')
print "number of cif files", len(allcifs)

n_success = 0
errors = {"KeyError": [],
          "JobKilled": [], 
          "CantConverge": [],
          "Other": [],
          "NotCalculated": []}

for cif in allcifs:
    icsdno = cif[-10:-4]
    if not os.path.exists('%s/%s'%(jobdir, icsdno)):
        errors["NotCalculated"].append([icsdno, None])
        continue
    subdirs = os.listdir('%s/%s'%(jobdir, icsdno))
    for subdir in subdirs:
        a = vasp.Extract('%s/%s/%s'%(jobdir,icsdno, subdir))
        if a.success:
            n_success += 1
        else:
            if os.path.isfile('%s/%s/%s/pbserr'%(jobdir,icsdno, subdir)):
                err = commands.getoutput('tail -1 %s/%s/%s/pbserr'%(jobdir,icsdno, subdir))
                if "KeyError" in err:
                    errors["KeyError"].append([icsdno, subdir, err])
                elif "job killed" in err:
                    err2 = commands.getoutput("grep -B 1 'F=' %s/%s/%s/relax_cellshape/*/OSZICAR"%(jobdir, icsdno, subdir)).split('\n')
                    err3 = [j.split()[1] for j in err2 if "RMM" in j]
                    errors["JobKilled"].append([icsdno, subdir, err3])
                elif err == '':
                    err2 = commands.getoutput('grep failed %s/%s/%s/relax_cellshape/0/stdout'%(jobdir, icsdno, subdir)).split('\n')[0]
                    errors["Other"].append([icsdno, subdir, err2])
                elif "Could not converge" in err:
                    errors["CantConverge"].append([icsdno, subdir, err])
                else:
                    print icsdno, err
            else:
                pass

print "number of success", n_success

for key, item in errors.items():
    for i in item:
        print key, i


json.dump(errors, open("errors.json", "w"))
