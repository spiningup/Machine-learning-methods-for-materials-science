import commands
import shutil

text = commands.getoutput("grep 'F m -3 m' /scratch/vstevano/RS/*/for_computing/*/icsd_*.cif")

pp = set()

for item in text.split('\n'):
    f = item.split(":")[0]
    icsdno = f[-10:-4]
    filename = f[:-4]
    pseudopot = commands.getoutput("grep '= PAW_PBE' %s/non-magnetic/POTCAR"%(filename)).split()[3]
    pp.add(pseudopot)

    oxidation_number = commands.getoutput("grep -A 2 _atom_type_oxidation_number  %s.cif"%(filename))
#    print oxidation_number

    shutil.copy2('%s/non-magnetic/CONTCAR'%(filename), 'structures/icsd_%s_CONTCAR'%(icsdno))
    shutil.copy2('%s.cif'%(filename), 'structures/icsd_%s.cif'%(icsdno))
    #print icsdno,",
#print pp
