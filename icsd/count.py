import numpy as np
import os
import commands

dirs = os.listdir('/home/jyan/icsd')
for dir in dirs:
    print dir, commands.getoutput("ls /home/jyan/icsd/%s/ | wc -l"%(dir))
