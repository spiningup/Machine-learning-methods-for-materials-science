import numpy as np
import os
import commands

count = 0
dirs = os.listdir('/home/jyan/icsd')
for dir in dirs:
    number = commands.getoutput("ls /home/jyan/icsd/%s/ | wc -l"%(dir))
    print dir, number
    count += int(number)
print count
