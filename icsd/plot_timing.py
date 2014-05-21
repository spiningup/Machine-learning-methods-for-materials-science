from pylab import *
import json
import commands

timing = json.load(open('timing.json', 'r'))
print len(timing)
timing2 = [i[2] for i in timing if i[2] < 200]
timing3 = [i for i in timing if i[2] > 40]

print len(timing2)
print len(timing3)

for item in timing3:
    print item, commands.getoutput("tail -1 %s/%s/pbserr"%(item[0], item[1]))

fig = figure(figsize=(7.5,5.2))
hist(timing2, 50)
xlabel("Computing time on 24 cores [minute]", fontsize=18)
ylabel("Number of calculations", fontsize=18)

savefig('timing-n10.pdf', bbox_inches='tight')                                  
#show()


