from pylab import *
import json

timing = json.load(open('timing.json', 'r'))
print len(timing)
timing2 = [i[2] for i in timing if i[2] < 200]
timing3 = [i[2] for i in timing if i[2] > 200]

print len(timing2)
print timing3

hist(timing2, 50)
show()

