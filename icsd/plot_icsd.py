import json
import numpy as np
from collections import defaultdict
import operator
from pylab import *
from filterdb import *
from visualize import *
from atomic_constants import *

d = json.load(open('icsd_includedcalced.json', 'r'))
d = filter_elements(d, "D", incl=False)
elnum, els, natoms  = get_elements_statistics(d)

n1 = {key: value for key, value in natoms.items() if value > 200}
print len(n1)/ float(len(natoms))

atomic_num = {el:atomic_number[el] for el in els.keys()}
atomic_num = sorted(atomic_num.iteritems(), key=operator.itemgetter(1))

fig = figure(figsize=(10,5))
els["O"] /= 2
elements = [(i[0], els[i[0]]) for i in atomic_num]
X = np.arange(len(elements))
Y = [i[1] for i in elements]
xll = [i[0] for i in elements]
bar(X, Y, align="center", width=0.5)
xticks([-2, max(X)+2])
xlim(-1,max(X)+1)
ylim(0, max(Y)+20)
for i in range(len(elements)):
    if i %2 == 0:
        text(i, -170, xll[i], verticalalignment="center", horizontalalignment="center")
    else:
        text(i, 7400, xll[i], verticalalignment="center", horizontalalignment="center")
text(7, 6800, "/2")
ylabel("Number of ICSD Entries")
gca().yaxis.labelpad = 10
savefig('els-hist.pdf', bbox_inches='tight')                                  
    

#print len(d)
#print elnum
#print sorted(els.iteritems(), key=operator.itemgetter(1))
#print sorted(natoms.iteritems(), key=operator.itemgetter(0))

fig = figure(figsize=(6,4))
bar(natoms.keys(), natoms.values())
xlim(0,200)
xlabel("Number of Atoms in Primitive Cell", fontsize=18)
ylabel("Number of ICSD Entries", fontsize=18)
savefig('natoms-hist.pdf', bbox_inches='tight')                                  

fig = figure(figsize=(6.5,4))
bar(np.array(elnum.keys())-0.4, elnum.values())
xticks(elnum.keys())
print elnum
#xticks([-1, 10])
#xlim(0, 9)
#for i in range(1, 10):
#    text(i, -1000, i)
xlabel("Number of Elements per ICSD Entry", fontsize=18)
ylabel("Number of ICSD Entries", fontsize=18)
savefig('nels-hist.pdf', bbox_inches='tight')                                  
#show()


