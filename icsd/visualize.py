import json
import numpy as np
from collections import defaultdict
import operator
from pylab import *
from filterdb import *

def get_elements_statistics(d):
    elnum = defaultdict(int)
    els = defaultdict(int)
    natoms = defaultdict(int)
    for key, item in d.items():
        formula = item[0].replace("'",'').split()
        elnum[len(formula)] += 1
        if item[2] == 0: print key, item
        natoms[item[2]] += 1

        for el in formula:
            for i in el:
                if i.isdigit():
                    el = el.replace(i, "")
            els[el] += 1

    return elnum, els, natoms

def plot_dict(a):
    X = np.arange(len(a))
    bar(X, a.values(), align="center", width=0.5)
    xticks(X, a.keys())
    ymax = max(a.values()) + 20
    ylim(0, ymax)
#    show()


if __name__ == "__main__":
    d = json.load(open('icsd_includedcalced.json', 'r'))
#    d = filter_natoms(d, nmax=10)
    elnum, els, natoms  = get_elements_statistics(d)
    print len(d)
    print elnum
    print sorted(els.iteritems(), key=operator.itemgetter(1))
    print sorted(natoms.iteritems(), key=operator.itemgetter(0))
    bar(natoms.keys(), natoms.values())
    show()

#    d = filter_nel(d, 1)

#    for el in els.keys():
#        dnew = filter_elements(d, el)
#        elnumnew, elsnew = get_elements_statistics(dnew)
#        print el, elnumnew
#        plot_dict(elnumnew)
