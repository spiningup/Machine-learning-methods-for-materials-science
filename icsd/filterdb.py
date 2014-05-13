import numpy as np
import json


def filter_elements(d, el="O"):
    dnew = {}
    for key, item in d.items():
        if el in item[0]:
            dnew[key] = item
    return dnew

def filter_nel(d, l="2"):
    dnew = {}
    for key, item in d.items():
        formula = item[0].split()
        if len(formula) == l:
            dnew[key] = item
    return dnew

def filter_calculated(d, d0):
    dnew = {}
    for key, item in d.items():
        if key in d0:
            print key
            continue
        else:
            dnew[key] = item
    return dnew

if __name__ == "__main__":

    d = json.load(open('icsd.json', 'r'))
    dtmp = json.load(open('icsdnum_NRELDB.json', 'r'))
    d0 = []
    for i in dtmp:
        icsdno = i.values()[0]
        icsdstr = "{0:06d}".format(int(icsdno))
        d0.append(icsdstr)

    d = filter_calculated(d,d0)
    print len(d)
    json.dump(d, open("icsd_notinDB.json", 'w'))


