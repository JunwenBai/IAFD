import numpy as np
import math
import glob

stickfiles = glob.glob("*+*.txt")
indexedfiles = []

for file in stickfiles:
    index = eval(file.split('+')[0])
    indexedfiles.append((index, file))

indexedfiles = sorted(indexedfiles, key=lambda x:x[0])
fw = open("sticks.txt", "w")

for indexedfile in indexedfiles:
    file = indexedfile[-1]
    index = indexedfile[0]
    fopen = open(file, "r")
    tmpQ = []
    tmpP = []
    for line in fopen:
        line = line.strip("\n\r ")
        if line[0:1] == 'Q':
            continue
        l = line.split(' ')
        tmpQ.append(eval(l[0]))
        tmpP.append(eval(l[-1]))
    
    print
    print "index: ", index
    print len(tmpQ)
    print "Q[-1]:", tmpQ[-1]

    fw.write("Q%d=" % index)
    lenQ = len(tmpQ)
    for i, ele in enumerate(tmpQ):
        fw.write("%f" % ele)
        if i == lenQ-1:
            fw.write("\n")
        else:
            fw.write(",")

    fw.write("P%d=" % index)
    lenP = len(tmpP)
    for i, ele in enumerate(tmpP):
        fw.write("%f" % ele)
        if i == lenP-1:
            fw.write("\n")
        else:
            fw.write(",")


