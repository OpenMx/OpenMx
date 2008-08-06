#!/usr/bin/env python

import sys

total = 0
argstr = sys.argv
i = 0
for a in argstr:
    if a == '-mods':
        total = int(argstr[i+1])
    i += 1
for perm in range(0, total):
    print("["+str(perm)+"]"+" "+str(perm)+"_model.rdata")

    



