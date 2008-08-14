#!/usr/bin/env python

import sys

def fact(x): return reduce(lambda x, y: x*y, xrange(2, x+1))

conn = 1
size = 1
total = 1
argstr = sys.argv

argfile = open("mapper.args", "w")


i = 0
for a in argstr:
    argfile.write(a+" ")    
    if a == '-conn':
        conn = int(argstr[i+1])
    elif a == '-size':
        size = int(argstr[i+1])
    i += 1

if (conn <2):
    total = size

else:
    total = fact(size)/(fact(conn)*fact(size-conn))    


for perm in range(0, total):
    print("["+str(perm)+"]"+" "+str(conn)+"_"+str(perm)+"_model.rdata")            


    



