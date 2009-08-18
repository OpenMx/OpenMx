#!/usr/bin/env python
import sys
import re

class MxMatrix:
    name = None
    type = None
    free = False
    unique = False
    nrow = None
    ncol = None
    values = None
    specification = None

def parseTitle( mxInput ):
    match = re.search("^Title(.*)$", mxInput, re.MULTILINE | re.IGNORECASE)
    if match == None:
        return None
    else:
        title = match.group(1).strip()
        title = title.replace('.','_')
        return title

def parseMatrices( mxInput ):
    retval = list()
    block = re.search("Begin Matrices;(.*?)End Matrices;", mxInput, re.MULTILINE | re.IGNORECASE | re.DOTALL)
    if block == None:
        return retval
    declareLines = block.group(1).strip().split('\n')
    for declare in declareLines:
        declare = declare.strip()
        pieces = re.search("(\w+)\s+(\w+)\s+(\d+)\s+(\d+)", declare)
        if pieces != None:
            matrix = MxMatrix()
            matrix.name = pieces.group(1)
            matrix.type = pieces.group(2)
            matrix.nrow = pieces.group(3)
            matrix.ncol = pieces.group(4)
            if re.search("Free", declare, re.IGNORECASE) != None:
                matrix.free = True
            elif re.search("Unique", declare, re.IGNORECASE) != None:
                matrix.free = True
                matrix.unique = True
            retval.append(matrix)
    return retval

def printMatrices( matrices ):
    for matrix in matrices:
        outstring = "matrix" + matrix.name + " <- "
        outstring += "mxMatrix(type = \"" + matrix.type + "\", "
        outstring += "nrow = " + matrix.nrow + ', '
        outstring += "ncol = " + matrix.ncol + ', '
        outstring += "free = " + str(matrix.free).upper() + ', '
        outstring += "name = \"" + matrix.name + "\")"
        print outstring

def parseModel( mxInput ):

    # Remove all comments from the file
    mxInput = re.sub('!.*\n', '\n', mxInput)

    # Find the title of the super-model
    title = parseTitle(mxInput)

    # Find all the matrix definitions
    matrices = parseMatrices(mxInput)

    # Print matrix declarations
    printMatrices(matrices)


parseModel(sys.stdin.read())
