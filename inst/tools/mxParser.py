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

def parseDefines( mxInput ):
    defines = {}
    matchIter = re.finditer("#define\s+(\S+)\s+(\S+)", mxInput)
    for match in matchIter:
        defines[match.group(1)] = match.group(2)
    return defines

def parseTitle( mxInput ):
    match = re.search("^Title(.*)$", mxInput, re.MULTILINE | re.IGNORECASE)
    if match == None:
        return None
    else:
        title = match.group(1).strip()
        title = title.replace('.','_')
        return title

def parseMatrices( mxInput, defines ):
    retval = {}
    block = re.search("Begin Matrices;(.*?)End Matrices;", 
    	mxInput, re.MULTILINE | re.IGNORECASE | re.DOTALL)
    if block == None:
        return retval
    declareLines = block.group(1).strip().split('\n')
    for declare in declareLines:
        declare = declare.strip()
        pieces = re.search("(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", declare)
        if pieces != None:
            matrix = MxMatrix()
            matrix.name = pieces.group(1)
            matrix.type = pieces.group(2)
            try:
                matrix.nrow = int(pieces.group(3))
            except:
                matrix.nrow = int(defines[pieces.group(3)])
            try:
                matrix.ncol = int(pieces.group(4))
            except:
                matrix.ncol = int(defines[pieces.group(4)])
            if re.search("Free", declare, re.IGNORECASE) != None:
                matrix.free = True
            elif re.search("Unique", declare, re.IGNORECASE) != None:
                matrix.free = True
                matrix.unique = True
            retval[matrix.name] = matrix
    return retval    

def printMatrices( matrices ):
    startCounter = 1
    for matrix in matrices.values():
        outstring = "matrix" + matrix.name + " <- "
        outstring += "mxMatrix(type = \"" + matrix.type + "\", "
        outstring += "nrow = " + str(matrix.nrow) + ", "
        outstring += "ncol = " + str(matrix.ncol) + ", "
        outstring += "free = " + str(matrix.free).upper() + ", "
        nextcounter = None
        if matrix.free and not matrix.unique and matrix.specification == None:
            endCounter = startCounter + specLength[matrix.type](matrix.nrow, matrix.ncol) - 1
            outstring += "labels = " + str(startCounter) + " : "
            outstring += str(endCounter) + ", "
            startCounter = endCounter + 1
        outstring += "byrow = TRUE, "
        outstring += "name = \"" + matrix.name + "\")"
        print outstring
        print

specLength = {   "Diag"  : lambda row, col: row,
                 "SDiag" : lambda row, col: row * (row - 1) / 2,
                 "Stand" : lambda row, col: row * (row - 1) / 2,
                 "Symm"  : lambda row, col: row * (row + 1) / 2,
                 "Lower" : lambda row, col: row * (row + 1) / 2,
                 "Full"  : lambda row, col: row * col }
    

def parseModel( mxInput ):

    # Remove all comments from the file
    mxInput = re.sub('!.*\n', '\n', mxInput)

    # Find the title of the super-model
    title = parseTitle(mxInput)

    # Find all the #define statements
    defines = parseDefines(mxInput)

    # Find all the matrix definitions
    matrices = parseMatrices(mxInput, defines)

    # Find all the matrix specifications
#    parseSpecifications(mxInput, matrices)

    # Print matrix declarations
    printMatrices(matrices)


parseModel(sys.stdin.read())
