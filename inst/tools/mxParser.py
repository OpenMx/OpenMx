#!/usr/bin/env python
import sys
import re
import string

matrices = {}
defines = {}
title = None
valBuffer = ""

class MxMatrix:
    name = None
    type = None
    free = False
    unique = False
    nrow = None
    ncol = None
    values = None
    specification = None

def parseDefine( mxInput ):
    global defines
    match = re.match("\s*#define\s+(\S+)\s+(\S+)\s*", mxInput)
    defines[match.group(1)] = match.group(2)
    return len(match.group(0))

def parseTitle( mxInput ):
    global title
    match = re.match("\s*Title(.*)$", mxInput, re.MULTILINE | re.IGNORECASE)
    title = match.group(1).strip()
    title = title.replace('.','_')
    return len(match.group(0)) + 1

def parseMatrices( mxInput ):
    global defines, matrices
    block = re.match("\s*Begin Matrices;(.*?)End Matrices;\s*", 
    	mxInput, re.MULTILINE | re.IGNORECASE | re.DOTALL)
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
            matrices[matrix.name] = matrix
    return len(block.group(0))

def parseSpecification( mxInput ):
    global defines, matrices
    match = re.match("\s*spec\S*\s+(\S+)", mxInput, re.IGNORECASE | re.MULTILINE)
    specification = list()
    matrixName = match.group(1)
    matrix = matrices[matrixName]
    specCount = specLength[matrix.type](matrix.nrow, matrix.ncol)
    specs = re.search("\s*spec\S*\s+" + matrixName + "((\s*\S+\s*)" + "{" + str(specCount) + "})\s*",
                       mxInput, re.IGNORECASE | re.MULTILINE)
    specIter = re.finditer("\S+", specs.group(1))
    for spec in specIter:
         try:
             specification.append(int(spec.group(0)))
         except:
             specification.append(int(defines[spec.group(0)]))
    matrices[matrixName].specification = specification
    return len(specs.group(0))

def parseStartOrValue( mxInput ):
    global valBuffer
    match = re.match("\s*(Start|Value)\s+(\S+)\s+(All|(\S+[ ]+)+)\s*", mxInput, re.IGNORECASE)
    if (match.group(3) == "All"):
        for matrix in matrices.values():
            mname = "matrix" + matrix.name
            if match.group(2) == "Value":
                direction = "!"
            else:
                direction = ""
            valBuffer += mname + "@values[" + direction + mname + "@free] <- " + match.group(2) + '\n'
    else:
       tripleIter = re.finditer("(\S+)\s*(\S+)\s*(\S+)\s*", match.group(3))
       for triple in tripleIter:
           mname = "matrix" + triple.group(1)
           if triple.group(2) in defines:
               row = defines[triple.group(2)]
           else:
               row = triple.group(2)
           if triple.group(3) in defines:
               col = defines[triple.group(3)]
           else:
               col = triple.group(3)
           valBuffer += mname + "@values[" + row + "," + col + "] <- " + match.group(2) + '\n'
    valBuffer += '\n'
    return len(match.group(0))
        

def printMatrices():
    global matrices
    startCounter = 1
    for matrix in matrices.values():
        outstring = "matrix" + matrix.name + " <- "
        outstring += "mxMatrix(type = \"" + matrix.type + "\", "
        outstring += "nrow = " + str(matrix.nrow) + ", "
        outstring += "ncol = " + str(matrix.ncol) + ", "
        if matrix.specification == None:
            outstring += "free = " + str(matrix.free).upper() + ", "
        else:
            freelist = map(lambda x: x > 0, matrix.specification)
            freestring = str(freelist).upper()
            freestring = string.replace(freestring, "[", "c(")
            freestring = string.replace(freestring, "]", ")")
            outstring += "free = " + freestring + ", "
        nextcounter = None
        if matrix.free and not matrix.unique and matrix.specification == None:
            endCounter = startCounter + specLength[matrix.type](matrix.nrow, matrix.ncol) - 1
            outstring += "labels = " + str(startCounter) + " : "
            outstring += str(endCounter) + ", "
            startCounter = endCounter + 1
        elif matrix.specification != None:
            specstring = str(matrix.specification)
            specstring = string.replace(specstring, "[", "c(")
            specstring = string.replace(specstring, "]", ")")
            outstring += "labels = " + specstring + ", "
        outstring += "byrow = TRUE, "
        outstring += "name = \"" + matrix.name + "\")"
        print outstring
        print

def printValueBuffer():
    global valBuffer
    print valBuffer

def printModel():
    global matrices
    print "model <- mxModel(name = \"" + title + "\")"
    matrixNames = str(map(lambda x: "matrix" + x, matrices.keys()))
    matrixNames = string.replace(matrixNames, "[", "")
    matrixNames = string.replace(matrixNames, "]", "")
    matrixNames = string.replace(matrixNames, "'", "")
    print "model <- mxModel(model, " + matrixNames + ")"
    print

specLength = {   "Diag"  : lambda row, col: row,
                 "SDiag" : lambda row, col: row * (row - 1) / 2,
                 "Stand" : lambda row, col: row * (row - 1) / 2,
                 "Symm"  : lambda row, col: row * (row + 1) / 2,
                 "Lower" : lambda row, col: row * (row + 1) / 2,
                 "Full"  : lambda row, col: row * col }

mxDirectives = {    "\s*title" : parseTitle,
                    "\s*#define" : parseDefine,
                    "\s*Begin Matrices" : parseMatrices,
                    "\s*spec\S*" : parseSpecification,
                    "\s*(start|value)" : parseStartOrValue }

def tryDirectives ( mxInput ):
    for directive in mxDirectives.keys():
        match = re.match(directive, mxInput, re.IGNORECASE)
        if (match != None):
            subtract = mxDirectives[directive](mxInput)
            mxInput = mxInput[subtract : ]
            return mxInput
    match = re.match(".+", mxInput)
    if (match != None):
        subtract = len(match.group(0)) + 1
        mxInput = mxInput[subtract : ]
    return mxInput
                
def parseModel( mxInput ):

    # Remove all comments from the file
    mxInput = re.sub(re.compile('!.*$', re.MULTILINE), '', mxInput)
    
    while len(mxInput) > 0:
        mxInput = tryDirectives(mxInput)

    # Print matrix declarations
    printMatrices()

    # Print any value assignments
    printValueBuffer()

    # Print model declaration
    printModel()

parseModel(sys.stdin.read())
