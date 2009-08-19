#!/usr/bin/env python
import sys
import re
import string

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

def parseSpecifications( mxInput, defines, matrices ):
    matchIter = re.finditer("^spec\S*\s+(\S+)", mxInput, re.IGNORECASE | re.MULTILINE)
    for match in matchIter:
        specification = list()
        matrixName = match.group(1)
        matrix = matrices[matrixName]
        specCount = specLength[matrix.type](matrix.nrow, matrix.ncol)
        specs = re.search("^spec\S*\s+" + matrixName + "((\s*\S+\s*)" + "{" + str(specCount) + "})",
                          mxInput, re.IGNORECASE | re.MULTILINE)
        specIter = re.finditer("\S+", specs.group(1))
        for spec in specIter:
            try:
                specification.append(int(spec.group(0)))
            except:
                specification.append(int(defines[spec.group(0)]))
        matrices[matrixName].specification = specification

        

def printMatrices( matrices ):
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

def printModel( title, matrices ):
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

def parseStartOrValue( mxInput, defines, matrices ):
    matchIter = re.finditer("(Start|Value)\s+(\S+)\s+(All|(\S+[ ]+)+)", mxInput, re.IGNORECASE)
    for match in matchIter:
        if (match.group(3) == "All"):
            for matrix in matrices.values():
                mname = "matrix" + matrix.name
                if match.group(2) == "Value":
                    direction = "!"
                else:
                    direction = ""
                print mname + "@values[" + direction + mname + "@free] <- " + match.group(2)
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
               print mname + "@values[" + row + "," + col + "] <- " + match.group(2)
        print    

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
    parseSpecifications(mxInput, defines, matrices)

    # Print matrix declarations
    printMatrices(matrices)

    #Find all Start or Value declarations and print them
    parseStartOrValue(mxInput, defines, matrices)


    # Print model declaration
    printModel(title, matrices)

parseModel(sys.stdin.read())
