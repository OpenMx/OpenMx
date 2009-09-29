import os
import sys
import antlr3

sys.path.append(os.path.abspath(os.path.dirname(sys.argv[0])))

import MxAlgebraParser
import MxAlgebraLexer

def generateAlgebraTree( string ):
    char_stream = antlr3.ANTLRStringStream(string)
    lexer = MxAlgebraLexer.MxAlgebraLexer(char_stream)
    tokens = antlr3.CommonTokenStream(lexer)
    parser = MxAlgebraParser.MxAlgebraParser(tokens)
    return(parser.expression())

def antlrParseAlgebra( string ):
    ast = generateAlgebraTree(string)
    return parseHelper(ast.tree)

operators = [MxAlgebraLexer.AMPERSAND, MxAlgebraLexer.ASTERISK, \
             MxAlgebraLexer.AT, MxAlgebraLexer.CARET, \
             MxAlgebraLexer.DOT, MxAlgebraLexer.LPAREN, \
             MxAlgebraLexer.MINUS, MxAlgebraLexer.PERCENT, \
             MxAlgebraLexer.PIPE, MxAlgebraLexer.PLUS, \
             MxAlgebraLexer.SQUOTE, MxAlgebraLexer.TILDE, \
             MxAlgebraLexer.UNDERSCORE]

cognates = ['det', 'tr', 'sum', 'prod', 'max', \
            'min', 'abs', 'cos', 'cosh', 'sin', \
            'sinh', 'tan', 'tanh', 'exp', 'sqrt']

def parseHelper( tree ):
    head = tree.getType()
    if head in operators:
        return(parseOperator(tree))
    elif head == MxAlgebraLexer.FNAME:
        return(parseFunction(tree))
    elif head == MxAlgebraLexer.REFERENCE:
        return(tree.getText())
    else:
        raise Exception("Error parsing algebra expression, text is '" + \
                        tree.getText() + "' and number of children is " + \
                        str(tree.getChildCount()))

def parseOperator( tree ):
    head = tree.getType()
    if head == MxAlgebraLexer.AMPERSAND:
        return(parseHelper(tree.getChild(0)) + \
               ' %&% ' + parseHelper(tree.getChild(1)))    
    elif head == MxAlgebraLexer.ASTERISK:
        return(parseHelper(tree.getChild(0)) + \
               ' %*% ' + parseHelper(tree.getChild(1)))
    elif head == MxAlgebraLexer.AT:
        return(parseHelper(tree.getChild(0)) + \
               ' %x% ' + parseHelper(tree.getChild(1)))
    elif head == MxAlgebraLexer.CARET:
        return(parseHelper(tree.getChild(0)) + \
               ' ^ ' + parseHelper(tree.getChild(1)))
    elif head == MxAlgebraLexer.DOT:
        return(parseHelper(tree.getChild(0)) + \
               ' * ' + parseHelper(tree.getChild(1)))
    elif head == MxAlgebraLexer.LPAREN:
        return('(' + parseHelper(tree.getChild(0)) + ')')    
    elif head == MxAlgebraLexer.MINUS and tree.getChildCount() == 1:
        return('-' + parseHelper(tree.getChild(0)))
    elif head == MxAlgebraLexer.MINUS and tree.getChildCount() == 2:
        return(parseHelper(tree.getChild(0)) + \
               ' - ' + parseHelper(tree.getChild(1)))
    elif head == MxAlgebraLexer.PERCENT:
        return(parseHelper(tree.getChild(0)) + \
               ' / ' + parseHelper(tree.getChild(1)))
    elif head == MxAlgebraLexer.PIPE:
        return('cbind(' + parseHelper(tree.getChild(0)) + \
               ', ' + parseHelper(tree.getChild(1)) + ')')
    elif head == MxAlgebraLexer.PLUS:
        return(parseHelper(tree.getChild(0)) + \
               ' + ' + parseHelper(tree.getChild(1)))
    elif head == MxAlgebraLexer.SQUOTE:
        return('t(' + parseHelper(tree.getChild(0)) + ')')
    elif head == MxAlgebraLexer.TILDE:
        return('solve(' + parseHelper(tree.getChild(0)) + ')')    
    elif head == MxAlgebraLexer.UNDERSCORE:
        return('rbind(' + parseHelper(tree.getChild(0)) + \
               ', ' + parseHelper(tree.getChild(1)) + ')')
    elif head == MxAlgebraLexer.REFERENCE:
        return(tree.getText())
    else:
        raise Exception("Operator has not been implemented: " + tree.getText())

def parseFunction( tree ):
    fname = tree.getText()[1:]
    if fname in cognates:
        return(parseFunctionHelper(tree, fname))
    elif fname == 'ln':
        return(parseFunctionHelper(tree, 'log'))
    else:
        raise Exception("Function has not been implemented: " + fname)

def parseFunctionHelper( tree, fname ):
    result = fname + '('
    count = tree.getChildCount()
    for i in range(count):
        result += parseHelper(tree.getChild(i))
        if i < (count - 1):
            result += ', '
    result += ')'
    return(result)
