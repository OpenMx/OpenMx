#!/usr/bin/env python

# Mx 1.0 algebras are parsed using the
# PLY implementation of lex and yacc parsing tools for Python.

cognates = ['det', 'tr', 'sum', 'prod', 'max', \
            'min', 'abs', 'cos', 'cosh', 'sin', \
            'sinh', 'tan', 'tanh', 'exp', 'sqrt']

tokens = (
    'NAME', 'FNAME',
    'AMPERSAND', 'ASTERISK', 'AT', 'CARET',
    'DOT', 'MINUS', 'PERCENT', 'PIPE',
    'PLUS', 'SQUOTE', 'TILDE', 'UNDERSCORE',
    'LPAREN','RPAREN', 'COMMA'
    )

# Tokens

t_AMPERSAND  = r'&'
t_ASTERISK   = r'\*'
t_AT         = r'@'
t_CARET      = r'\^'
t_DOT        = r'\.'
t_MINUS      = r'-'
t_PERCENT    = r'%'
t_PIPE       = r'\|'
t_PLUS       = r'\+'
t_SQUOTE     = r'\''
t_TILDE      = r'~'
t_UNDERSCORE = r'_'
t_LPAREN     = r'\('
t_RPAREN     = r'\)'
t_COMMA      = r','
t_FNAME      = r'\\[a-zA-Z0-9]+'
t_NAME       = r'[a-zA-Z0-9]+'

# Ignored characters
t_ignore = " \t\n"
    
def t_error(t):
    raise Exception("Illegal character " + str(t.value[0]))
    
# Build the lexer
import ply.lex as lex
lex.lex()

# Parsing rules

precedence = (
    ('left','UNDERSCORE'),
    ('left','PIPE'),
    ('left','PLUS','MINUS'),
    ('left','ASTERISK','DOT','AT','AMPERSAND','PERCENT'),
    ('right','CARET'),
    ('right','UMINUS'),    
    ('left','TILDE','SQUOTE')
    )

def p_statement_expr(t):
    'statement : expression'
    t[0] = t[1]

def p_expression_binop(t):
    '''expression : expression PLUS expression
                  | expression MINUS expression
                  | expression PIPE expression
                  | expression UNDERSCORE expression
                  | expression ASTERISK expression
                  | expression DOT expression
                  | expression AT expression                  
                  | expression AMPERSAND expression
                  | expression PERCENT expression
                  | expression CARET expression'''
    if t[2] == '+'  : t[0] = t[1] + ' + ' + t[3]
    elif t[2] == '-': t[0] = t[1] + ' - ' + t[3]
    elif t[2] == '|': t[0] = 'cbind(' + t[1] + ', ' + t[3] + ')'
    elif t[2] == '_': t[0] = 'rbind(' + t[1] + ', ' + t[3] + ')'
    elif t[2] == '*': t[0] = t[1] + ' %*% ' + t[3]  
    elif t[2] == '.': t[0] = '(' + t[1] + ' * ' + t[3] + ')'
    elif t[2] == '@': t[0] = t[1] + ' %x% ' + t[3]    
    elif t[2] == '&': t[0] = t[1] + ' %&% ' + t[3]
    elif t[2] == '%': t[0] = '(' + t[1] + ' / ' + t[3] + ')'    
    elif t[2] == '^': t[0] = t[1] + ' ^ ' + t[3]
    

def p_expression_unaryop(t):
    '''expression : expression TILDE
                  | expression SQUOTE'''
    if t[2] == '~' : t[0] = 'solve(' + t[1] + ')'
    elif t[2] == '\'' : t[0] = 't(' + t[1] + ')'

def p_expression_uminus(t):
    'expression : MINUS expression %prec UMINUS'
    t[0] = "(-" + t[2] + ")"

def p_expression_group(t):
    'expression : LPAREN expression RPAREN'
    t[0] = "(" + t[2] + ")"

def p_expression_function(t):
    'expression : FNAME LPAREN arglist RPAREN'
    funcname = convertFunctionName(t[1][1:])
    t[0] = funcname + '('
    for i in range(len(t[3])):
        t[0] += t[3][i]
        if i < len(t[3]) - 1: t[0] += ', '
    t[0] += ')'

def convertFunctionName(name):
    if name in cognates:
        return(name)
    elif fname == 'ln':
        return('log')
    else:
        raise Exception("Function has not been implemented: " + fname)    
    

def p_arglist(t):
    '''arglist : expression args
               | empty'''
    if len(t) == 3: t[0] = [t[1]] + t[2]
    else: t[0] = []

def p_args(t):
    '''args : COMMA expression args
            | empty'''
    if len(t) == 4: t[0] = [t[2]] + t[3]
    else: t[0] = []
    
def p_empty(t):
    'empty :'
    t[0] = []

def p_expression_name(t):
    'expression : NAME'
    t[0] = t[1]

def p_error(t):
    raise Exception("Syntax error on token " + str(t.value) +\
                    " at line " + str(t.lineno) + " at position " +\
                    str(t.lexpos))

import ply.yacc as parser
parser.yacc(write_tables=0,debug=0)

if __name__ == "__main__":
    import sys
    lines = sys.stdin.readlines()
    input = "".join(lines)
    print(parser.parse(input))
