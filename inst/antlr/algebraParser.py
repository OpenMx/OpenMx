#!/usr/bin/env python
from lepl import *

class Primary(Node): pass
class Secondary(Node): pass
class Tertiary(Node): pass
class Quaternary(Node): pass
class Expression(Node): pass

secondary    = Delayed()
expression   = Delayed()
namebase     = Word(Or(Letter(), Digit()))
reference    = namebase                        > "reference"
fname        = Drop(Literal("\\")) & namebase  > "function"
spaces       = Drop(Regexp(r"\s*"))
unaryprefix  = Any("+-")                       > "operator"
unarypostfix = Any("~'")                       > "operator"
exponent     = Literal("^")                    > "operator"
tightoprs    = Any('*.@&%')                    > "operator"
looseoprs    = Any('|_+-')                     > "operator"

with Separator(spaces):
    parens      = "(" & expression & ")"
    arguments   = expression & Star(Drop(Literal(",")) & expression)
    function    = fname & Drop(Literal("(")) & arguments & Drop(Literal(")"))
    atom        = reference | parens | function
    primary     = Star(unaryprefix) & atom & Star(unarypostfix)      > Primary
    secondary  += primary & (exponent & secondary)[:1]               > Secondary
    tertiary    = secondary & Star(tightoprs & secondary)            > Tertiary
    quaternary  = tertiary & Star(looseoprs & tertiary)              > Quaternary
    expression += quaternary                                         > Expression
	

#expression :    quaternary ;	
#quaternary :    tertiary ( PIPE^ tertiary | UNDERSCORE^ tertiary | PLUS^ tertiary | MINUS^ tertiary )* ;
#tertiary   :    secondary ( ASTERISK^ secondary | DOT^ secondary | AT^ secondary | AMPERSAND^ secondary | PERCENT^ secondary )* ;
#secondary  :    primary (CARET^ secondary)? ;

#primary	   :    (MINUS^)* atom ( TILDE^ | SQUOTE^ )* ;

#atom       :    REFERENCE^ | parens | function ;
#parens     :    LPAREN^ expression RPAREN!;
	
#function   :    FNAME^ LPAREN! arguments RPAREN! ;
#arguments  :    expression (COMMA! expression)* ;
