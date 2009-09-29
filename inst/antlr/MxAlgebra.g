grammar MxAlgebra;
options {
     language=Python;
}

// BEGIN MATRIX ALGEBRA FORMULAS
// Comment -- I am assuming all operators are left-associative,
// except for exponentiation which is right-associative. 

expression :    quaternary ;	
quaternary :    tertiary ( PIPE^ tertiary | UNDERSCORE^ tertiary | PLUS^ tertiary | MINUS^ tertiary )* ;
tertiary   :    secondary ( ASTERISK^ secondary | DOT^ secondary | AT^ secondary | AMPERSAND^ secondary | PERCENT^ secondary )* ;
secondary  :    primary (CARET^ secondary)? ;

primary	   :    (MINUS^)* atom ( TILDE^ | SQUOTE^ )* ;

atom       :    REFERENCE^ | parens | function ;
parens     :    LPAREN^ expression RPAREN!;
	
function   :    FNAME^ LPAREN! arguments RPAREN! ;
arguments  :    expression (COMMA! expression)* ;

AMPERSAND  : '&' ;
ASTERISK   : '*' ;
AT         : '@' ;
CARET      : '^' ;
COMMA      : ',' ;
DOT        : '.' ;
LPAREN     : '(' ;
MINUS      : '-' ;
PERCENT    : '%' ;
PIPE       : '|' ;	
PLUS       : '+' ;
RPAREN     : ')' ;
SQUOTE 	   : '\'';             // single quote
TILDE      : '~' ;
UNDERSCORE : '_' ;

fragment ID  :	('a'..'z'|'A'..'Z'|'0'..'9')+ ;
REFERENCE    :	 ID ;
FNAME        :	'\\' ID ;
	
WS  :   ( ' '
        | '\t'
        | '\r'
        | '\n'
        ) {$channel=HIDDEN;}
    ;