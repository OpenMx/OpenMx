grammar MxAlgebra;

// BEGIN MATRIX ALGEBRA FORMULAS
// Comment -- I am assuming all operators are left-associative,
// except for exponentiation which is right-associative. 

expression :	quaternary ;
	
quaternary :	tertiary ( hAdhere | vAdhere  | plus | minus )* ;
hAdhere	   :	PIPE! tertiary ;
vAdhere	   :	UNDERSCORE! tertiary ;
plus	   :	PLUS! tertiary ;
minus	   :	MINUS! tertiary ;

tertiary   :	secondary ( mult | dot | kronecker | quadratic | div )* ;
mult	   :	ASTERISK! secondary ;
dot 	   :	DOT! secondary ;
kronecker  :	AT! secondary ;
quadratic  :	AMPERSAND! secondary ;
div	   :	PERCENT! secondary ;

exponent   :	CARET! ;
secondary  :	primary (exponent secondary)? ;

primary	   :	mAtom ( inverse | transpose )* ;
inverse	   :	TILDE! ;
transpose  :	SQUOTE! ;
	
mAtom	   :	REFERENCE | parens | function ;
parens     :	LPAREN! expression RPAREN! ;
	
function   :	FNAME LPAREN! arguments RPAREN! ;
arguments  :	expression (COMMA! expression)? ;	

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
PLUS 	   : '+' ;
RPAREN	   : ')' ;
SQUOTE 	   : '\'';             // single quote
TILDE      : '~' ;
UNDERSCORE : '_' ;

fragment ID  :	('a'..'z'|'A'..'Z'|'0'..'9')+ ;
REFERENCE    :	 ID ;
FNAME	     :	'\\' ID ;
	
WS  :   ( ' '
        | '\t'
        | '\r'
        | '\n'
        ) {$channel=HIDDEN;}
    ;

