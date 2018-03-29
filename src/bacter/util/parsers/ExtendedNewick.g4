grammar ExtendedNewick;


// Parser rules:

tree: node ';'? EOF;

node: ('(' node (',' node)* ')')? post ;

post: label? hybrid? meta? (':' length=number)? ;

label: number | string ;

hybrid: '#' type=RECTYPE? INT ;

meta: '[&' attrib (',' attrib)* ']' ;

attrib: attribKey=string '=' attribValue ;

attribValue: number | string | vector;


number: INT | FLOAT ;

vector: '{' attribValue (',' attribValue)* '}' ;

string : RECTYPE | STRINGPRIM ;


// Lexer rules:

FLOAT : '-'? NNINT ('.' D*) ([eE] ('-'|'+')? D+)? ;
INT : '-'? NNINT;
fragment NNINT : '0' | NZD D* ;
fragment NZD : [1-9] ;
fragment D : [0-9] ;

RECTYPE: 'R' | 'H' | 'LGT' ;

STRINGPRIM :
    [a-zA-Z0-9|*%/.\-+_&]+  // these chars don't need quotes
    | '"' .*? '"'
    | '\'' .*? '\''
    ;

WHITESPACE : [ \t\r\n]+ -> skip ;