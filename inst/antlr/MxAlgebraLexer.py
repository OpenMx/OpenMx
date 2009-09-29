# $ANTLR 3.1.2 MxAlgebra.g 2009-09-29 15:16:15

import sys
from antlr3 import *
from antlr3.compat import set, frozenset


# for convenience in actions
HIDDEN = BaseRecognizer.HIDDEN

# token types
FNAME=19
PERCENT=12
AMPERSAND=11
REFERENCE=16
UNDERSCORE=5
MINUS=7
SQUOTE=15
ID=21
EOF=-1
ASTERISK=8
LPAREN=17
AT=10
RPAREN=18
WS=22
COMMA=20
CARET=13
TILDE=14
PIPE=4
PLUS=6
DOT=9


class MxAlgebraLexer(Lexer):

    grammarFileName = "MxAlgebra.g"
    antlr_version = version_str_to_tuple("3.1.2")
    antlr_version_str = "3.1.2"

    def __init__(self, input=None, state=None):
        if state is None:
            state = RecognizerSharedState()
        Lexer.__init__(self, input, state)






    # $ANTLR start "AMPERSAND"
    def mAMPERSAND(self, ):

        try:
            _type = AMPERSAND
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:23:12: ( '&' )
            # MxAlgebra.g:23:14: '&'
            pass 
            self.match(38)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "AMPERSAND"



    # $ANTLR start "ASTERISK"
    def mASTERISK(self, ):

        try:
            _type = ASTERISK
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:24:12: ( '*' )
            # MxAlgebra.g:24:14: '*'
            pass 
            self.match(42)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "ASTERISK"



    # $ANTLR start "AT"
    def mAT(self, ):

        try:
            _type = AT
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:25:12: ( '@' )
            # MxAlgebra.g:25:14: '@'
            pass 
            self.match(64)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "AT"



    # $ANTLR start "CARET"
    def mCARET(self, ):

        try:
            _type = CARET
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:26:12: ( '^' )
            # MxAlgebra.g:26:14: '^'
            pass 
            self.match(94)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "CARET"



    # $ANTLR start "COMMA"
    def mCOMMA(self, ):

        try:
            _type = COMMA
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:27:12: ( ',' )
            # MxAlgebra.g:27:14: ','
            pass 
            self.match(44)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "COMMA"



    # $ANTLR start "DOT"
    def mDOT(self, ):

        try:
            _type = DOT
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:28:12: ( '.' )
            # MxAlgebra.g:28:14: '.'
            pass 
            self.match(46)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "DOT"



    # $ANTLR start "LPAREN"
    def mLPAREN(self, ):

        try:
            _type = LPAREN
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:29:12: ( '(' )
            # MxAlgebra.g:29:14: '('
            pass 
            self.match(40)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "LPAREN"



    # $ANTLR start "MINUS"
    def mMINUS(self, ):

        try:
            _type = MINUS
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:30:12: ( '-' )
            # MxAlgebra.g:30:14: '-'
            pass 
            self.match(45)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "MINUS"



    # $ANTLR start "PERCENT"
    def mPERCENT(self, ):

        try:
            _type = PERCENT
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:31:12: ( '%' )
            # MxAlgebra.g:31:14: '%'
            pass 
            self.match(37)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "PERCENT"



    # $ANTLR start "PIPE"
    def mPIPE(self, ):

        try:
            _type = PIPE
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:32:12: ( '|' )
            # MxAlgebra.g:32:14: '|'
            pass 
            self.match(124)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "PIPE"



    # $ANTLR start "PLUS"
    def mPLUS(self, ):

        try:
            _type = PLUS
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:33:12: ( '+' )
            # MxAlgebra.g:33:14: '+'
            pass 
            self.match(43)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "PLUS"



    # $ANTLR start "RPAREN"
    def mRPAREN(self, ):

        try:
            _type = RPAREN
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:34:12: ( ')' )
            # MxAlgebra.g:34:14: ')'
            pass 
            self.match(41)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "RPAREN"



    # $ANTLR start "SQUOTE"
    def mSQUOTE(self, ):

        try:
            _type = SQUOTE
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:35:12: ( '\\'' )
            # MxAlgebra.g:35:14: '\\''
            pass 
            self.match(39)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "SQUOTE"



    # $ANTLR start "TILDE"
    def mTILDE(self, ):

        try:
            _type = TILDE
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:36:12: ( '~' )
            # MxAlgebra.g:36:14: '~'
            pass 
            self.match(126)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "TILDE"



    # $ANTLR start "UNDERSCORE"
    def mUNDERSCORE(self, ):

        try:
            _type = UNDERSCORE
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:37:12: ( '_' )
            # MxAlgebra.g:37:14: '_'
            pass 
            self.match(95)



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "UNDERSCORE"



    # $ANTLR start "ID"
    def mID(self, ):

        try:
            # MxAlgebra.g:39:14: ( ( 'a' .. 'z' | 'A' .. 'Z' | '0' .. '9' )+ )
            # MxAlgebra.g:39:16: ( 'a' .. 'z' | 'A' .. 'Z' | '0' .. '9' )+
            pass 
            # MxAlgebra.g:39:16: ( 'a' .. 'z' | 'A' .. 'Z' | '0' .. '9' )+
            cnt1 = 0
            while True: #loop1
                alt1 = 2
                LA1_0 = self.input.LA(1)

                if ((48 <= LA1_0 <= 57) or (65 <= LA1_0 <= 90) or (97 <= LA1_0 <= 122)) :
                    alt1 = 1


                if alt1 == 1:
                    # MxAlgebra.g:
                    pass 
                    if (48 <= self.input.LA(1) <= 57) or (65 <= self.input.LA(1) <= 90) or (97 <= self.input.LA(1) <= 122):
                        self.input.consume()
                    else:
                        mse = MismatchedSetException(None, self.input)
                        self.recover(mse)
                        raise mse



                else:
                    if cnt1 >= 1:
                        break #loop1

                    eee = EarlyExitException(1, self.input)
                    raise eee

                cnt1 += 1






        finally:

            pass

    # $ANTLR end "ID"



    # $ANTLR start "REFERENCE"
    def mREFERENCE(self, ):

        try:
            _type = REFERENCE
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:40:14: ( ID )
            # MxAlgebra.g:40:17: ID
            pass 
            self.mID()



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "REFERENCE"



    # $ANTLR start "FNAME"
    def mFNAME(self, ):

        try:
            _type = FNAME
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:41:14: ( '\\\\' ID )
            # MxAlgebra.g:41:16: '\\\\' ID
            pass 
            self.match(92)
            self.mID()



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "FNAME"



    # $ANTLR start "WS"
    def mWS(self, ):

        try:
            _type = WS
            _channel = DEFAULT_CHANNEL

            # MxAlgebra.g:43:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            # MxAlgebra.g:43:9: ( ' ' | '\\t' | '\\r' | '\\n' )
            pass 
            if (9 <= self.input.LA(1) <= 10) or self.input.LA(1) == 13 or self.input.LA(1) == 32:
                self.input.consume()
            else:
                mse = MismatchedSetException(None, self.input)
                self.recover(mse)
                raise mse

            #action start
            _channel=HIDDEN;
            #action end



            self._state.type = _type
            self._state.channel = _channel

        finally:

            pass

    # $ANTLR end "WS"



    def mTokens(self):
        # MxAlgebra.g:1:8: ( AMPERSAND | ASTERISK | AT | CARET | COMMA | DOT | LPAREN | MINUS | PERCENT | PIPE | PLUS | RPAREN | SQUOTE | TILDE | UNDERSCORE | REFERENCE | FNAME | WS )
        alt2 = 18
        LA2 = self.input.LA(1)
        if LA2 == 38:
            alt2 = 1
        elif LA2 == 42:
            alt2 = 2
        elif LA2 == 64:
            alt2 = 3
        elif LA2 == 94:
            alt2 = 4
        elif LA2 == 44:
            alt2 = 5
        elif LA2 == 46:
            alt2 = 6
        elif LA2 == 40:
            alt2 = 7
        elif LA2 == 45:
            alt2 = 8
        elif LA2 == 37:
            alt2 = 9
        elif LA2 == 124:
            alt2 = 10
        elif LA2 == 43:
            alt2 = 11
        elif LA2 == 41:
            alt2 = 12
        elif LA2 == 39:
            alt2 = 13
        elif LA2 == 126:
            alt2 = 14
        elif LA2 == 95:
            alt2 = 15
        elif LA2 == 48 or LA2 == 49 or LA2 == 50 or LA2 == 51 or LA2 == 52 or LA2 == 53 or LA2 == 54 or LA2 == 55 or LA2 == 56 or LA2 == 57 or LA2 == 65 or LA2 == 66 or LA2 == 67 or LA2 == 68 or LA2 == 69 or LA2 == 70 or LA2 == 71 or LA2 == 72 or LA2 == 73 or LA2 == 74 or LA2 == 75 or LA2 == 76 or LA2 == 77 or LA2 == 78 or LA2 == 79 or LA2 == 80 or LA2 == 81 or LA2 == 82 or LA2 == 83 or LA2 == 84 or LA2 == 85 or LA2 == 86 or LA2 == 87 or LA2 == 88 or LA2 == 89 or LA2 == 90 or LA2 == 97 or LA2 == 98 or LA2 == 99 or LA2 == 100 or LA2 == 101 or LA2 == 102 or LA2 == 103 or LA2 == 104 or LA2 == 105 or LA2 == 106 or LA2 == 107 or LA2 == 108 or LA2 == 109 or LA2 == 110 or LA2 == 111 or LA2 == 112 or LA2 == 113 or LA2 == 114 or LA2 == 115 or LA2 == 116 or LA2 == 117 or LA2 == 118 or LA2 == 119 or LA2 == 120 or LA2 == 121 or LA2 == 122:
            alt2 = 16
        elif LA2 == 92:
            alt2 = 17
        elif LA2 == 9 or LA2 == 10 or LA2 == 13 or LA2 == 32:
            alt2 = 18
        else:
            nvae = NoViableAltException("", 2, 0, self.input)

            raise nvae

        if alt2 == 1:
            # MxAlgebra.g:1:10: AMPERSAND
            pass 
            self.mAMPERSAND()


        elif alt2 == 2:
            # MxAlgebra.g:1:20: ASTERISK
            pass 
            self.mASTERISK()


        elif alt2 == 3:
            # MxAlgebra.g:1:29: AT
            pass 
            self.mAT()


        elif alt2 == 4:
            # MxAlgebra.g:1:32: CARET
            pass 
            self.mCARET()


        elif alt2 == 5:
            # MxAlgebra.g:1:38: COMMA
            pass 
            self.mCOMMA()


        elif alt2 == 6:
            # MxAlgebra.g:1:44: DOT
            pass 
            self.mDOT()


        elif alt2 == 7:
            # MxAlgebra.g:1:48: LPAREN
            pass 
            self.mLPAREN()


        elif alt2 == 8:
            # MxAlgebra.g:1:55: MINUS
            pass 
            self.mMINUS()


        elif alt2 == 9:
            # MxAlgebra.g:1:61: PERCENT
            pass 
            self.mPERCENT()


        elif alt2 == 10:
            # MxAlgebra.g:1:69: PIPE
            pass 
            self.mPIPE()


        elif alt2 == 11:
            # MxAlgebra.g:1:74: PLUS
            pass 
            self.mPLUS()


        elif alt2 == 12:
            # MxAlgebra.g:1:79: RPAREN
            pass 
            self.mRPAREN()


        elif alt2 == 13:
            # MxAlgebra.g:1:86: SQUOTE
            pass 
            self.mSQUOTE()


        elif alt2 == 14:
            # MxAlgebra.g:1:93: TILDE
            pass 
            self.mTILDE()


        elif alt2 == 15:
            # MxAlgebra.g:1:99: UNDERSCORE
            pass 
            self.mUNDERSCORE()


        elif alt2 == 16:
            # MxAlgebra.g:1:110: REFERENCE
            pass 
            self.mREFERENCE()


        elif alt2 == 17:
            # MxAlgebra.g:1:120: FNAME
            pass 
            self.mFNAME()


        elif alt2 == 18:
            # MxAlgebra.g:1:126: WS
            pass 
            self.mWS()







 



def main(argv, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr):
    from antlr3.main import LexerMain
    main = LexerMain(MxAlgebraLexer)
    main.stdin = stdin
    main.stdout = stdout
    main.stderr = stderr
    main.execute(argv)


if __name__ == '__main__':
    main(sys.argv)
