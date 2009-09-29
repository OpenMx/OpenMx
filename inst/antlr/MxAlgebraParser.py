# $ANTLR 3.1.2 MxAlgebra.g 2009-09-29 15:16:14

import sys
from antlr3 import *
from antlr3.compat import set, frozenset

from antlr3.tree import *



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
PLUS=6
PIPE=4
DOT=9

# token names
tokenNames = [
    "<invalid>", "<EOR>", "<DOWN>", "<UP>", 
    "PIPE", "UNDERSCORE", "PLUS", "MINUS", "ASTERISK", "DOT", "AT", "AMPERSAND", 
    "PERCENT", "CARET", "TILDE", "SQUOTE", "REFERENCE", "LPAREN", "RPAREN", 
    "FNAME", "COMMA", "ID", "WS"
]




class MxAlgebraParser(Parser):
    grammarFileName = "MxAlgebra.g"
    antlr_version = version_str_to_tuple("3.1.2")
    antlr_version_str = "3.1.2"
    tokenNames = tokenNames

    def __init__(self, input, state=None):
        if state is None:
            state = RecognizerSharedState()

        Parser.__init__(self, input, state)







                
        self._adaptor = CommonTreeAdaptor()


        
    def getTreeAdaptor(self):
        return self._adaptor

    def setTreeAdaptor(self, adaptor):
        self._adaptor = adaptor

    adaptor = property(getTreeAdaptor, setTreeAdaptor)


    class expression_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "expression"
    # MxAlgebra.g:10:1: expression : quaternary ;
    def expression(self, ):

        retval = self.expression_return()
        retval.start = self.input.LT(1)

        root_0 = None

        quaternary1 = None



        try:
            try:
                # MxAlgebra.g:10:12: ( quaternary )
                # MxAlgebra.g:10:17: quaternary
                pass 
                root_0 = self._adaptor.nil()

                self._state.following.append(self.FOLLOW_quaternary_in_expression31)
                quaternary1 = self.quaternary()

                self._state.following.pop()
                self._adaptor.addChild(root_0, quaternary1.tree)



                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "expression"

    class quaternary_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "quaternary"
    # MxAlgebra.g:11:1: quaternary : tertiary ( PIPE tertiary | UNDERSCORE tertiary | PLUS tertiary | MINUS tertiary )* ;
    def quaternary(self, ):

        retval = self.quaternary_return()
        retval.start = self.input.LT(1)

        root_0 = None

        PIPE3 = None
        UNDERSCORE5 = None
        PLUS7 = None
        MINUS9 = None
        tertiary2 = None

        tertiary4 = None

        tertiary6 = None

        tertiary8 = None

        tertiary10 = None


        PIPE3_tree = None
        UNDERSCORE5_tree = None
        PLUS7_tree = None
        MINUS9_tree = None

        try:
            try:
                # MxAlgebra.g:11:12: ( tertiary ( PIPE tertiary | UNDERSCORE tertiary | PLUS tertiary | MINUS tertiary )* )
                # MxAlgebra.g:11:17: tertiary ( PIPE tertiary | UNDERSCORE tertiary | PLUS tertiary | MINUS tertiary )*
                pass 
                root_0 = self._adaptor.nil()

                self._state.following.append(self.FOLLOW_tertiary_in_quaternary43)
                tertiary2 = self.tertiary()

                self._state.following.pop()
                self._adaptor.addChild(root_0, tertiary2.tree)
                # MxAlgebra.g:11:26: ( PIPE tertiary | UNDERSCORE tertiary | PLUS tertiary | MINUS tertiary )*
                while True: #loop1
                    alt1 = 5
                    LA1 = self.input.LA(1)
                    if LA1 == PIPE:
                        alt1 = 1
                    elif LA1 == UNDERSCORE:
                        alt1 = 2
                    elif LA1 == PLUS:
                        alt1 = 3
                    elif LA1 == MINUS:
                        alt1 = 4

                    if alt1 == 1:
                        # MxAlgebra.g:11:28: PIPE tertiary
                        pass 
                        PIPE3=self.match(self.input, PIPE, self.FOLLOW_PIPE_in_quaternary47)

                        PIPE3_tree = self._adaptor.createWithPayload(PIPE3)
                        root_0 = self._adaptor.becomeRoot(PIPE3_tree, root_0)

                        self._state.following.append(self.FOLLOW_tertiary_in_quaternary50)
                        tertiary4 = self.tertiary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, tertiary4.tree)


                    elif alt1 == 2:
                        # MxAlgebra.g:11:45: UNDERSCORE tertiary
                        pass 
                        UNDERSCORE5=self.match(self.input, UNDERSCORE, self.FOLLOW_UNDERSCORE_in_quaternary54)

                        UNDERSCORE5_tree = self._adaptor.createWithPayload(UNDERSCORE5)
                        root_0 = self._adaptor.becomeRoot(UNDERSCORE5_tree, root_0)

                        self._state.following.append(self.FOLLOW_tertiary_in_quaternary57)
                        tertiary6 = self.tertiary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, tertiary6.tree)


                    elif alt1 == 3:
                        # MxAlgebra.g:11:68: PLUS tertiary
                        pass 
                        PLUS7=self.match(self.input, PLUS, self.FOLLOW_PLUS_in_quaternary61)

                        PLUS7_tree = self._adaptor.createWithPayload(PLUS7)
                        root_0 = self._adaptor.becomeRoot(PLUS7_tree, root_0)

                        self._state.following.append(self.FOLLOW_tertiary_in_quaternary64)
                        tertiary8 = self.tertiary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, tertiary8.tree)


                    elif alt1 == 4:
                        # MxAlgebra.g:11:85: MINUS tertiary
                        pass 
                        MINUS9=self.match(self.input, MINUS, self.FOLLOW_MINUS_in_quaternary68)

                        MINUS9_tree = self._adaptor.createWithPayload(MINUS9)
                        root_0 = self._adaptor.becomeRoot(MINUS9_tree, root_0)

                        self._state.following.append(self.FOLLOW_tertiary_in_quaternary71)
                        tertiary10 = self.tertiary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, tertiary10.tree)


                    else:
                        break #loop1





                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "quaternary"

    class tertiary_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "tertiary"
    # MxAlgebra.g:12:1: tertiary : secondary ( ASTERISK secondary | DOT secondary | AT secondary | AMPERSAND secondary | PERCENT secondary )* ;
    def tertiary(self, ):

        retval = self.tertiary_return()
        retval.start = self.input.LT(1)

        root_0 = None

        ASTERISK12 = None
        DOT14 = None
        AT16 = None
        AMPERSAND18 = None
        PERCENT20 = None
        secondary11 = None

        secondary13 = None

        secondary15 = None

        secondary17 = None

        secondary19 = None

        secondary21 = None


        ASTERISK12_tree = None
        DOT14_tree = None
        AT16_tree = None
        AMPERSAND18_tree = None
        PERCENT20_tree = None

        try:
            try:
                # MxAlgebra.g:12:12: ( secondary ( ASTERISK secondary | DOT secondary | AT secondary | AMPERSAND secondary | PERCENT secondary )* )
                # MxAlgebra.g:12:17: secondary ( ASTERISK secondary | DOT secondary | AT secondary | AMPERSAND secondary | PERCENT secondary )*
                pass 
                root_0 = self._adaptor.nil()

                self._state.following.append(self.FOLLOW_secondary_in_tertiary87)
                secondary11 = self.secondary()

                self._state.following.pop()
                self._adaptor.addChild(root_0, secondary11.tree)
                # MxAlgebra.g:12:27: ( ASTERISK secondary | DOT secondary | AT secondary | AMPERSAND secondary | PERCENT secondary )*
                while True: #loop2
                    alt2 = 6
                    LA2 = self.input.LA(1)
                    if LA2 == ASTERISK:
                        alt2 = 1
                    elif LA2 == DOT:
                        alt2 = 2
                    elif LA2 == AT:
                        alt2 = 3
                    elif LA2 == AMPERSAND:
                        alt2 = 4
                    elif LA2 == PERCENT:
                        alt2 = 5

                    if alt2 == 1:
                        # MxAlgebra.g:12:29: ASTERISK secondary
                        pass 
                        ASTERISK12=self.match(self.input, ASTERISK, self.FOLLOW_ASTERISK_in_tertiary91)

                        ASTERISK12_tree = self._adaptor.createWithPayload(ASTERISK12)
                        root_0 = self._adaptor.becomeRoot(ASTERISK12_tree, root_0)

                        self._state.following.append(self.FOLLOW_secondary_in_tertiary94)
                        secondary13 = self.secondary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, secondary13.tree)


                    elif alt2 == 2:
                        # MxAlgebra.g:12:51: DOT secondary
                        pass 
                        DOT14=self.match(self.input, DOT, self.FOLLOW_DOT_in_tertiary98)

                        DOT14_tree = self._adaptor.createWithPayload(DOT14)
                        root_0 = self._adaptor.becomeRoot(DOT14_tree, root_0)

                        self._state.following.append(self.FOLLOW_secondary_in_tertiary101)
                        secondary15 = self.secondary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, secondary15.tree)


                    elif alt2 == 3:
                        # MxAlgebra.g:12:68: AT secondary
                        pass 
                        AT16=self.match(self.input, AT, self.FOLLOW_AT_in_tertiary105)

                        AT16_tree = self._adaptor.createWithPayload(AT16)
                        root_0 = self._adaptor.becomeRoot(AT16_tree, root_0)

                        self._state.following.append(self.FOLLOW_secondary_in_tertiary108)
                        secondary17 = self.secondary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, secondary17.tree)


                    elif alt2 == 4:
                        # MxAlgebra.g:12:84: AMPERSAND secondary
                        pass 
                        AMPERSAND18=self.match(self.input, AMPERSAND, self.FOLLOW_AMPERSAND_in_tertiary112)

                        AMPERSAND18_tree = self._adaptor.createWithPayload(AMPERSAND18)
                        root_0 = self._adaptor.becomeRoot(AMPERSAND18_tree, root_0)

                        self._state.following.append(self.FOLLOW_secondary_in_tertiary115)
                        secondary19 = self.secondary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, secondary19.tree)


                    elif alt2 == 5:
                        # MxAlgebra.g:12:107: PERCENT secondary
                        pass 
                        PERCENT20=self.match(self.input, PERCENT, self.FOLLOW_PERCENT_in_tertiary119)

                        PERCENT20_tree = self._adaptor.createWithPayload(PERCENT20)
                        root_0 = self._adaptor.becomeRoot(PERCENT20_tree, root_0)

                        self._state.following.append(self.FOLLOW_secondary_in_tertiary122)
                        secondary21 = self.secondary()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, secondary21.tree)


                    else:
                        break #loop2





                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "tertiary"

    class secondary_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "secondary"
    # MxAlgebra.g:13:1: secondary : primary ( CARET secondary )? ;
    def secondary(self, ):

        retval = self.secondary_return()
        retval.start = self.input.LT(1)

        root_0 = None

        CARET23 = None
        primary22 = None

        secondary24 = None


        CARET23_tree = None

        try:
            try:
                # MxAlgebra.g:13:12: ( primary ( CARET secondary )? )
                # MxAlgebra.g:13:17: primary ( CARET secondary )?
                pass 
                root_0 = self._adaptor.nil()

                self._state.following.append(self.FOLLOW_primary_in_secondary137)
                primary22 = self.primary()

                self._state.following.pop()
                self._adaptor.addChild(root_0, primary22.tree)
                # MxAlgebra.g:13:25: ( CARET secondary )?
                alt3 = 2
                LA3_0 = self.input.LA(1)

                if (LA3_0 == CARET) :
                    alt3 = 1
                if alt3 == 1:
                    # MxAlgebra.g:13:26: CARET secondary
                    pass 
                    CARET23=self.match(self.input, CARET, self.FOLLOW_CARET_in_secondary140)

                    CARET23_tree = self._adaptor.createWithPayload(CARET23)
                    root_0 = self._adaptor.becomeRoot(CARET23_tree, root_0)

                    self._state.following.append(self.FOLLOW_secondary_in_secondary143)
                    secondary24 = self.secondary()

                    self._state.following.pop()
                    self._adaptor.addChild(root_0, secondary24.tree)






                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "secondary"

    class primary_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "primary"
    # MxAlgebra.g:15:1: primary : ( MINUS )* atom ( TILDE | SQUOTE )* ;
    def primary(self, ):

        retval = self.primary_return()
        retval.start = self.input.LT(1)

        root_0 = None

        MINUS25 = None
        TILDE27 = None
        SQUOTE28 = None
        atom26 = None


        MINUS25_tree = None
        TILDE27_tree = None
        SQUOTE28_tree = None

        try:
            try:
                # MxAlgebra.g:15:12: ( ( MINUS )* atom ( TILDE | SQUOTE )* )
                # MxAlgebra.g:15:17: ( MINUS )* atom ( TILDE | SQUOTE )*
                pass 
                root_0 = self._adaptor.nil()

                # MxAlgebra.g:15:17: ( MINUS )*
                while True: #loop4
                    alt4 = 2
                    LA4_0 = self.input.LA(1)

                    if (LA4_0 == MINUS) :
                        alt4 = 1


                    if alt4 == 1:
                        # MxAlgebra.g:15:18: MINUS
                        pass 
                        MINUS25=self.match(self.input, MINUS, self.FOLLOW_MINUS_in_primary161)

                        MINUS25_tree = self._adaptor.createWithPayload(MINUS25)
                        root_0 = self._adaptor.becomeRoot(MINUS25_tree, root_0)



                    else:
                        break #loop4


                self._state.following.append(self.FOLLOW_atom_in_primary166)
                atom26 = self.atom()

                self._state.following.pop()
                self._adaptor.addChild(root_0, atom26.tree)
                # MxAlgebra.g:15:32: ( TILDE | SQUOTE )*
                while True: #loop5
                    alt5 = 3
                    LA5_0 = self.input.LA(1)

                    if (LA5_0 == TILDE) :
                        alt5 = 1
                    elif (LA5_0 == SQUOTE) :
                        alt5 = 2


                    if alt5 == 1:
                        # MxAlgebra.g:15:34: TILDE
                        pass 
                        TILDE27=self.match(self.input, TILDE, self.FOLLOW_TILDE_in_primary170)

                        TILDE27_tree = self._adaptor.createWithPayload(TILDE27)
                        root_0 = self._adaptor.becomeRoot(TILDE27_tree, root_0)



                    elif alt5 == 2:
                        # MxAlgebra.g:15:43: SQUOTE
                        pass 
                        SQUOTE28=self.match(self.input, SQUOTE, self.FOLLOW_SQUOTE_in_primary175)

                        SQUOTE28_tree = self._adaptor.createWithPayload(SQUOTE28)
                        root_0 = self._adaptor.becomeRoot(SQUOTE28_tree, root_0)



                    else:
                        break #loop5





                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "primary"

    class atom_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "atom"
    # MxAlgebra.g:17:1: atom : ( REFERENCE | parens | function );
    def atom(self, ):

        retval = self.atom_return()
        retval.start = self.input.LT(1)

        root_0 = None

        REFERENCE29 = None
        parens30 = None

        function31 = None


        REFERENCE29_tree = None

        try:
            try:
                # MxAlgebra.g:17:12: ( REFERENCE | parens | function )
                alt6 = 3
                LA6 = self.input.LA(1)
                if LA6 == REFERENCE:
                    alt6 = 1
                elif LA6 == LPAREN:
                    alt6 = 2
                elif LA6 == FNAME:
                    alt6 = 3
                else:
                    nvae = NoViableAltException("", 6, 0, self.input)

                    raise nvae

                if alt6 == 1:
                    # MxAlgebra.g:17:17: REFERENCE
                    pass 
                    root_0 = self._adaptor.nil()

                    REFERENCE29=self.match(self.input, REFERENCE, self.FOLLOW_REFERENCE_in_atom197)

                    REFERENCE29_tree = self._adaptor.createWithPayload(REFERENCE29)
                    root_0 = self._adaptor.becomeRoot(REFERENCE29_tree, root_0)



                elif alt6 == 2:
                    # MxAlgebra.g:17:30: parens
                    pass 
                    root_0 = self._adaptor.nil()

                    self._state.following.append(self.FOLLOW_parens_in_atom202)
                    parens30 = self.parens()

                    self._state.following.pop()
                    self._adaptor.addChild(root_0, parens30.tree)


                elif alt6 == 3:
                    # MxAlgebra.g:17:39: function
                    pass 
                    root_0 = self._adaptor.nil()

                    self._state.following.append(self.FOLLOW_function_in_atom206)
                    function31 = self.function()

                    self._state.following.pop()
                    self._adaptor.addChild(root_0, function31.tree)


                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "atom"

    class parens_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "parens"
    # MxAlgebra.g:18:1: parens : LPAREN expression RPAREN ;
    def parens(self, ):

        retval = self.parens_return()
        retval.start = self.input.LT(1)

        root_0 = None

        LPAREN32 = None
        RPAREN34 = None
        expression33 = None


        LPAREN32_tree = None
        RPAREN34_tree = None

        try:
            try:
                # MxAlgebra.g:18:12: ( LPAREN expression RPAREN )
                # MxAlgebra.g:18:17: LPAREN expression RPAREN
                pass 
                root_0 = self._adaptor.nil()

                LPAREN32=self.match(self.input, LPAREN, self.FOLLOW_LPAREN_in_parens221)

                LPAREN32_tree = self._adaptor.createWithPayload(LPAREN32)
                root_0 = self._adaptor.becomeRoot(LPAREN32_tree, root_0)

                self._state.following.append(self.FOLLOW_expression_in_parens224)
                expression33 = self.expression()

                self._state.following.pop()
                self._adaptor.addChild(root_0, expression33.tree)
                RPAREN34=self.match(self.input, RPAREN, self.FOLLOW_RPAREN_in_parens226)



                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "parens"

    class function_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "function"
    # MxAlgebra.g:20:1: function : FNAME LPAREN arguments RPAREN ;
    def function(self, ):

        retval = self.function_return()
        retval.start = self.input.LT(1)

        root_0 = None

        FNAME35 = None
        LPAREN36 = None
        RPAREN38 = None
        arguments37 = None


        FNAME35_tree = None
        LPAREN36_tree = None
        RPAREN38_tree = None

        try:
            try:
                # MxAlgebra.g:20:12: ( FNAME LPAREN arguments RPAREN )
                # MxAlgebra.g:20:17: FNAME LPAREN arguments RPAREN
                pass 
                root_0 = self._adaptor.nil()

                FNAME35=self.match(self.input, FNAME, self.FOLLOW_FNAME_in_function241)

                FNAME35_tree = self._adaptor.createWithPayload(FNAME35)
                root_0 = self._adaptor.becomeRoot(FNAME35_tree, root_0)

                LPAREN36=self.match(self.input, LPAREN, self.FOLLOW_LPAREN_in_function244)
                self._state.following.append(self.FOLLOW_arguments_in_function247)
                arguments37 = self.arguments()

                self._state.following.pop()
                self._adaptor.addChild(root_0, arguments37.tree)
                RPAREN38=self.match(self.input, RPAREN, self.FOLLOW_RPAREN_in_function249)



                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "function"

    class arguments_return(ParserRuleReturnScope):
        def __init__(self):
            ParserRuleReturnScope.__init__(self)

            self.tree = None




    # $ANTLR start "arguments"
    # MxAlgebra.g:21:1: arguments : expression ( COMMA expression )* ;
    def arguments(self, ):

        retval = self.arguments_return()
        retval.start = self.input.LT(1)

        root_0 = None

        COMMA40 = None
        expression39 = None

        expression41 = None


        COMMA40_tree = None

        try:
            try:
                # MxAlgebra.g:21:12: ( expression ( COMMA expression )* )
                # MxAlgebra.g:21:17: expression ( COMMA expression )*
                pass 
                root_0 = self._adaptor.nil()

                self._state.following.append(self.FOLLOW_expression_in_arguments262)
                expression39 = self.expression()

                self._state.following.pop()
                self._adaptor.addChild(root_0, expression39.tree)
                # MxAlgebra.g:21:28: ( COMMA expression )*
                while True: #loop7
                    alt7 = 2
                    LA7_0 = self.input.LA(1)

                    if (LA7_0 == COMMA) :
                        alt7 = 1


                    if alt7 == 1:
                        # MxAlgebra.g:21:29: COMMA expression
                        pass 
                        COMMA40=self.match(self.input, COMMA, self.FOLLOW_COMMA_in_arguments265)
                        self._state.following.append(self.FOLLOW_expression_in_arguments268)
                        expression41 = self.expression()

                        self._state.following.pop()
                        self._adaptor.addChild(root_0, expression41.tree)


                    else:
                        break #loop7





                retval.stop = self.input.LT(-1)


                retval.tree = self._adaptor.rulePostProcessing(root_0)
                self._adaptor.setTokenBoundaries(retval.tree, retval.start, retval.stop)


            except RecognitionException, re:
                self.reportError(re)
                self.recover(self.input, re)
                retval.tree = self._adaptor.errorNode(self.input, retval.start, self.input.LT(-1), re)
        finally:

            pass

        return retval

    # $ANTLR end "arguments"


    # Delegated rules


 

    FOLLOW_quaternary_in_expression31 = frozenset([1])
    FOLLOW_tertiary_in_quaternary43 = frozenset([1, 4, 5, 6, 7])
    FOLLOW_PIPE_in_quaternary47 = frozenset([7, 16, 17, 19])
    FOLLOW_tertiary_in_quaternary50 = frozenset([1, 4, 5, 6, 7])
    FOLLOW_UNDERSCORE_in_quaternary54 = frozenset([7, 16, 17, 19])
    FOLLOW_tertiary_in_quaternary57 = frozenset([1, 4, 5, 6, 7])
    FOLLOW_PLUS_in_quaternary61 = frozenset([7, 16, 17, 19])
    FOLLOW_tertiary_in_quaternary64 = frozenset([1, 4, 5, 6, 7])
    FOLLOW_MINUS_in_quaternary68 = frozenset([7, 16, 17, 19])
    FOLLOW_tertiary_in_quaternary71 = frozenset([1, 4, 5, 6, 7])
    FOLLOW_secondary_in_tertiary87 = frozenset([1, 8, 9, 10, 11, 12])
    FOLLOW_ASTERISK_in_tertiary91 = frozenset([7, 16, 17, 19])
    FOLLOW_secondary_in_tertiary94 = frozenset([1, 8, 9, 10, 11, 12])
    FOLLOW_DOT_in_tertiary98 = frozenset([7, 16, 17, 19])
    FOLLOW_secondary_in_tertiary101 = frozenset([1, 8, 9, 10, 11, 12])
    FOLLOW_AT_in_tertiary105 = frozenset([7, 16, 17, 19])
    FOLLOW_secondary_in_tertiary108 = frozenset([1, 8, 9, 10, 11, 12])
    FOLLOW_AMPERSAND_in_tertiary112 = frozenset([7, 16, 17, 19])
    FOLLOW_secondary_in_tertiary115 = frozenset([1, 8, 9, 10, 11, 12])
    FOLLOW_PERCENT_in_tertiary119 = frozenset([7, 16, 17, 19])
    FOLLOW_secondary_in_tertiary122 = frozenset([1, 8, 9, 10, 11, 12])
    FOLLOW_primary_in_secondary137 = frozenset([1, 13])
    FOLLOW_CARET_in_secondary140 = frozenset([7, 16, 17, 19])
    FOLLOW_secondary_in_secondary143 = frozenset([1])
    FOLLOW_MINUS_in_primary161 = frozenset([7, 16, 17, 19])
    FOLLOW_atom_in_primary166 = frozenset([1, 14, 15])
    FOLLOW_TILDE_in_primary170 = frozenset([1, 14, 15])
    FOLLOW_SQUOTE_in_primary175 = frozenset([1, 14, 15])
    FOLLOW_REFERENCE_in_atom197 = frozenset([1])
    FOLLOW_parens_in_atom202 = frozenset([1])
    FOLLOW_function_in_atom206 = frozenset([1])
    FOLLOW_LPAREN_in_parens221 = frozenset([7, 16, 17, 19])
    FOLLOW_expression_in_parens224 = frozenset([18])
    FOLLOW_RPAREN_in_parens226 = frozenset([1])
    FOLLOW_FNAME_in_function241 = frozenset([17])
    FOLLOW_LPAREN_in_function244 = frozenset([7, 16, 17, 19])
    FOLLOW_arguments_in_function247 = frozenset([18])
    FOLLOW_RPAREN_in_function249 = frozenset([1])
    FOLLOW_expression_in_arguments262 = frozenset([1, 20])
    FOLLOW_COMMA_in_arguments265 = frozenset([7, 16, 17, 19])
    FOLLOW_expression_in_arguments268 = frozenset([1, 20])



def main(argv, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr):
    from antlr3.main import ParserMain
    main = ParserMain("MxAlgebraLexer", MxAlgebraParser)
    main.stdin = stdin
    main.stdout = stdout
    main.stderr = stderr
    main.execute(argv)


if __name__ == '__main__':
    main(sys.argv)
