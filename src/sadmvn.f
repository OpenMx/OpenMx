      SUBROUTINE RANMVN( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM, TID )
*
*     A subroutine for computing multivariate normal probabilities.
*     This subroutine uses the Monte-Carlo algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*     MAXPTS INTEGER, maximum number of function values allowed. This
*            parameter can be used to limit the time taken. A
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS
*                           function vaules used; increase MAXPTS to
*                           decrease ERROR;
*            if INFORM = 2, N > 100 or N < 1.
*
      EXTERNAL MVNFNC
      INTEGER N, INFIN(*), MAXPTS, MPT, 
     &     INFORM, INFIS, IVLS, TID
      DOUBLE PRECISION
     &     CORREL(*), LOWER(*), UPPER(*), MVNFNC,
     &     ABSEPS, RELEPS, ERROR, VALUE, D, E, EPS, MVNNIT
      IF ( N .GT. 100 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
         RETURN
      ENDIF
      INFORM = MVNNIT(N, CORREL, LOWER, UPPER, 
     &                INFIN, INFIS, D, E, TID)
      IF ( N-INFIS .EQ. 0 ) THEN
         VALUE = 1
         ERROR = 0
      ELSE IF ( N-INFIS .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call then Monte-Carlo integration subroutine
*
         MPT = 25 + 10*N
         CALL RCRUDE(N-INFIS-1,MPT,MVNFNC,ERROR,VALUE,0,TID)
         IVLS = MPT
 10      EPS = MAX( ABSEPS, RELEPS*ABS(VALUE) )
         IF ( ERROR .GT. EPS .AND. IVLS .LT. MAXPTS ) THEN
            MPT = MAX( MIN( INT(MPT*(ERROR/(EPS))**2),
     &                      MAXPTS-IVLS ), 10 )
            CALL RCRUDE(N-INFIS-1,MPT,MVNFNC,ERROR,VALUE,1,TID)
            IVLS = IVLS + MPT
            GO TO 10
         ENDIF
         IF ( ERROR. GT. EPS .AND. IVLS .GE. MAXPTS ) INFORM = 1
      ENDIF
      END
      SUBROUTINE RCRUDE(NDIM,MAXPTS,FUNCTN,ABSEST,FINEST,IR,TID)
*
*     Crude Monte-Carlo Algorithm with simple antithetic variates
*      and weighted results on restart
*
      EXTERNAL FUNCTN
      INTEGER NDIM, MAXPTS, M, K, IR, NPTS, TID
      DOUBLE PRECISION FINEST, ABSEST, X(100), FUN, FUNCTN, UNI,
     &     VARSQR, VAREST, VARPRD, FINDIF, FINVAL
      SAVE VAREST
      IF ( IR .LE. 0 ) THEN
         VAREST = 0
         FINEST = 0
      ENDIF
      FINVAL = 0
      VARSQR = 0
      NPTS = MAXPTS/2
      DO M = 1,NPTS
         DO K = 1,NDIM
            X(K) = UNI()
         END DO
         FUN = FUNCTN(NDIM, X, TID)
         DO K = 1,NDIM
            X(K) = 1 - X(K)
         END DO
         FUN = ( FUNCTN(NDIM, X, TID) + FUN )/2
         FINDIF = ( FUN - FINVAL )/M
         VARSQR = ( M - 2 )*VARSQR/M + FINDIF**2
         FINVAL = FINVAL + FINDIF
      END DO
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/(1 + VARPRD)
      IF ( VARSQR .GT. 0 ) VAREST = (1 + VARPRD)/VARSQR
      ABSEST = 3*SQRT( VARSQR/( 1 + VARPRD ) )
      END
      SUBROUTINE SPHMVN(N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     &     ABSEPS, RELEPS, ERROR, VALUE, INFORM)
*
*     A subroutine for computing multivariate normal probabilities.
*     This subroutine uses a Mont-Carlo algorithm given in the paper
*       "Three Digit Accurate Multiple Normal Probabilities",
*          pp. 369-380, Numer. Math. 35(1980), by I. Deak
*
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*     MAXPTS INTEGER, maximum number of function values allowed. This
*            parameter can be used to limit the time. A sensible
*            strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL, estimated absolute error, with 99% confidence level.
*     VALUE  REAL, estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS
*                           function vaules used; increase MAXPTS to
*                           decrease ERROR;
*            if INFORM = 2, N > 100.
*
      EXTERNAL SPNRML
      INTEGER N, INFIS, INFIN(*), MAXPTS, MPT, INFORM, NS, IVLS
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*),
     &     ABSEPS, RELEPS, ERROR, VALUE, D, E, EPS, SPNRNT
      IF ( N .GT. 100 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
         RETURN
      ENDIF
      INFORM = SPNRNT(N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E, NS)
      IF ( N-INFIS .EQ. 0 ) THEN
         VALUE = 1
         ERROR = 0
      ELSE IF ( N-INFIS .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call then Monte-Carlo integration subroutine
*
         MPT = 25 + NS/N**3
         CALL SCRUDE( N-INFIS, MPT, ERROR, VALUE, 0 )
         IVLS = MPT*NS
 10      EPS = MAX( ABSEPS, RELEPS*ABS(VALUE) )
         IF ( ERROR .GT. EPS .AND. IVLS .LT. MAXPTS ) THEN
            MPT = MAX( MIN( INT(MPT*(ERROR/(EPS))**2),
     &                      ( MAXPTS - IVLS )/NS ), 10 )
            CALL SCRUDE( N-INFIS, MPT, ERROR, VALUE, 1 )
            IVLS = IVLS + MPT*NS
            GO TO 10
         ENDIF
         IF ( ERROR. GT. EPS .AND. IVLS .GE. MAXPTS ) INFORM = 1
      ENDIF
      END
      DOUBLE PRECISION FUNCTION SPNRML(N)
*
*     Integrand subroutine
*
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL(*), D, E, ZERO
      INTEGER N, INFIN(*), INFIS
      INTEGER NL, IJ, I, J, K, NS, NSO, ND
      PARAMETER ( NL = 100, ND = 3, ZERO = 0 )
      DOUBLE PRECISION A(NL), B(NL), U(NL,NL), Y(NL)
      INTEGER INFI(NL), IS(NL), IC(NL)
      DOUBLE PRECISION RS, TMP, BT, RNOR, SPHLIM, SPNRNT
      SAVE A, B, INFI, U
*
*    First generate U = COV*(random orthogonal matrix)
*
      DO K = N-1, 1, -1
         TMP = 0
         DO J = K, N
            Y(J) = RNOR()
            TMP = TMP + Y(J)**2
         END DO
         TMP = -SQRT(TMP)
         BT = 1/( TMP*( Y(K) + TMP ) )
         Y(K) = Y(K) + TMP
         DO I = 1, N
            TMP = 0
            DO J = K, N
               TMP = TMP + U(I,J)*Y(J)
            END DO
            TMP = BT*TMP
            DO J = K, N
               U(I,J) = U(I,J) - TMP*Y(J)
            END DO
         END DO
      END DO
*
*     Compute integrand average
*
      RS = SQRT( DBLE(ND) )
      DO I = 1,ND
         IC(I) = I
      END DO
      IC(ND+1) = N+1
      SPNRML = 0
      NS = 0
 10   DO I = 1,ND
         IS(I) = -1
      END DO
 20   DO I = 1, N
         TMP = 0
         DO J = 1,ND
            TMP = TMP + IS(J)*U(I,IC(J))
         END DO
         Y(I) = TMP/RS
      END DO
      NS = NS + 1
      SPNRML = SPNRML + ( SPHLIM( N, A, B, INFI, Y ) - SPNRML )/NS
      DO I = 1, ND
         IS(I) = IS(I) + 2
         IF ( IS(I) .LT. 2 ) GO TO 20
         IS(I) = -1
      END DO
      DO I = 1, ND
         IC(I) = IC(I) + 1
         IF ( IC(I) .LT. IC(I+1)  ) GO TO 10
         IC(I) = I
      END DO
      SPNRML = SPNRML/2
      RETURN
      ENTRY SPNRNT( N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E, NSO )
      SPNRNT = 0
*
*     Initialisation
*
      IJ = 0
      INFIS = 0
      DO I = 1, N
         INFI(I) = INFIN(I)
         IF ( INFI(I) .LT. 0 ) THEN
            INFIS = INFIS + 1
         ELSE
            A(I) = 0
            B(I) = 0
            IF ( INFI(I) .NE. 0 ) A(I) = LOWER(I)
            IF ( INFI(I) .NE. 1 ) B(I) = UPPER(I)
         ENDIF
         DO J = 1, I-1
            IJ = IJ + 1
            U(I,J) = CORREL(IJ)
            U(J,I) = 0
         END DO
         U(I,I) = 1
      END DO
      NSO = 1
      DO I = 1,ND
         NSO = 2*NSO*( N - INFIS - I + 1 )/I
      END DO
*
*     First move any doubly infinite limits to innermost positions
*
      IF ( INFIS .LT. N ) THEN
         outer: DO I = N, N-INFIS+1, -1
            IF ( INFI(I) .GE. 0 ) THEN
               DO J = 1,I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     DO K = 1, J-1
                        TMP = U(J,K)
                        U(J,K) = U(I,K)
                        U(I,K) = TMP
                     END DO
                     DO K = J+1, I-1
                        TMP = U(I,K)
                        U(I,K) = U(K,J)
                        U(K,J) = TMP
                     END DO
                     DO K = I+1, N
                        TMP = U(K,J)
                        U(K,J) = U(K,I)
                        U(K,I) = TMP
                     END DO
                     TMP = A(J)
                     A(J) = A(I)
                     A(I) = TMP
                     TMP = B(J)
                     B(J) = B(I)
                     B(I) = TMP
                     TMP = INFI(J)
                     INFI(J) = INFI(I)
                     INFI(I) = TMP
                     CYCLE outer
                  ENDIF
               END DO
            ENDIF
         END DO outer
      ENDIF
*
*     Determine Cholesky decomposition
*
      DO J = 1, N-INFIS
         DO I = J, N-INFIS
            TMP = U(I,J)
            DO K = 1, J-1
               TMP = TMP - U(I,K)*U(J,K)
            END DO
            IF ( I .EQ. J ) THEN
               U(J,J) = SQRT( MAX( TMP, ZERO ) )
            ELSE IF ( U(I,I) .GT. 0 ) THEN
               U(I,J) = TMP/U(J,J)
            ELSE
               U(I,J) = 0
            END IF
         END DO
      END DO
      DO I = 1, N-INFIS
         IF ( U(I,I) .GT. 0 ) THEN
            IF ( INFI(I) .NE. 0 ) A(I) = A(I)/U(I,I)
            IF ( INFI(I) .NE. 1 ) B(I) = B(I)/U(I,I)
            DO J = 1,I
               U(I,J) = U(I,J)/U(I,I)
            END DO
         ENDIF
      END DO
      CALL LIMITS( A(1), B(1), INFI(1), D, E )
      END
      DOUBLE PRECISION FUNCTION SPHLIM( N, A, B, INFI, Y )
      DOUBLE PRECISION A(*), B(*), Y(*), CMN, CMX, SPHINC
      INTEGER INFI(*), I, N
      CMN = -10*N
      CMX =  10*N
      DO I = 1,N
         IF ( Y(I) .GT. 0 ) THEN
            IF ( INFI(I) .NE. 1 ) CMX = MIN( CMX, B(I)/Y(I) )
            IF ( INFI(I) .NE. 0 ) CMN = MAX( CMN, A(I)/Y(I) )
         ELSE
            IF ( INFI(I) .NE. 1 ) CMN = MAX( CMN, B(I)/Y(I) )
            IF ( INFI(I) .NE. 0 ) CMX = MIN( CMX, A(I)/Y(I) )
         ENDIF
      END DO
      IF ( CMN .LT. CMX ) THEN
         IF ( CMN .GE. 0 .AND. CMX .GE. 0 ) THEN
            SPHLIM = SPHINC( N,  CMX ) - SPHINC( N,  CMN )
         ELSEIF ( CMN .LT. 0 .AND. CMX .GE. 0 ) THEN
            SPHLIM = SPHINC( N, -CMN ) + SPHINC( N,  CMX )
         ELSE
            SPHLIM = SPHINC( N, -CMN ) - SPHINC( N, -CMX )
         ENDIF
      ELSE
         SPHLIM = 0
      ENDIF
      END
      SUBROUTINE SCRUDE( NDIM, MAXPTS, ABSEST, FINEST, IR )
*
*     Crude Monte-Carlo Algorithm for Deak method with
*      weighted results on restart
*
      INTEGER NDIM, MAXPTS, M, K, IR, NPTS
      DOUBLE PRECISION FINEST, ABSEST, SPNRML, UNI,
     &     VARSQR, VAREST, VARPRD, FINDIF, FINVAL
      SAVE VAREST
      IF ( IR .LE. 0 ) THEN
         VAREST = 0
         FINEST = 0
      ENDIF
      FINVAL = 0
      VARSQR = 0
      DO M = 1,MAXPTS
         FINDIF = ( SPNRML(NDIM) - FINVAL )/M
         FINVAL = FINVAL + FINDIF
         VARSQR = ( M - 2 )*VARSQR/M + FINDIF**2
      END DO
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/(1 + VARPRD)
      IF ( VARSQR .GT. 0 ) VAREST = (1 + VARPRD)/VARSQR
      ABSEST = 3*SQRT( VARSQR/( 1 + VARPRD ) )
      END
      DOUBLE PRECISION FUNCTION SPHINC( N, R )
*
*                   R
*     SPHINC =  K  I  exp(-t*t/2) t**(N-1) dt, for N > 1.
*                N  0
*
      INTEGER I, N
      DOUBLE PRECISION R, RR, RP, PF, ET, PHI
      PARAMETER ( RP = 2.5066 28274 63100 04D0 )
      IF ( R .GT. 0 ) THEN
         RR = R*R
         PF = 1
         DO I = N-2, 2, -2
            PF = 1 + RR*PF/I
         END DO
         IF ( MOD( N, 2 ) .EQ. 0 ) THEN
            ET = LOG(PF) - RR/2
            IF ( ET .GT. -40 ) THEN
               SPHINC = 1 - EXP( ET )
            ELSE
               SPHINC = 1
            END IF
         ELSE
            SPHINC = 1  - 2*PHI(-R)
            ET = LOG(R*PF) - RR/2
            IF ( ET .GT. -40 ) SPHINC = SPHINC - 2*EXP( ET )/RP
         ENDIF
      ELSE
         SPHINC = 0
      ENDIF
      END
      DOUBLE PRECISION FUNCTION RNOR()
*
*     RNOR generates normal random numbers with zero mean and unit
*     standard deviation, often denoted N(0,1),adapted from G. Marsaglia
*     and W. W. Tsang: "A Fast, Easily Implemented Method for Sampling
*     from Decreasing or Symmetric Unimodal Density Functions"
*      SIAM J. Sci. Stat. Comput. 5(1984), pp. 349-359.
*
      INTEGER J, N, TN
      DOUBLE PRECISION TWOPIS, AA, B, C, XDN
      PARAMETER ( N = 64, TN = 2*N, TWOPIS = TN/2.506628274631000D0 )
      PARAMETER ( XDN = 0.3601015713011893D0, B = 0.4878991777603940D0 )
      PARAMETER (  AA =  12.37586029917064D0, C =  12.67705807886560D0 )
      DOUBLE PRECISION XT, XX, Y, UNI
      DOUBLE PRECISION X(0:N)
      SAVE X
      DATA ( X(J), J = 0, 31 ) /
     &  0.3409450287039653D+00,  0.4573145918669259D+00,
     &  0.5397792816116612D+00,  0.6062426796530441D+00,
     &  0.6631690627645207D+00,  0.7136974590560222D+00,
     &  0.7596124749339174D+00,  0.8020356003555283D+00,
     &  0.8417226679789527D+00,  0.8792102232083114D+00,
     &  0.9148948043867484D+00,  0.9490791137530882D+00,
     &  0.9820004812398864D+00,  0.1013849238029940D+01,
     &  0.1044781036740172D+01,  0.1074925382028552D+01,
     &  0.1104391702268125D+01,  0.1133273776243940D+01,
     &  0.1161653030133931D+01,  0.1189601040838737D+01,
     &  0.1217181470700870D+01,  0.1244451587898246D+01,
     &  0.1271463480572119D+01,  0.1298265041883197D+01,
     &  0.1324900782180860D+01,  0.1351412509933371D+01,
     &  0.1377839912870011D+01,  0.1404221063559975D+01,
     &  0.1430592868502691D+01,  0.1456991476137671D+01,
     &  0.1483452656603219D+01,  0.1510012164318519D+01 /
      DATA ( X(J), J = 32, 64 ) /
     &  0.1536706093359520D+01,  0.1563571235037691D+01,
     &  0.1590645447014253D+01,  0.1617968043674446D+01,
     &  0.1645580218369081D+01,  0.1673525509567038D+01,
     &  0.1701850325062740D+01,  0.1730604541317782D+01,
     &  0.1759842199038300D+01,  0.1789622321566574D+01,
     &  0.1820009890130691D+01,  0.1851077020230275D+01,
     &  0.1882904397592872D+01,  0.1915583051943031D+01,
     &  0.1949216574916360D+01,  0.1983923928905685D+01,
     &  0.2019843052906235D+01,  0.2057135559990095D+01,
     &  0.2095992956249391D+01,  0.2136645022544389D+01,
     &  0.2179371340398135D+01,  0.2224517507216017D+01,
     &  0.2272518554850147D+01,  0.2323933820094302D+01,
     &  0.2379500774082828D+01,  0.2440221797979943D+01,
     &  0.2507511701865317D+01,  0.2583465835225429D+01,
     &  0.2671391590320836D+01,4*0.2776994269662875D+01 /
      Y = UNI()
      J = MOD( INT( TN*UNI() ), N )
      XT = X(J+1)
      RNOR = ( Y + Y - 1 )*XT
      IF ( ABS(RNOR) .GT. X(J) ) THEN
         XX = B*( XT - ABS(RNOR) )/( XT - X(J) )
         Y = UNI()
         IF ( Y .GT. C - AA*EXP( -XX**2/2 ) ) THEN
            RNOR = SIGN( XX, RNOR )
         ELSE
            IF ( EXP(-XT**2/2)+Y/(TWOPIS*XT).GT.EXP(-RNOR**2/2) ) THEN
 10            XX = XDN*LOG( UNI() )
               IF ( -2*LOG( UNI() ) .LE. XX**2 ) GO TO 10
               RNOR = SIGN( X(N) - XX, RNOR )
            END IF
         END IF
      END IF
      END
*
      DOUBLE PRECISION FUNCTION UNI()
*
*     Uniform (0, 1) random number generator
*
*     Reference:
*     L'Ecuyer, Pierre (1996),
*     "Combined Multiple Recursive Random Number Generators"
*     Operations Research 44, pp. 816-822.
*
*
      INTEGER A12, A13, A21, A23, P12, P13, P21, P23
      INTEGER Q12, Q13, Q21, Q23, R12, R13, R21, R23
      INTEGER X10, X11, X12, X20, X21, X22, Z, M1, M2, H
      DOUBLE PRECISION INVMP1
      PARAMETER (  M1 = 2147483647,  M2 = 2145483479 )
      PARAMETER ( A12 =   63308,    Q12 = 33921, R12 = 12979 )
      PARAMETER ( A13 = -183326,    Q13 = 11714, R13 =  2883 )
      PARAMETER ( A21 =   86098,    Q21 = 24919, R21 =  7417 )
      PARAMETER ( A23 = -539608,    Q23 =  3976, R23 =  2071 )
      PARAMETER ( INVMP1 = 4.656612873077392578125D-10 )
*                 INVMP1 = 1/(M1+1)
      SAVE X10, X11, X12, X20, X21, X22
      DATA       X10,      X11,      X12,      X20,      X21,      X22
     &    / 11111111, 22222223, 33333335, 44444447, 55555559, 66666661 /
*
*     Component 1
*
      H = X10/Q13
      P13 = -A13*( X10 - H*Q13 ) - H*R13
      H = X11/Q12
      P12 =  A12*( X11 - H*Q12 ) - H*R12
      IF ( P13 .LT. 0 ) P13 = P13 + M1
      IF ( P12 .LT. 0 ) P12 = P12 + M1
      X10 = X11
      X11 = X12
      X12 = P12 - P13
      IF ( X12 .LT. 0 ) X12 = X12 + M1
*
*     Component 2
*
      H = X20/Q23
      P23 = -A23*( X20 - H*Q23 ) - H*R23
      H = X22/Q21
      P21 =  A21*( X22 - H*Q21 ) - H*R21
      IF ( P23 .LT. 0 ) P23 = P23 + M2
      IF ( P21 .LT. 0 ) P21 = P21 + M2
      X20 = X21
      X21 = X22
      X22 = P21 - P23
      IF ( X22 .LT. 0 ) X22 = X22 + M2
*
*     Combination
*
      Z = X12 - X22
      IF ( Z .LE. 0 ) Z = Z + M1
      UNI = Z*INVMP1
      END
      SUBROUTINE SADMVN( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM, TID )
*
*     A subroutine for computing multivariate normal probabilities.
*     This subroutine uses an algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*     MAXPTS INTEGER, maximum number of function values allowed. This
*            parameter can be used to limit the time taken. A
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS
*                           function vaules used; increase MAXPTS to
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*
      EXTERNAL MVNFNC
      INTEGER N, NL, M, INFIN(*), LENWRK, MAXPTS, INFORM, INFIS,
     &     RULCLS, TOTCLS, NEWCLS, MAXCLS, TID, I, NNOTINF
      DOUBLE PRECISION
     &     CORREL(*), LOWER(*), UPPER(*), ABSEPS, RELEPS, ERROR, VALUE,
     &     OLDVAL, D, E, MVNNIT, MVNFNC
      PARAMETER ( NL = 20 )
      PARAMETER ( LENWRK = 20*NL**2 )
      PARAMETER ( NTHREADS = 64 )
      DOUBLE PRECISION WORK(LENWRK, NTHREADS)
*
*	MCN change test to number of dimensions with INFIN(I) not <0
*	 
	NNOTINF = 0
	DO I = 1, N 
		IF (INFIN(I) .GE. 0) THEN
			NNOTINF = NNOTINF + 1
		ENDIF
	END DO
*	WRITE(6,*) 'NNotinf = ',NNOTINF
      IF ( NNOTINF .GT. 20 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
         RETURN
      ENDIF
      INFORM = MVNNIT( N, CORREL, LOWER, UPPER, 
     &                 INFIN, INFIS, D, E, TID )
      M = N - INFIS
      IF ( M .EQ. 0 ) THEN
         VALUE = 1
         ERROR = 0
      ELSE IF ( M .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call the subregion adaptive integration subroutine
*
         M = M - 1
         RULCLS = 1
         CALL ADAPT( M, RULCLS, 0, MVNFNC, ABSEPS, RELEPS,
     &               LENWRK, NTHREADS, WORK, ERROR, VALUE, INFORM, TID )
         MAXCLS = MIN( 10*RULCLS, MAXPTS )
         TOTCLS = 0
         CALL ADAPT(M, TOTCLS, MAXCLS, MVNFNC, ABSEPS, RELEPS,
     &        LENWRK, NTHREADS, WORK, ERROR, VALUE, INFORM, TID)
         IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
 10         OLDVAL = VALUE
            MAXCLS = MAX( 2*RULCLS, MIN( 3*MAXCLS/2, MAXPTS - TOTCLS ) )
            NEWCLS = -1
            CALL ADAPT(M, NEWCLS, MAXCLS, MVNFNC, ABSEPS, RELEPS,
     &           LENWRK, NTHREADS, WORK, ERROR, VALUE, INFORM, TID)
            TOTCLS = TOTCLS + NEWCLS
            ERROR = ABS(VALUE-OLDVAL) + SQRT(RULCLS*ERROR**2/TOTCLS)
            IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
               IF ( MAXPTS - TOTCLS .GT. 2*RULCLS ) GO TO 10
            ELSE
               INFORM = 0
            END IF
         ENDIF
      ENDIF
      END
*
      SUBROUTINE KROMVN( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM, TID )
*
*     A subroutine for computing multivariate normal probabilities.
*     This subroutine uses an algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : AlanGenz@wsu.edu
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*     MAXPTS INTEGER, maximum number of function values allowed. This
*            parameter can be used to limit the time. A sensible
*            strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS
*                           function vaules used; increase MAXPTS to
*                           decrease ERROR;
*            if INFORM = 2, N > 100 or N < 1.
*
      EXTERNAL MVNFNC
      INTEGER N, INFIN(*), MAXPTS, INFORM, INFIS, IVLS, TID
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*), RELEPS, ABSEPS,
     &       ERROR, VALUE, E, D, MVNNIT, MVNFNC
      IF ( N .GT. 100 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
      ELSE
         INFORM = MVNNIT(N, CORREL, LOWER, UPPER, 
     &                   INFIN, INFIS, D, E, TID)
         IF ( N-INFIS .EQ. 0 ) THEN
            VALUE = 1
            ERROR = 0
         ELSE IF ( N-INFIS .EQ. 1 ) THEN
            VALUE = E - D
            ERROR = 2E-16
         ELSE
*
*        Call the lattice rule integration subroutine
*
            IVLS = 0
            CALL KROBOV( N-INFIS-1, IVLS, MAXPTS, MVNFNC,
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM, TID )
         ENDIF
      ENDIF
      END
      SUBROUTINE KROBOV( NDIM, MINVLS, MAXVLS, FUNCTN, ABSEPS, RELEPS,
     &                   ABSERR, FINEST, INFORM, TID )
*
*  Automatic Multidimensional Integration Subroutine
*
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 4/15/98
*
*  KROBOV computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  KROBOV uses randomized Korobov rules. The primary references are
*  "Randomization of Number Theoretic Methods for Multiple Integration"
*   R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and
*   "Optimal Parameters for Multidimensional Integration",
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*
***************  Parameters ********************************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 40
*  MINVLS  Integer minimum number of function evaluations allowed.
*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
*          routine assumes a previous call has been made with
*          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
*          It must have parameters (NDIM,Z), where Z is a real array
*          of dimension NDIM.
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
****** Output parameters
*  MINVLS  Actual number of function evaluations used.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and
*                     INTVLS <= MAXCLS.
*          INFORM = 1 If MAXVLS was too small to obtain the required
*          accuracy. In this case a value FINEST is returned with
*          estimated absolute accuracy ABSERR.
************************************************************************
      EXTERNAL FUNCTN
      INTEGER NDIM, MINVLS, MAXVLS, INFORM, NP, PLIM, NLIM,
     &        SAMPLS, I, INTVLS, MINSMP, TID
      PARAMETER ( PLIM = 20, NLIM = 100, MINSMP = 6 )
      INTEGER C(PLIM,NLIM), P(PLIM)
      DOUBLE PRECISION FUNCTN, ABSEPS, RELEPS, FINEST, ABSERR, DIFINT,
     &       FINVAL, VARSQR, VAREST, VARPRD, VALUE
      DOUBLE PRECISION ALPHA(NLIM), X(NLIM), VK(NLIM), ONE
      PARAMETER ( ONE = 1 )
      SAVE P, C, SAMPLS, NP, VAREST
      INFORM = 1
      INTVLS = 0
      IF ( MINVLS .GE. 0 ) THEN
         FINEST = 0
         VAREST = 0
         SAMPLS = MINSMP
         DO I = 1, PLIM
            NP = I
            IF ( MINVLS .LT. 2*SAMPLS*P(I) ) GO TO 10
         END DO
         SAMPLS = MAX( MINSMP, MINVLS/( 2*P(NP) ) )
      ENDIF
 10   VK(1) = ONE/P(NP)
      DO I = 2, NDIM
         VK(I) = MOD( C(NP,NDIM-1)*VK(I-1), ONE )
      END DO
      FINVAL = 0
      VARSQR = 0
      DO I = 1, SAMPLS
         CALL KROSUM( NDIM, VALUE, P(NP), VK, FUNCTN, ALPHA, X, TID )
         DIFINT = ( VALUE - FINVAL )/I
         FINVAL = FINVAL + DIFINT
         VARSQR = ( I - 2 )*VARSQR/I + DIFINT**2
      END DO
      INTVLS = INTVLS + 2*SAMPLS*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1 + VARPRD )
      IF ( VARSQR .GT. 0 ) VAREST = ( 1 + VARPRD )/VARSQR
      ABSERR = 3*SQRT( VARSQR/( 1 + VARPRD ) )
      IF ( ABSERR .GT. MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN
         IF ( NP .LT. PLIM ) THEN
            NP = NP + 1
         ELSE
            SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P(NP) ) )
            SAMPLS = MAX( MINSMP, SAMPLS )
         ENDIF
         IF ( INTVLS + 2*SAMPLS*P(NP) .LE. MAXVLS ) GO TO 10
      ELSE
         INFORM = 0
      ENDIF
      MINVLS = INTVLS
      DATA P( 1), ( C( 1,I), I = 1, 99 ) /    113,
     &     42,    54,    55,    32,    13,    26,    26,    13,    26,
     &     14,    13,    26,    35,     2,     2,     2,     2,    56,
     &     28,     7,     7,    28,     4,    49,     4,    40,    48,
     &      5,    35,    27,    16,    16,     2,     2,     7,    28,
     &      4,    49,     4,    56,     8,     2,     2,    56,     7,
     &     16,    28,     7,     7,    28,     4,    49,     4,    37,
     &     55,    21,    33,    40,    16,    16,    28,     7,    16,
     &     28,     4,    49,     4,    56,    35,     2,     2,     2,
     &     16,    16,    28,     4,    16,    28,     4,    49,     4,
     &     40,    40,     5,    42,    27,    16,    16,    28,     4,
     &     16,    28,     4,    49,     4,     8,     8,     2,     2/
      DATA P( 2), ( C( 2,I), I = 1, 99 ) /    173,
     &     64,    34,    57,     9,    72,    86,    16,    75,    75,
     &     70,    42,     2,    86,    62,    62,    30,    30,     5,
     &     42,    70,    70,    70,    53,    70,    70,    53,    42,
     &     62,    53,    53,    53,    69,    75,     5,    53,    86,
     &      2,     5,    30,    75,    59,     2,    69,     5,     5,
     &     63,    62,     5,    69,    30,    44,    30,    86,    86,
     &      2,    69,     5,     5,     2,     2,    61,    69,    17,
     &      2,     2,     2,    53,    69,     2,     2,    86,    69,
     &     13,     2,     2,    37,    43,    65,     2,     2,    30,
     &     86,    45,    16,    32,    18,    86,    86,    86,     9,
     &     63,    63,    11,    76,    76,    76,    63,    60,    70/
      DATA P( 3), ( C( 3,I), I = 1, 99 ) /    263,
     &    111,    67,    98,    36,    48,   110,     2,   131,     2,
     &      2,   124,   124,    48,     2,     2,   124,   124,    70,
     &     70,    48,   126,    48,   126,    56,    65,    48,    48,
     &     70,     2,    92,   124,    92,   126,   131,   124,    70,
     &     70,    70,    20,   105,    70,     2,     2,    27,   108,
     &     27,    39,     2,   131,   131,    92,    92,    48,     2,
     &    126,    20,   126,     2,     2,   131,    38,   117,     2,
     &    131,    68,    58,    38,    90,    38,   108,    38,     2,
     &    131,   131,   131,    68,    14,    94,   131,   131,   131,
     &    108,    18,   131,    56,    85,   117,   117,     9,   131,
     &    131,    55,    92,    92,    92,   131,   131,    48,    48/
      DATA P( 4), ( C( 4,I), I = 1, 99 ) /    397,
     &    151,   168,    46,   197,    69,    64,     2,   198,   191,
     &    134,   134,   167,   124,    16,   124,   124,   124,   124,
     &    141,   134,   128,     2,     2,    32,    32,    32,    31,
     &     31,    64,    64,    99,     4,     4,   167,   124,   124,
     &    124,   124,   124,   124,   107,    85,    79,    85,   111,
     &     85,   128,    31,    31,    31,    31,    64,   167,     4,
     &    107,   167,   124,   124,   124,   124,   124,   124,   107,
     &    183,     2,     2,     2,    62,    32,    31,    31,    31,
     &     31,    31,   167,     4,   107,   167,   124,   124,   124,
     &    124,   124,   124,   107,   142,   184,   184,    65,    65,
     &    183,    31,    31,    31,    31,    31,   167,     4,   107/
      DATA P( 5), ( C( 5,I), I = 1, 99 ) /    593,
     &    229,    40,   268,    42,   153,   294,    71,     2,   130,
     &    199,   199,   199,   149,   199,   149,   153,   130,   149,
     &    149,    15,   119,   294,    31,    82,   260,   122,   209,
     &    209,   122,   296,   130,   130,   260,   260,    30,   206,
     &     94,   209,    94,   122,   209,   209,   122,   122,   209,
     &    130,     2,   130,   130,    38,    38,    79,    82,    94,
     &     82,   122,   122,   209,   209,   122,   122,   168,   220,
     &     62,    60,   168,   282,   282,    82,   209,   122,    94,
     &    209,   122,   122,   122,   122,   258,   148,   286,   256,
     &    256,    62,    62,    82,   122,    82,    82,   122,   122,
     &    122,   209,   122,    15,    79,    79,    79,    79,   168/
      DATA P( 6), ( C( 6,I), I = 1, 99 ) /    907,
     &    264,   402,   406,   147,   452,   153,   224,     2,     2,
     &    224,   224,   449,   101,   182,   449,   101,   451,   181,
     &    181,   101,   101,   377,    85,   453,   453,   453,    85,
     &    197,   451,     2,     2,   101,   449,   449,   449,   173,
     &    173,     2,   453,   453,     2,   426,    66,   367,   426,
     &    101,   453,     2,    32,    32,    32,   101,     2,     2,
     &    453,   223,   147,   449,   290,     2,   453,     2,    83,
     &    223,   101,   453,     2,    83,    83,   147,     2,   453,
     &    147,   147,   147,   147,   147,   147,   147,   453,   153,
     &    153,   147,     2,   224,   290,   320,   453,   147,   431,
     &    383,   290,   290,     2,   162,   162,   147,     2,   162/
      DATA P( 7), ( C( 7,I), I = 1, 99 ) /   1361,
     &    505,   220,   195,   410,   199,   248,   460,   471,     2,
     &    331,   662,   547,   209,   547,   547,   209,     2,   680,
     &    680,   629,   370,   574,    63,    63,   259,   268,   259,
     &    547,   209,   209,   209,   547,   547,   209,   209,   547,
     &    547,   108,    63,    63,   108,    63,    63,   108,   259,
     &    268,   268,   547,   209,   209,   209,   209,   547,   209,
     &    209,   209,   547,   108,    63,    63,    63,   405,   285,
     &    234,   259,   259,   259,   259,   209,   209,   209,   209,
     &    209,   209,   209,   209,   547,   289,   289,   234,   285,
     &    316,     2,   410,   259,   259,   259,   268,   209,   209,
     &    209,   209,   547,   547,   209,   209,   209,   285,   316/
      DATA P( 8), ( C( 8,I), I = 1, 99 ) /   2053,
     &    468,   635,   849,   687,   948,    37,  1014,   513,     2,
     &      2,     2,     2,     2,  1026,     2,     2,  1026,   201,
     &    201,     2,  1026,   413,  1026,  1026,     2,     2,   703,
     &    703,     2,     2,   393,   393,   678,   413,  1026,     2,
     &      2,  1026,  1026,     2,   405,   953,     2,  1026,   123,
     &    123,   953,   953,   123,   405,   794,   123,   647,   613,
     &   1026,   647,   768,   953,   405,   953,   405,   918,   918,
     &    123,   953,   953,   918,   953,   536,   405,    70,   124,
     &   1005,   529,   207,   405,   405,   953,   953,   123,   918,
     &    918,   953,   405,   918,   953,   468,   405,   794,   794,
     &    647,   613,   548,   405,   953,   405,   953,   123,   918/
      DATA P( 9), ( C( 9,I), I = 1, 99 ) /   3079,
     &   1189,  1423,   287,   186,   341,    77,   733,   733,  1116,
     &      2,  1539,     2,     2,     2,     2,     2,  1116,   847,
     &   1174,     2,   827,   713,   910,   944,   139,  1174,  1174,
     &   1539,  1397,  1397,  1174,   370,    33,  1210,     2,   370,
     &   1423,   370,   370,  1423,  1423,  1423,   434,  1423,   901,
     &    139,  1174,   427,   427,   200,  1247,   114,   114,  1441,
     &    139,   728,  1116,  1174,   139,   113,   113,   113,  1406,
     &   1247,   200,   200,   200,   200,  1247,  1247,    27,   427,
     &    427,  1122,  1122,   696,   696,   427,  1539,   435,  1122,
     &    758,  1247,  1247,  1247,   200,   200,   200,  1247,   114,
     &     27,   118,   118,   113,   118,   453,   453,  1084,  1406/
      DATA P(10), ( C(10,I), I = 1, 99 ) /   4621,
     &   1764,  1349,  1859,   693,    78,   438,   531,    68,  2234,
     &   2310,  2310,  2310,     2,  2310,  2310,  2102,  2102,   178,
     &    314,   921,  1074,  1074,  1074,  2147,   314,  1869,   178,
     &    178,  1324,  1324,   510,  2309,  1541,  1541,  1541,  1541,
     &    342,  1324,  1324,  1324,  1324,   510,   570,   570,  2197,
     &    173,  1202,   998,  1324,  1324,   178,  1324,  1324,  1541,
     &   1541,  1541,   342,  1541,   886,   178,  1324,  1324,  1324,
     &    510,   784,   784,   501,   652,  1541,  1541,  1324,   178,
     &   1324,   178,  1324,  1541,   342,  1541,  2144,   784,  2132,
     &   1324,  1324,  1324,  1324,   510,   652,  1804,  1541,  1541,
     &   1541,  2132,  1324,  1324,  1324,   178,   510,  1541,   652/
      DATA P(11), ( C(11,I), I = 1, 99 ) /   6947,
     &   2872,  1238,   387,  2135,   235,  1565,   221,  1515,  2950,
     &    486,  3473,     2,  2950,   982,  2950,  3122,  2950,  3172,
     &   2091,  2091,     9,  3449,  3122,  2846,  3122,  3122,  1947,
     &   2846,  3122,   772,  1387,  2895,  1387,     3,     3,     3,
     &   1320,  1320,  2963,  2963,  1320,  1320,  2380,   108,  1284,
     &    702,  1429,   907,  3220,  3125,  1320,  2963,  1320,  1320,
     &   2963,  1320,  1639,  3168,  1660,  2895,  2895,  2895,  2895,
     &   1639,  1297,  1639,   404,  3168,  2963,  2943,  2943,   550,
     &   1387,  1387,  2895,  2895,  2895,  1387,  2895,  1387,  2895,
     &   1320,  1320,  2963,  1320,  1320,  1320,  2963,  1320,     2,
     &   3473,     2,  3473,   772,  2550,     9,  1320,  2963,  1320/
      DATA P(12), ( C(12,I), I = 1, 99 ) /  10427,
     &   4309,  2339,  4154,  4480,  4967,   630,  5212,  2592,  4715,
     &   1808,  1808,  5213,     2,   216,  4014,  3499,  3499,  4204,
     &   2701,  2701,  5213,  4157,  1209,  4157,  4460,   335,  4460,
     &   1533,  4575,  4013,  4460,  1881,  2701,  4030,  4030,  1881,
     &   4030,  1738,   249,   335,    57,  2561,  2561,  2561,  1533,
     &   1533,  1533,  4013,  4013,  4013,  4013,  4013,  1533,   856,
     &    856,   468,   468,   468,  2561,   468,  2022,  2022,  2434,
     &    138,  4605,  1100,  2561,  2561,    57,    57,  3249,   468,
     &    468,   468,    57,   468,  1738,   313,   856,     6,  3877,
     &    468,   557,   468,    57,   468,  4605,  2022,     2,  4605,
     &    138,  1100,    57,  2561,    57,    57,  2022,  5213,  3249/
      DATA P(13), ( C(13,I), I = 1, 99 ) /  15641,
     &   6610,  1658,  3022,  2603,  5211,   265,  4985,     3,  4971,
     &   2127,  1877,  1877,     2,  2925,  3175,  3878,  1940,  1940,
     &   1940,  5117,  5117,  5771,  5117,  5117,  5117,  5117,  5117,
     &   5771,  5771,  5117,  3658,  3658,  3658,  3658,  3658,  3658,
     &   5255,  2925,  2619,  1714,  4100,  6718,  6718,  4100,  2322,
     &    842,  4100,  6718,  5119,  4728,  5255,  5771,  5771,  5771,
     &   5117,  5771,  5117,  5117,  5117,  5117,  5117,  5117,  5771,
     &   5771,  1868,  4483,  4728,  3658,  5255,  3658,  5255,  3658,
     &   3658,  5255,  5255,  3658,  6718,  6718,   842,  2322,  6718,
     &   4100,  6718,  4100,  4100,  5117,  5771,  5771,  5117,  5771,
     &   5771,  5771,  5771,  5117,  5117,  5117,  5771,  5771,  1868/
      DATA P(14), ( C(14,I), I = 1, 99 ) /  23473,
     &   9861,  7101,  6257,  7878, 11170, 11638,  7542,  2592,  2591,
     &   6074,  1428,  8925, 11736,  8925,  5623,  5623,  1535,  6759,
     &   9953,  9953, 11459,  9953,  7615,  7615, 11377, 11377,  2762,
     &  11734, 11459,  6892,  1535,  6759,  4695,  1535,  6892,     2,
     &      2,  6892,  6892,  4177,  4177,  6339,  6950,  1226,  1226,
     &   1226,  4177,  6892,  6890,  3640,  3640,  1226, 10590, 10590,
     &   6950,  6950,  6950,  1226,  6950,  6950,  7586,  7586,  7565,
     &   7565,  3640,  3640,  6950,  7565,  6950,  3599,  3599,  3599,
     &   2441,  4885,  4885,  4885,  7565,  7565,  1226,  1226,  1226,
     &   6950,  7586,  1346,  2441,  6339,  3640,  6950, 10590,  6339,
     &   6950,  6950,  6950,  1226,  1226,  6950,   836,  6891,  7565/
      DATA P(15), ( C(15,I), I = 1, 99 ) /  35221,
     &  13482,  5629,  6068, 11974,  4732, 14946, 12097, 17609, 11740,
     &  15170, 10478, 10478, 17610,     2,     2,  7064,  7064,  7064,
     &   5665,  1771,  2947,  4453, 12323, 17610, 14809, 14809,  5665,
     &   5665,  2947,  2947,  2947,  2947, 12323, 12323,  4453,  4453,
     &   2026, 11772,  2026, 11665, 12323, 12323,  3582,  2940,  2940,
     &   6654,  4449,  9254, 11470,   304,   304, 11470,   304, 11470,
     &   6156,  9254, 11772,  6654, 11772,  6156, 11470, 11470, 11772,
     &  11772, 11772, 11470, 11470,   304, 11470, 11470,   304, 11470,
     &    304, 11470,   304,   304,   304,  6654, 11508,   304,   304,
     &   6156,  3582, 11470, 11470, 11470, 17274,  6654,  6654,  6744,
     &   6711,  6654,  6156,  3370,  6654, 12134,  3370,  6654,  3582/
      DATA P(16), ( C(16,I), I = 1, 99 ) /  52837,
     &  13482,  5629,  6068, 11974,  4732, 14946, 12097, 17609, 11740,
     &  15170, 10478, 10478, 17610,     2,     2,  7064,  7064,  7064,
     &   5665,  1771,  2947,  4453, 12323, 17610, 14809, 14809,  5665,
     &   5665,  2947,  2947,  2947,  2947, 12323, 12323,  4453,  4453,
     &   2026, 11772,  2026, 11665, 12323, 12323,  3582,  2940,  2940,
     &   6654,  4449,  9254, 11470,   304,   304, 11470,   304, 11470,
     &   6156,  9254, 11772,  6654, 11772,  6156, 11470, 11470, 11772,
     &  11772, 11772, 11470, 11470,   304, 11470, 11470,   304, 11470,
     &    304, 11470,   304,   304,   304,  6654, 11508,   304,   304,
     &   6156,  3582, 11470, 11470, 11470, 17274,  6654,  6654,  6744,
     &   6711,  6654,  6156,  3370,  6654, 12134,  3370,  6654,  3582/
      DATA P(17), ( C(17,I), I = 1, 99 ) /  79259,
     &  34566, 38838, 23965, 17279, 35325, 33471,   330, 36050, 26419,
     &   3012, 38428, 36430, 36430, 36755, 39629,  5749,  5749, 36755,
     &   5749, 14353, 14353, 14353, 32395, 32395, 32395, 32395, 32396,
     &  32396, 32396, 32396, 27739, 14353, 36430, 36430, 36430, 15727,
     &  38428, 28987, 28987, 27739, 38428, 27739, 18786, 14353, 15727,
     &  28987, 19151, 19757, 19757, 19757, 14353, 22876, 19151, 24737,
     &  24737,  4412, 30567, 30537, 19757, 30537, 19757, 30537, 30537,
     &   4412, 24737, 28987, 19757, 19757, 19757, 30537, 30537, 33186,
     &   4010,  4010,  4010, 17307, 15217, 32789, 37709,  4010,  4010,
     &   4010, 33186, 33186,  4010, 11057, 39388, 33186,  1122, 15089,
     &  39629,     2,     2, 23899, 16466, 16466, 17038,  9477,  9260/
      DATA P(18), ( C(18,I), I = 1, 99 ) / 118891,
     &  31929, 40295,  2610,  5177, 17271, 23770,  9140,   952, 39631,
     &      3, 11424, 49719, 38267, 25172,     2,     2, 59445,     2,
     &  59445, 38267, 44358, 14673, 53892, 14674, 14673, 14674, 41368,
     &  17875, 17875, 30190, 20444, 55869, 15644, 25499, 15644, 20983,
     &  44358, 15644, 15644,   485, 41428,   485,   485,   485, 41428,
     &  53798, 50230, 53798, 50253, 50253, 35677, 35677, 17474,  7592,
     &   4098, 17474,   485, 41428,   485, 41428,   485, 41428,   485,
     &  41428, 41428, 41428, 41428, 41428,  9020, 22816,  4098,  4098,
     &   4098,  7592, 42517,   485, 50006, 50006, 22816, 22816,  9020,
     &    485, 41428, 41428, 41428, 41428, 50006,   485, 41428, 41428,
     &  41428, 41428, 22816, 41428, 41428,   485,   485,   485,  9020/
      DATA P(19), ( C(19,I), I = 1, 99 ) / 178349,
     &  73726, 16352, 16297, 74268, 60788,  8555,  1077, 25486, 86595,
     &  59450, 19958, 62205, 62205,  4825,  4825, 89174, 89174, 62205,
     &  19958, 62205, 19958, 27626, 63080, 62205, 62205, 62205, 19958,
     &   8914, 83856, 30760, 47774, 47774, 19958, 62205, 39865, 39865,
     &  74988, 75715, 75715, 74988, 34522, 74988, 74988, 25101, 44621,
     &  44621, 44621, 25101, 25101, 25101, 44621, 47768, 41547, 44621,
     &  10273, 74988, 74988, 74988, 74988, 74988, 74988, 34522, 34522,
     &  67796, 67796, 30208,     2, 67062, 18500, 29251, 29251,     2,
     &  67796, 67062, 38649, 59302,  6225, 67062,  6475,  6225, 46772,
     &  38649, 67062, 46772, 46772, 67062, 46772, 25372, 67062,  6475,
     &  25372, 67062, 67062, 67062,  6225, 67062, 67062, 68247, 80676/
      DATA P(20), ( C(20,I), I = 1, 99 )/ 267523,
     & 103650, 50089, 70223, 41805, 74847,112775, 40889, 64866, 44053,
     &   1754,129471, 13630, 53467, 53467, 61378,133761,     2,133761,
     &      2,133761,133761, 65531, 65531, 65531, 38080,133761,133761,
     & 131061,  5431, 65531, 78250, 11397, 38841, 38841,107233,107233,
     & 111286, 19065, 38841, 19065, 19065, 16099,127638, 82411, 96659,
     &  96659, 82411, 96659, 82411, 51986,101677, 39264, 39264,101677,
     &  39264, 39264, 47996, 96659, 82411, 47996, 10971, 10004, 82411,
     &  96659, 82411, 82411, 82411, 96659, 96659, 96659, 82411, 96659,
     &  51986,110913, 51986, 51986,110913, 82411, 54713, 54713, 22360,
     & 117652, 22360, 78250, 78250, 91996, 22360, 91996, 97781, 91996,
     &  97781, 91996, 97781, 97781, 91996, 97781, 97781, 36249, 39779/
      END
*
      SUBROUTINE KROSUM( NDIM, SUMKRO, PRIME, VK,
     &                   FUNCTN, ALPHA, X, TID )
      EXTERNAL FUNCTN
      INTEGER NDIM, PRIME, K, J, TID
      DOUBLE PRECISION SUMKRO, VK(*), FUNCTN, ALPHA(*), X(*), ONE, UNI
      PARAMETER ( ONE = 1 )
      SUMKRO = 0
      DO J = 1, NDIM
         ALPHA(J) = UNI()
      END DO
      DO K = 1, PRIME
         DO J = 1, NDIM
            X(J) = MOD( K*VK(J) + ALPHA(J), ONE )
            X(J) = ABS( 2*X(J) - 1 )
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X,TID) - SUMKRO )/( 2*K - 1 )
         DO J = 1, NDIM
            X(J) = 1 - X(J)
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X,TID) - SUMKRO )/( 2*K )
      END DO
      END
*
      DOUBLE PRECISION FUNCTION MVNFNC(N, W, TID)
*
*     Integrand subroutine
*
      INTEGER N, INFIN(*), INFIS, TID
      DOUBLE PRECISION W(*), LOWER(*), UPPER(*), CORREL(*), ONE
      INTEGER NL, IJ, I, J
      PARAMETER ( NL = 100, NTHREADS = 64, ONE = 1 )
      DOUBLE PRECISION COV((NL*(NL+1))/2,NTHREADS), A(NL,NTHREADS) 
      DOUBLE PRECISION B(NL,NTHREADS), Y(NL), BVN
      INTEGER INFI(NL,NTHREADS)
      DOUBLE PRECISION PROD, D1(NTHREADS), E1(NTHREADS), DI, EI 
      DOUBLE PRECISION SUM, PHINV, D, E, MVNNIT
      SAVE D1, E1, A, B, INFI, COV
      DI = D1(TID)
      EI = E1(TID)
      PROD = E1(TID) - D1(TID)
      IJ = 1
      DO I = 1,N
         Y(I) = PHINV( DI + W(I)*(EI-DI) )
         SUM = 0
         DO J = 1,I
            IJ = IJ + 1
            SUM = SUM + COV(IJ,TID)*Y(J)
         END DO
         IJ = IJ + 1
         IF ( COV(IJ,TID) .GT. 0 ) THEN
            CALL LIMITS( A(I+1,TID)-SUM, B(I+1,TID)-SUM, 
     &                   INFI(I+1,TID), DI, EI )
         ELSE
            DI = ( 1 + SIGN( ONE, A(I+1,TID)-SUM ) )/2
            EI = ( 1 + SIGN( ONE, B(I+1,TID)-SUM ) )/2
         ENDIF
         PROD = PROD*(EI-DI)
      END DO
      MVNFNC = PROD
      RETURN
*
*     Entry point for intialization.
*
      ENTRY MVNNIT(N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E, TID)
      MVNNIT = 0
*
*     Initialization and computation of covariance Cholesky factor.
*
      CALL NCVSRT(N,LOWER,UPPER,CORREL,INFIN,Y,INFIS,
     &            A(:,TID),B(:,TID),INFI(:,TID),COV(:,TID),D,E)
      D1(TID) = D
      E1(TID) = E
      IF ( N - INFIS .EQ. 2 ) THEN
         D = SQRT( 1 + COV(2,TID)**2 )
         A(2,TID) = A(2,TID)/D
         B(2,TID) = B(2,TID)/D
         E = BVN(A(:,TID),B(:,TID),INFI(:,TID),COV(2,TID)/D)
         D = 0
         INFIS = INFIS + 1
      END IF
      END
      SUBROUTINE LIMITS( A, B, INFIN, LOWER, UPPER )
      DOUBLE PRECISION A, B, LOWER, UPPER, PHI
      INTEGER INFIN
      LOWER = 0
      UPPER = 1
      IF ( INFIN .GE. 0 ) THEN
         IF ( INFIN .NE. 0 ) LOWER = PHI(A)
         IF ( INFIN .NE. 1 ) UPPER = PHI(B)
      ENDIF
      END
      SUBROUTINE NCVSRT( N, LOWER, UPPER, CORREL, INFIN, Y, INFIS,
     &                   A, B, INFI, COV, D, E )
*
*     Subroutine to sort integration limits.
*
      INTEGER N, INFI(*), INFIN(*), INFIS
      DOUBLE PRECISION
     &     A(*), B(*), COV(*), LOWER(*), UPPER(*), CORREL(*), Y(*), D, E
      INTEGER I, J, K, IJ, II, JMIN
      DOUBLE PRECISION SUMSQ, ZERO
      PARAMETER ( ZERO = 0 )
      DOUBLE PRECISION AJ, BJ, SUM, SQTWPI
      DOUBLE PRECISION CVDIAG, AMIN, BMIN, DMIN, EMIN, YL, YU
      PARAMETER ( SQTWPI = 2.50662 82746 31000 50240 )
      IJ = 0
      II = 0
      INFIS = 0
      DO I = 1,N
         INFI(I) = INFIN(I)
         IF ( INFI(I) .LT. 0 ) THEN
            INFIS = INFIS + 1
         ELSE
            A(I) = 0
            B(I) = 0
            IF ( INFI(I) .NE. 0 ) A(I) = LOWER(I)
            IF ( INFI(I) .NE. 1 ) B(I) = UPPER(I)
         ENDIF
         DO J = 1,I-1
            IJ = IJ + 1
            II = II + 1
            COV(IJ) = CORREL(II)
         END DO
         IJ = IJ + 1
         COV(IJ) = 1
      END DO
*
*     First move any doubly infinite limits to innermost positions
*
      IF ( INFIS .LT. N ) THEN
         outer: DO I = N,N-INFIS+1,-1
            IF ( INFI(I) .GE. 0 ) THEN
               inner: DO J = 1,I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     CALL RCSWAP(J, I, A, B, INFI, N, COV)
                     CYCLE outer
                  ENDIF
               END DO inner
            ENDIF
         END DO outer
*
*     Sort remaining limits and determine Cholesky decomposition
*
         II = 0
         DO I = 1,N-INFIS
*
*     Determine the integration limits for variable with minimum
*      expected probability and interchange that variable with Ith.
*
            EMIN = 1
            DMIN = 0
            JMIN = I
            CVDIAG = 0
            IJ = II
            DO J = I, N-INFIS
               SUM = 0
               SUMSQ = 0
               DO K = 1, I-1
                  SUM = SUM + COV(IJ+K)*Y(K)
                  SUMSQ = SUMSQ + COV(IJ+K)**2
               END DO
               IJ = IJ + J
               SUMSQ = SQRT( MAX( COV(IJ)-SUMSQ, ZERO ) )
               IF ( SUMSQ .GT. 0 ) THEN
                  IF ( INFI(J) .NE. 0 ) AJ = ( A(J) - SUM )/SUMSQ
                  IF ( INFI(J) .NE. 1 ) BJ = ( B(J) - SUM )/SUMSQ
                  CALL LIMITS( AJ, BJ, INFI(J), D, E )
                  IF ( EMIN - DMIN .GE. E - D ) THEN
                     JMIN = J
                     IF ( INFI(J) .NE. 0 ) AMIN = AJ
                     IF ( INFI(J) .NE. 1 ) BMIN = BJ
                     DMIN = D
                     EMIN = E
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
            END DO
            IF ( JMIN .NE. I) CALL RCSWAP(I, JMIN, A, B, INFI, N, COV)
*
*     Compute Ith column of Cholesky factor.
*
            IJ = II + I
            COV(IJ) = CVDIAG
            DO J = I+1, N-INFIS
               IF ( CVDIAG .GT. 0 ) THEN
                  SUM = COV(IJ+I)
                  DO K = 1, I-1
                     SUM = SUM - COV(II+K)*COV(IJ+K)
                  END DO
                  COV(IJ+I) = SUM/CVDIAG
               ELSE
                  COV(IJ+I) = 0
               ENDIF
               IJ = IJ + J
            END DO
*
*     Compute expected value for Ith integration variable and
*     scale Ith covariance matrix row and limits.
*
            IF ( CVDIAG .GT. 0 ) THEN
               IF ( EMIN .GT. DMIN + 1D-8 ) THEN
                  YL = 0
                  YU = 0
                  IF ( INFI(I) .NE. 0 ) YL = -EXP( -AMIN**2/2 )/SQTWPI
                  IF ( INFI(I) .NE. 1 ) YU = -EXP( -BMIN**2/2 )/SQTWPI
                  Y(I) = ( YU - YL )/( EMIN - DMIN )
               ELSE
                  IF ( INFI(I) .EQ. 0 ) Y(I) = BMIN
                  IF ( INFI(I) .EQ. 1 ) Y(I) = AMIN
                  IF ( INFI(I) .EQ. 2 ) Y(I) = ( AMIN + BMIN )/2
               END IF
               DO J = 1,I
                  II = II + 1
                  COV(II) = COV(II)/CVDIAG
               END DO
               IF ( INFI(I) .NE. 0 ) A(I) = A(I)/CVDIAG
               IF ( INFI(I) .NE. 1 ) B(I) = B(I)/CVDIAG
            ELSE
               Y(I) = 0
               II = II + I
            ENDIF
         END DO
         CALL LIMITS( A(1), B(1), INFI(1), D, E)
      ENDIF
      END
      DOUBLE PRECISION FUNCTION CONDIT( N, SYMIN )
*
*     Computes condition number of symmetric matix in situ
*
      INTEGER NL, N
      PARAMETER ( NL = 100 )
      DOUBLE PRECISION DET, SYMIN(*), SUM, ROWMX, ROWMXI,
     & SYM(NL*(NL+1)/2)
      INTEGER II, IJ, I, J, IM
      ROWMX = 0
      IJ = 0
      DO I = 1,N
         SUM = 0
         IM = (I-2)*(I-1)/2
         DO J = 1,I-1
            IM = IM + 1
            SUM = SUM + ABS(SYMIN(IM))
            IJ = IJ + 1
            SYM(IJ) = SYMIN(IM)
         END DO
         SUM = SUM + 1
         IJ = IJ + 1
         SYM(IJ) = 1
         IM = IM + I
         DO J = I,N-1
            SUM = SUM + ABS(SYMIN(IM))
            IM = IM + J
         END DO
         ROWMX = MAX( SUM, ROWMX )
      END DO
      CALL SYMINV2(N, SYM, DET)
      ROWMXI = 0
      II = 0
      DO I = 1,N
         SUM = 0
         IJ = II
         DO J = 1,I
            IJ = IJ + 1
            SUM = SUM + ABS(SYM(IJ))
         END DO
         DO J = I,N-1
            IJ = IJ + J
            SUM = SUM + ABS(SYM(IJ))
         END DO
         ROWMXI = MAX( SUM, ROWMXI )
         II = II + I
      END DO
      CONDIT = ROWMX*ROWMXI
      END
      SUBROUTINE SYMINV2(N, LOWINV, DET)
*
*     Computes lower symmetric inverse and determinant in situ
*
      INTEGER I, II, N
      DOUBLE PRECISION LOWINV(*), DET
      CALL CHOLSK(N, LOWINV)
      DET = 1
      II = 0
      DO I = 1,N
         II = II + I
         DET = DET*LOWINV(II)
      END DO
      DET = DET*DET
      CALL CHOLNV(N, LOWINV)
      CALL CHOLPI(N, LOWINV)
      END
      SUBROUTINE CHOLSK(N, CHOFAC)
*
*     Computes Choleski factor in situ
*
      INTEGER I, II, J, JJ, K, N
      DOUBLE PRECISION CHOFAC(*), T
      DOUBLE PRECISION S, ZERO
      PARAMETER ( ZERO = 0 )
      JJ = 0
      DO J = 1,N
         II = JJ
         DO I = J,N
            S = CHOFAC(II+J)
            DO K = 1,J-1
               S = S - CHOFAC(II+K)*CHOFAC(JJ+K)
            END DO
            IF ( I .EQ. J ) THEN
               T = SQRT( MAX( S, ZERO ) )
               CHOFAC(II+J) = T
            ELSE
               CHOFAC(II+J) = S/T
            ENDIF
            II = II + I
         END DO
         JJ = JJ + J
      END DO
      END
      SUBROUTINE CHOLNV(N, CHOINV)
*
*     Inverts a lower triangular matrix in situ
*
      INTEGER I, II, J, JJ, K, KK, N
      DOUBLE PRECISION CHOINV(*), T
      DOUBLE PRECISION S
      II = 0
      DO I = 1,N
         T = 1/CHOINV(II+I)
         JJ = 0
         DO J = 1,I-1
            S = 0
            JJ = JJ + J
            KK = JJ
            DO K = J,I-1
               S = S + CHOINV(II+K)*CHOINV(KK)
               KK = KK + K
            END DO
            CHOINV(II+J) = -S*T
         END DO
         II = II + I
         CHOINV(II) = T
      END DO
      END
      SUBROUTINE CHOLPI(N, CHOPDI)
*
*     Multiplies Choleski inverse factors in situ
*
      INTEGER I, II, J, JJ, K, KK, N
      DOUBLE PRECISION CHOPDI(*)
      DOUBLE PRECISION S
      II = 0
      DO I = 1,N
         DO J = 1,I
            S = 0
            JJ = II + I
            KK = II + J
            DO K = I,N
               S = S + CHOPDI(KK)*CHOPDI(JJ)
               JJ = JJ + K
               KK = KK + K
            END DO
            CHOPDI(II+J) = S
         END DO
         II = II + I
      END DO
      END
      SUBROUTINE RCSWAP(P, Q, A, B, INFIN, N, C)
*
*     Swaps rows and columns P and Q in situ.
*
      DOUBLE PRECISION A(*), B(*), C(*), T
      INTEGER INFIN(*), P, Q, N, I, J, II, JJ
      T = A(P)
      A(P) = A(Q)
      A(Q) = T
      T = B(P)
      B(P) = B(Q)
      B(Q) = T
      J = INFIN(P)
      INFIN(P) = INFIN(Q)
      INFIN(Q) = J
      JJ = (P*(P-1))/2
      II = (Q*(Q-1))/2
      T = C(JJ+P)
      C(JJ+P) = C(II+Q)
      C(II+Q) = T
      DO J = 1, P-1
         T = C(JJ+J)
         C(JJ+J) = C(II+J)
         C(II+J) = T
      END DO
      JJ = JJ + P
      DO I = P+1, Q-1
         T = C(JJ+P)
         C(JJ+P) = C(II+I)
         C(II+I) = T
         JJ = JJ + I
      END DO
      II = II + Q
      DO I = Q+1, N
         T = C(II+P)
         C(II+P) = C(II+Q)
         C(II+Q) = T
         II = II + I
      END DO
      END
      DOUBLE PRECISION FUNCTION PHI(Z)
*
*     Normal distribution probabilities accurate to 1.e-15.
*     Z = no. of standard deviations from the mean.
*
*     Based upon algorithm 5666 for the error function, from:
*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*
*     Programmer: Alan Miller
*
*     Latest revision - 30 March 1986
*
      DOUBLE PRECISION P0, P1, P2, P3, P4, P5, P6,
     &     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     &     Z, P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(
     &     P0 = 220.20 68679 12376 1D0,
     &     P1 = 221.21 35961 69931 1D0,
     &     P2 = 112.07 92914 97870 9D0,
     &     P3 = 33.912 86607 83830 0D0,
     &     P4 = 6.3739 62203 53165 0D0,
     &     P5 = .70038 30644 43688 1D0,
     &     P6 = .035262 49659 98910 9D0)
      PARAMETER(
     &     Q0 = 440.41 37358 24752 2D0,
     &     Q1 = 793.82 65125 19948 4D0,
     &     Q2 = 637.33 36333 78831 1D0,
     &     Q3 = 296.56 42487 79673 7D0,
     &     Q4 = 86.780 73220 29460 8D0,
     &     Q5 = 16.064 17757 92069 5D0,
     &     Q6 = 1.7556 67163 18264 2D0,
     &     Q7 = .088388 34764 83184 4D0)
      PARAMETER(ROOTPI = 2.5066 28274 63100 1D0)
      PARAMETER(CUTOFF = 7.0710 67811 86547 5D0)
*
      ZABS = ABS(Z)
*
*     |Z| > 37
*
      IF (ZABS .GT. 37) THEN
         P = 0
      ELSE
*
*     |Z| <= 37
*
         EXPNTL = EXP(-ZABS**2/2)
*
*     |Z| < CUTOFF = 10/SQRT(2)
*
         IF (ZABS .LT. CUTOFF) THEN
            P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS
     &           + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS
     &           + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS
     &           + Q0)
*
*     |Z| >= CUTOFF.
*
         ELSE
            P = EXPNTL/(ZABS + 1/(ZABS + 2/(ZABS + 3/(ZABS + 4/
     &           (ZABS + 0.65D0)))))/ROOTPI
         END IF
      END IF
      IF (Z .GT. 0) P = 1 - P
      PHI = P
      END
      DOUBLE PRECISION FUNCTION PHINV(P)
*
*	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*
*	Produces the normal deviate Z corresponding to a given lower
*	tail area of P.
*
*	The hash sums below are the sums of the mantissas of the
*	coefficients.   They are included for use in checking
*	transcription.
*
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2,
     &     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7,
     &     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7,
     &     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7,
     &     P, Q, R
      PARAMETER (SPLIT1 = 0.425, SPLIT2 = 5,
     &     CONST1 = 0.180625D0, CONST2 = 1.6D0)
*
*     Coefficients for P close to 0.5
*
      PARAMETER (
     &     A0 = 3.38713 28727 96366 6080D0,
     &     A1 = 1.33141 66789 17843 7745D+2,
     &     A2 = 1.97159 09503 06551 4427D+3,
     &     A3 = 1.37316 93765 50946 1125D+4,
     &     A4 = 4.59219 53931 54987 1457D+4,
     &     A5 = 6.72657 70927 00870 0853D+4,
     &     A6 = 3.34305 75583 58812 8105D+4,
     &     A7 = 2.50908 09287 30122 6727D+3,
     &     B1 = 4.23133 30701 60091 1252D+1,
     &     B2 = 6.87187 00749 20579 0830D+2,
     &     B3 = 5.39419 60214 24751 1077D+3,
     &     B4 = 2.12137 94301 58659 5867D+4,
     &     B5 = 3.93078 95800 09271 0610D+4,
     &     B6 = 2.87290 85735 72194 2674D+4,
     &     B7 = 5.22649 52788 52854 5610D+3)
*     HASH SUM AB    55.88319 28806 14901 4439
*
*     Coefficients for P not close to 0, 0.5 or 1.
*
      PARAMETER (
     &     C0 = 1.42343 71107 49683 57734D0,
     &     C1 = 4.63033 78461 56545 29590D0,
     &     C2 = 5.76949 72214 60691 40550D0,
     &     C3 = 3.64784 83247 63204 60504D0,
     &     C4 = 1.27045 82524 52368 38258D0,
     &     C5 = 2.41780 72517 74506 11770D-1,
     &     C6 = 2.27238 44989 26918 45833D-2,
     &     C7 = 7.74545 01427 83414 07640D-4,
     &     D1 = 2.05319 16266 37758 82187D0,
     &     D2 = 1.67638 48301 83803 84940D0,
     &     D3 = 6.89767 33498 51000 04550D-1,
     &     D4 = 1.48103 97642 74800 74590D-1,
     &     D5 = 1.51986 66563 61645 71966D-2,
     &     D6 = 5.47593 80849 95344 94600D-4,
     &     D7 = 1.05075 00716 44416 84324D-9)
*     HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (
     &     E0 = 6.65790 46435 01103 77720D0,
     &     E1 = 5.46378 49111 64114 36990D0,
     &     E2 = 1.78482 65399 17291 33580D0,
     &     E3 = 2.96560 57182 85048 91230D-1,
     &     E4 = 2.65321 89526 57612 30930D-2,
     &     E5 = 1.24266 09473 88078 43860D-3,
     &     E6 = 2.71155 55687 43487 57815D-5,
     &     E7 = 2.01033 43992 92288 13265D-7,
     &     F1 = 5.99832 20655 58879 37690D-1,
     &     F2 = 1.36929 88092 27358 05310D-1,
     &     F3 = 1.48753 61290 85061 48525D-2,
     &     F4 = 7.86869 13114 56132 59100D-4,
     &     F5 = 1.84631 83175 10054 68180D-5,
     &     F6 = 1.42151 17583 16445 88870D-7,
     &     F7 = 2.04426 31033 89939 78564D-15)
*     HASH SUM EF    47.52583 31754 92896 71629
*
      Q = ( 2*P - 1 )/2
      IF ( ABS(Q) .LE. SPLIT1 ) THEN
         R = CONST1 - Q*Q
         PHINV = Q*(((((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     &        *R + A2)*R + A1)*R + A0) /
     &        (((((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     &        *R + B2)*R + B1)*R + 1)
      ELSE
         R = MIN( P, 1 - P )
         IF (R .GT. 0) THEN
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               PHINV = (((((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     &              *R + C2)*R + C1)*R + C0) /
     &              (((((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     &              *R + D2)*R + D1)*R + 1)
            ELSE
               R = R - SPLIT2
               PHINV = (((((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     &              *R + E2)*R + E1)*R + E0) /
     &              (((((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     &              *R + F2)*R + F1)*R + 1)
            END IF
         ELSE
            PHINV = 9
         END IF
         IF ( Q .LT. 0 ) PHINV = - PHINV
      END IF
      END
      DOUBLE PRECISION FUNCTION BVN ( LOWER, UPPER, INFIN, CORREL )
*
*     A function for computing bivariate normal probabilities.
*
*  Parameters
*
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, correlation coefficient.
*
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL, BVNU
      INTEGER INFIN(*)
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
     +        - BVNU ( UPPER(1), LOWER(2), CORREL )
     +        - BVNU ( LOWER(1), UPPER(2), CORREL )
     +        + BVNU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
     +        - BVNU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
     +        - BVNU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         BVN =  BVNU ( -UPPER(1), -UPPER(2), CORREL )
     +        - BVNU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         BVN =  BVNU ( -UPPER(1), -UPPER(2), CORREL )
     +        - BVNU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         BVN =  BVNU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         BVN =  BVNU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         BVN =  BVNU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      END
      DOUBLE PRECISION FUNCTION BVNU( SH, SK, R )
*
*     A function for computing bivariate normal probabilities.
*
*       Yihong Ge
*       Department of Computer Science and Electrical Engineering
*       Washington State University
*       Pullman, WA 99164-2752
*       Email : yge@eecs.wsu.edu
*     and
*       Alan Genz
*       Department of Mathematics
*       Washington State University
*       Pullman, WA 99164-3113
*       Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK.
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION BVN, SH, SK, R, ZERO, TWOPI
      INTEGER I, LG, NG
      PARAMETER ( ZERO = 0, TWOPI = 6.2831 85307 179586 )
      DOUBLE PRECISION X(10,3), W(10,3), AS, A, B, C, D, RS, XS
      DOUBLE PRECISION PHI, SN, ASR, H, K, BS, HS, HK
*     Gauss Legendre Points and Weights, N =  6
      DATA ( W(I,1), X(I,1), I = 1,3) /
     &  0.1713244923791705D+00,-0.9324695142031522D+00,
     &  0.3607615730481384D+00,-0.6612093864662647D+00,
     &  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     &  0.4717533638651177D-01,-0.9815606342467191D+00,
     &  0.1069393259953183D+00,-0.9041172563704750D+00,
     &  0.1600783285433464D+00,-0.7699026741943050D+00,
     &  0.2031674267230659D+00,-0.5873179542866171D+00,
     &  0.2334925365383547D+00,-0.3678314989981802D+00,
     &  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1,10) /
     &  0.1761400713915212D-01,-0.9931285991850949D+00,
     &  0.4060142980038694D-01,-0.9639719272779138D+00,
     &  0.6267204833410906D-01,-0.9122344282513259D+00,
     &  0.8327674157670475D-01,-0.8391169718222188D+00,
     &  0.1019301198172404D+00,-0.7463319064601508D+00,
     &  0.1181945319615184D+00,-0.6360536807265150D+00,
     &  0.1316886384491766D+00,-0.5108670019508271D+00,
     &  0.1420961093183821D+00,-0.3737060887154196D+00,
     &  0.1491729864726037D+00,-0.2277858511416451D+00,
     &  0.1527533871307259D+00,-0.7652652113349733D-01/
      SAVE X, W
      IF ( ABS(R) .LT. 0.3 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75 ) THEN
         NG = 2
         LG = 6
      ELSE
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK
      HK = H*K
      BVN = 0
      IF ( ABS(R) .LT. 0.925 ) THEN
         HS = ( H*H + K*K )/2
         ASR = ASIN(R)
         DO 10  I = 1, LG
            SN = SIN(ASR*( X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
            SN = SIN(ASR*(-X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
 10      CONTINUE
         BVN = BVN*ASR/(2*TWOPI) + PHI(-H)*PHI(-K)
      ELSE
         IF ( R .LT. 0 ) THEN
            K = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. 1 ) THEN
            AS = ( 1 - R )*( 1 + R )
            A = SQRT(AS)
            BS = ( H - K )**2
            C = ( 4 - HK )/8
            D = ( 12 - HK )/16
            BVN = A*EXP( -(BS/AS + HK)/2 )
     +             *( 1 - C*(BS - AS)*(1 - D*BS/5)/3 + C*D*AS*AS/5 )
            IF ( HK .GT. -160 ) THEN
               B = SQRT(BS)
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*PHI(-B/A)*B
     +                    *( 1 - C*BS*( 1 - D*BS/5 )/3 )
            ENDIF
            A = A/2
            DO 20 I = 1, LG
               XS = ( A*(X(I,NG)+1) )**2
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*
     +              ( EXP( -BS/(2*XS) - HK/(1+RS) )/RS
     +              - EXP( -(BS/XS+HK)/2 )*( 1 + C*XS*( 1 + D*XS ) ) )
               XS = AS*(-X(I,NG)+1)**2/4
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*EXP( -(BS/XS + HK)/2 )
     +                    *( EXP( -HK*(1-RS)/(2*(1+RS)) )/RS
     +                       - ( 1 + C*XS*( 1 + D*XS ) ) )
 20         CONTINUE
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R .GT. 0 ) BVN =  BVN + PHI( -MAX( H, K ) )
         IF ( R .LT. 0 ) BVN = -BVN + MAX( ZERO, PHI(-H) - PHI(-K) )
      ENDIF
      BVNU = BVN
      END
*
      SUBROUTINE ADAPT(NDIM, MINCLS, MAXCLS, FUNCTN,
     &     ABSREQ, RELREQ, LENWRK, NTHREADS, WORK, 
     &     ABSEST, FINEST, INFORM, TID)
*
*   Adaptive Multidimensional Integration Subroutine
*
*   Author: Alan Genz
*           Department of Mathematics
*           Washington State University
*           Pullman, WA 99164-3113 USA
*
*  This subroutine computes an approximation to the integral
*
*      1 1     1
*     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1)
*      0 0     0
*
***************  Parameters for ADAPT  ********************************
*
****** Input Parameters
*
*  NDIM    Integer number of integration variables.
*  MINCLS  Integer minimum number of FUNCTN calls to be allowed; MINCLS
*          must not exceed MAXCLS. If MINCLS < 0, then ADAPT assumes
*          that a previous call of ADAPT has been made with the same
*          integrand and continues that calculation.
*  MAXCLS  Integer maximum number of FUNCTN calls to be used; MAXCLS
*          must be >= RULCLS, the number of function calls required for
*          one application of the basic integration rule.
*           IF ( NDIM .EQ. 1 ) THEN
*              RULCLS = 11
*           ELSE IF ( NDIM .LT. 15 ) THEN
*              RULCLS = 2**NDIM + 2*NDIM*(NDIM+3) + 1
*           ELSE
*              RULCLS = 1 + NDIM*(24-NDIM*(6-NDIM*4))/3
*           ENDIF
*  FUNCTN  Externally declared real user defined integrand. Its
*          parameters must be (NDIM, Z), where Z is a real array of
*          length NDIM.
*  ABSREQ  Real required absolute accuracy.
*  RELREQ  Real required relative accuracy.
*  LENWRK  Integer length of real array WORK (working storage); ADAPT
*          needs LENWRK >= 16*NDIM + 27. For maximum efficiency LENWRK
*          should be about 2*NDIM*MAXCLS/RULCLS if MAXCLS FUNCTN
*          calls are needed. If LENWRK is significantly less than this,
*          ADAPT may be less efficient.
*
****** Output Parameters
*
*  MINCLS  Actual number of FUNCTN calls used by ADAPT.
*  WORK    Real array (length LENWRK) of working storage. This contains
*          information that is needed for additional calls of ADAPT
*          using the same integrand (input MINCLS < 0).
*  ABSEST  Real estimated absolute accuracy.
*  FINEST  Real estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when ABSEST <= ABSREQ or
*                     ABSEST <= |FINEST|*RELREQ with MINCLS <= MAXCLS.
*          INFORM = 1 if MAXCLS was too small for ADAPT to obtain the
*                     result FINEST to within the requested accuracy.
*          INFORM = 2 if MINCLS > MAXCLS, LENWRK < 16*NDIM + 27 or
*                     RULCLS > MAXCLS.
*
************************************************************************
*
*     Begin driver routine. This routine partitions the working storage
*      array and then calls the main subroutine ADBASE.
*
      EXTERNAL FUNCTN
      INTEGER NDIM, MINCLS, MAXCLS, LENWRK, NTHREADS, INFORM
      DOUBLE PRECISION
     &     FUNCTN, ABSREQ, RELREQ, WORK(LENWRK, NTHREADS), 
     &     ABSEST, FINEST
      INTEGER SBRGNS, MXRGNS, RULCLS, LENRUL,
     & INERRS, INVALS, INPTRS, INLWRS, INUPRS, INMSHS, INPNTS, INWGTS,
     & INLOWR, INUPPR, INWDTH, INMESH, INWORK, TID
      IF ( NDIM .EQ. 1 ) THEN
         LENRUL = 5
         RULCLS = 9
      ELSE IF ( NDIM .LT. 12 ) THEN
         LENRUL = 6
         RULCLS = 2**NDIM + 2*NDIM*(NDIM+2) + 1
      ELSE
         LENRUL = 6
         RULCLS = 1 + 2*NDIM*(1+2*NDIM)
      ENDIF
      IF ( LENWRK .GE. LENRUL*(NDIM+4) + 10*NDIM + 3 .AND.
     &     RULCLS. LE. MAXCLS .AND. MINCLS .LE. MAXCLS ) THEN
        MXRGNS = ( LENWRK - LENRUL*(NDIM+4) - 7*NDIM )/( 3*NDIM + 3 )
        INERRS = 1
        INVALS = INERRS + MXRGNS
        INPTRS = INVALS + MXRGNS
        INLWRS = INPTRS + MXRGNS
        INUPRS = INLWRS + MXRGNS*NDIM
        INMSHS = INUPRS + MXRGNS*NDIM
        INWGTS = INMSHS + MXRGNS*NDIM
        INPNTS = INWGTS + LENRUL*4
        INLOWR = INPNTS + LENRUL*NDIM
        INUPPR = INLOWR + NDIM
        INWDTH = INUPPR + NDIM
        INMESH = INWDTH + NDIM
        INWORK = INMESH + NDIM
        IF ( MINCLS .LT. 0 ) SBRGNS = WORK(LENWRK,TID)
        CALL ADBASE(NDIM, MINCLS, MAXCLS, FUNCTN, ABSREQ, RELREQ,
     &       ABSEST, FINEST, SBRGNS, MXRGNS, RULCLS, LENRUL,
     &       WORK(INERRS,TID), WORK(INVALS,TID), 
     &       WORK(INPTRS,TID), WORK(INLWRS,TID),
     &       WORK(INUPRS,TID), WORK(INMSHS,TID), 
     &       WORK(INWGTS,TID), WORK(INPNTS,TID),
     &       WORK(INLOWR,TID), WORK(INUPPR,TID), 
     &       WORK(INWDTH,TID), WORK(INMESH,TID),
     &       WORK(INWORK,TID), INFORM, TID)
        WORK(LENWRK,TID) = SBRGNS
       ELSE
        INFORM = 2
        MINCLS = RULCLS
      ENDIF
      END
      SUBROUTINE BSINIT(NDIM, W, LENRUL, G)
*
*     For initializing basic rule weights and symmetric sum parameters.
*
      INTEGER NDIM, LENRUL, RULPTS(6), I, J, NUMNUL, SDIM
      PARAMETER ( NUMNUL = 4, SDIM = 12 )
      DOUBLE PRECISION W(LENRUL,4), G(NDIM,LENRUL)
      DOUBLE PRECISION LAM1, LAM2, LAM3, LAMP, RULCON
*
*     The following code determines rule parameters and weights for a
*      degree 7 rule (W(1,1),...,W(5,1)), two degree 5 comparison rules
*      (W(1,2),...,W(5,2) and W(1,3),...,W(5,3)) and a degree 3
*      comparison rule (W(1,4),...W(5,4)).
*
*       If NDIM = 1, then LENRUL = 5 and total points = 9.
*       If NDIM < SDIM, then LENRUL = 6 and
*                      total points = 1+2*NDIM*(NDIM+2)+2**NDIM.
*       If NDIM > = SDIM, then LENRUL = 6 and
*                      total points = 1+2*NDIM*(1+2*NDIM).
*
      DO I = 1,LENRUL
         DO J = 1,NDIM
            G(J,I) = 0
         END DO
         DO J = 1,NUMNUL
            W(I,J) = 0
         END DO
      END DO
      RULPTS(5) = 2*NDIM*(NDIM-1)
      RULPTS(4) = 2*NDIM
      RULPTS(3) = 2*NDIM
      RULPTS(2) = 2*NDIM
      RULPTS(1) = 1
      LAMP = 0.85
      LAM3 = 0.4707
      LAM2 = 4/(15 - 5/LAM3)
      W(5,1) = ( 3 - 5*LAM3 )/( 180*(LAM2-LAM3)*LAM2**2 )
      IF ( NDIM .LT. SDIM ) THEN
         LAM1 = 8*LAM3*(31*LAM3-15)/( (3*LAM3-1)*(5*LAM3-3)*35 )
         W(LENRUL,1) = 1/(3*LAM3)**3/2**NDIM
      ELSE
         LAM1 = ( LAM3*(15 - 21*LAM2) + 35*(NDIM-1)*(LAM2-LAM3)/9 )
     &       /  ( LAM3*(21 - 35*LAM2) + 35*(NDIM-1)*(LAM2/LAM3-1)/9 )
         W(6,1) = 1/(4*(3*LAM3)**3)
      ENDIF
      W(3,1) = ( 15 - 21*(LAM3+LAM1) + 35*LAM3*LAM1 )
     &     /( 210*LAM2*(LAM2-LAM3)*(LAM2-LAM1) ) - 2*(NDIM-1)*W(5,1)
      W(2,1) = ( 15 - 21*(LAM3+LAM2) + 35*LAM3*LAM2 )
     &     /( 210*LAM1*(LAM1-LAM3)*(LAM1-LAM2) )
      IF ( NDIM .LT. SDIM ) THEN
         RULPTS(LENRUL) = 2**NDIM
         LAM3 = SQRT(LAM3)
         DO I = 1,NDIM
            G(I,LENRUL) = LAM3
         END DO
      ELSE
         W(6,1) = 1/(4*(3*LAM3)**3)
         RULPTS(6) = 2*NDIM*(NDIM-1)
         LAM3 = SQRT(LAM3)
         DO I = 1,2
            G(I,6) = LAM3
         END DO
      ENDIF
      IF ( NDIM .GT. 1 ) THEN
         W(5,2) = 1/(6*LAM2)**2
         W(5,3) = 1/(6*LAM2)**2
      ENDIF
      W(3,2) = ( 3 - 5*LAM1 )/( 30*LAM2*(LAM2-LAM1) )
     &     - 2*(NDIM-1)*W(5,2)
      W(2,2) = ( 3 - 5*LAM2 )/( 30*LAM1*(LAM1-LAM2) )
      W(4,3) = ( 3 - 5*LAM2 )/( 30*LAMP*(LAMP-LAM2) )
      W(3,3) = ( 3 - 5*LAMP )/( 30*LAM2*(LAM2-LAMP) )
     &     - 2*(NDIM-1)*W(5,3)
      W(2,4) = 1/(6*LAM1)
      LAMP = SQRT(LAMP)
      LAM2 = SQRT(LAM2)
      LAM1 = SQRT(LAM1)
      G(1,2) = LAM1
      G(1,3) = LAM2
      G(1,4) = LAMP
      IF ( NDIM .GT. 1 ) THEN
         G(1,5) = LAM2
         G(2,5) = LAM2
      ENDIF
      DO J = 1, NUMNUL
         W(1,J) = 1
         DO I = 2,LENRUL
            W(1,J) = W(1,J) - RULPTS(I)*W(I,J)
         END DO
      END DO
      RULCON = 2
      CALL RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      END
      SUBROUTINE RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      INTEGER LENRUL, NUMNUL, I, J, K, RULPTS(*)
      DOUBLE PRECISION ALPHA, NORMCF, NORMNL, W(LENRUL, *), RULCON
*
*     Compute orthonormalized null rules.
*
      NORMCF = 0
      DO I = 1,LENRUL
         NORMCF = NORMCF + RULPTS(I)*W(I,1)*W(I,1)
      END DO
      DO K = 2,NUMNUL
         DO I = 1,LENRUL
            W(I,K) = W(I,K) - W(I,1)
         END DO
         DO J = 2,K-1
            ALPHA = 0
            DO I = 1,LENRUL
               ALPHA = ALPHA + RULPTS(I)*W(I,J)*W(I,K)
            END DO
            ALPHA = -ALPHA/NORMCF
            DO I = 1,LENRUL
               W(I,K) = W(I,K) + ALPHA*W(I,J)
            END DO
         END DO
         NORMNL = 0
         DO I = 1,LENRUL
            NORMNL = NORMNL + RULPTS(I)*W(I,K)*W(I,K)
         END DO
         ALPHA = SQRT(NORMCF/NORMNL)
         DO I = 1,LENRUL
            W(I,K) = ALPHA*W(I,K)
         END DO
      END DO
      DO J = 2, NUMNUL
         DO I = 1,LENRUL
            W(I,J) = W(I,J)/RULCON
         END DO
      END DO
      END
      SUBROUTINE ADBASE(NDIM, MINCLS, MAXCLS, FUNCTN, ABSREQ, RELREQ,
     &     ABSEST, FINEST, SBRGNS, MXRGNS, RULCLS, LENRUL,
     &     ERRORS, VALUES, PONTRS, LOWERS,
     &     UPPERS, MESHES, WEGHTS, POINTS,
     &     LOWER, UPPER, WIDTH, MESH, WORK, INFORM, TID)
*
*        Main adaptive integration subroutine
*
      EXTERNAL FUNCTN
      INTEGER I, J, NDIM, MINCLS, MAXCLS, SBRGNS, MXRGNS,
     &     RULCLS, LENRUL, INFORM, NWRGNS, TID
      DOUBLE PRECISION FUNCTN, ABSREQ, RELREQ, ABSEST, FINEST,
     &     ERRORS(*), VALUES(*), PONTRS(*),
     &     LOWERS(NDIM,*), UPPERS(NDIM,*),
     &     MESHES(NDIM,*), WEGHTS(*), POINTS(*),
     &     LOWER(*), UPPER(*), WIDTH(*), MESH(*), WORK(*)
      INTEGER DIVAXN, TOP, RGNCLS, FUNCLS, DIFCLS
 
*
*     Initialization of subroutine
*
      INFORM = 2
      FUNCLS = 0
      CALL BSINIT(NDIM, WEGHTS, LENRUL, POINTS)
      IF ( MINCLS .GE. 0) THEN
*
*       When MINCLS >= 0 determine initial subdivision of the
*       integration region and apply basic rule to each subregion.
*
         SBRGNS = 0
         DO I = 1,NDIM
            LOWER(I) = 0
            MESH(I) = 1
            WIDTH(I) = 1/(2*MESH(I))
            UPPER(I) = 1
         END DO
         DIVAXN = 0
         RGNCLS = RULCLS
         NWRGNS = 1
 10      CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK, 
     &        WORK(NDIM+1), FUNCTN, DIVAXN, DIFCLS, TID)
         FUNCLS = FUNCLS + DIFCLS
         IF ( FUNCLS + RGNCLS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
     &        .LE. MINCLS ) THEN
            RGNCLS = RGNCLS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
            NWRGNS = NWRGNS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
            MESH(DIVAXN) = MESH(DIVAXN) + 1
            WIDTH(DIVAXN) = 1/( 2*MESH(DIVAXN) )
            GO TO 10
         ENDIF
         IF ( NWRGNS .LE. MXRGNS ) THEN
            DO I = 1,NDIM
               UPPER(I) = LOWER(I) + 2*WIDTH(I)
               MESH(I) = 1
            END DO
         ENDIF
*
*     Apply basic rule to subregions and store results in heap.
*
 20      SBRGNS = SBRGNS + 1
         CALL BASRUL(NDIM, LOWER, UPPER, WIDTH, FUNCTN,
     &        WEGHTS, LENRUL, POINTS, WORK, 
     &        WORK(NDIM+1), ERRORS(SBRGNS), 
     &        VALUES(SBRGNS), TID)
         CALL TRESTR(SBRGNS, SBRGNS, PONTRS, ERRORS)
         DO I = 1,NDIM
            LOWERS(I,SBRGNS) = LOWER(I)
            UPPERS(I,SBRGNS) = UPPER(I)
            MESHES(I,SBRGNS) = MESH(I)
         END DO
         DO I = 1,NDIM
            LOWER(I) = UPPER(I)
            UPPER(I) = LOWER(I) + 2*WIDTH(I)
            IF ( LOWER(I)+WIDTH(I) .LT. 1 )  GO TO 20
            LOWER(I) = 0
            UPPER(I) = LOWER(I) + 2*WIDTH(I)
         END DO
         FUNCLS = FUNCLS + SBRGNS*RULCLS
      ENDIF
*
*     Check for termination
*
 30   FINEST = 0
      ABSEST = 0
      DO I = 1, SBRGNS
         FINEST = FINEST + VALUES(I)
         ABSEST = ABSEST + ERRORS(I)
      END DO
      IF ( ABSEST .GT. MAX( ABSREQ, RELREQ*ABS(FINEST) )
     &     .OR. FUNCLS .LT. MINCLS ) THEN
*
*     Prepare to apply basic rule in (parts of) subregion with
*     largest error.
*
         TOP = PONTRS(1)
         RGNCLS = RULCLS
         DO I = 1,NDIM
            LOWER(I) = LOWERS(I,TOP)
            UPPER(I) = UPPERS(I,TOP)
            MESH(I) = MESHES(I,TOP)
            WIDTH(I) = (UPPER(I)-LOWER(I))/(2*MESH(I))
            RGNCLS = RGNCLS*MESH(I)
         END DO
         CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK, 
     &        WORK(NDIM+1), FUNCTN, DIVAXN, DIFCLS, TID)
         FUNCLS = FUNCLS + DIFCLS
         RGNCLS = RGNCLS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
         IF ( FUNCLS + RGNCLS .LE. MAXCLS ) THEN
            IF ( SBRGNS + 1 .LE. MXRGNS ) THEN
*
*     Prepare to subdivide into two pieces.
*
               NWRGNS = 1
               WIDTH(DIVAXN) = WIDTH(DIVAXN)/2
            ELSE
               NWRGNS = 0
               WIDTH(DIVAXN) = WIDTH(DIVAXN)
     &                        *MESH(DIVAXN)/( MESH(DIVAXN) + 1 )
               MESHES(DIVAXN,TOP) = MESH(DIVAXN) + 1
            ENDIF
            IF ( NWRGNS .GT. 0 ) THEN
*
*     Only allow local subdivision when space is available.
*
               DO J = SBRGNS+1,SBRGNS+NWRGNS
                  DO I = 1,NDIM
                     LOWERS(I,J) = LOWER(I)
                     UPPERS(I,J) = UPPER(I)
                     MESHES(I,J) = MESH(I)
                  END DO
               END DO
               UPPERS(DIVAXN,TOP) = LOWER(DIVAXN) + 2*WIDTH(DIVAXN)
               LOWERS(DIVAXN,SBRGNS+1) = UPPERS(DIVAXN,TOP)
            ENDIF
            FUNCLS = FUNCLS + RGNCLS
            CALL BASRUL(NDIM, LOWERS(1,TOP), UPPERS(1,TOP), WIDTH,
     &           FUNCTN, WEGHTS, LENRUL, POINTS, WORK, 
     &           WORK(NDIM+1), ERRORS(TOP), VALUES(TOP), TID)
            CALL TRESTR(TOP, SBRGNS, PONTRS, ERRORS)
            DO I = SBRGNS+1, SBRGNS+NWRGNS
*
*     Apply basic rule and store results in heap.
*
               CALL BASRUL(NDIM, LOWERS(1,I), UPPERS(1,I), WIDTH,
     &              FUNCTN, WEGHTS, LENRUL, POINTS, WORK, 
     &              WORK(NDIM+1), ERRORS(I), VALUES(I), TID)
               CALL TRESTR(I, I, PONTRS, ERRORS)
            END DO
            SBRGNS = SBRGNS + NWRGNS
            GO TO 30
         ELSE
            INFORM = 1
         ENDIF
      ELSE
         INFORM = 0
      ENDIF
      MINCLS = FUNCLS
      END
      SUBROUTINE BASRUL( NDIM, A, B, WIDTH, FUNCTN, W, LENRUL, G,
     &     CENTER, Z, RGNERT, BASEST, TID )
*
*     For application of basic integration rule
*
      EXTERNAL FUNCTN
      INTEGER I, LENRUL, NDIM, TID
      DOUBLE PRECISION
     &     A(NDIM), B(NDIM), WIDTH(NDIM), FUNCTN, W(LENRUL,4),
     &     G(NDIM,LENRUL), CENTER(NDIM), Z(NDIM), RGNERT, BASEST
      DOUBLE PRECISION
     &     FULSUM, FSYMSM, RGNCMP, RGNVAL, RGNVOL, RGNCPT, RGNERR
*
*     Compute Volume and Center of Subregion
*
      RGNVOL = 1
      DO I = 1,NDIM
         RGNVOL = 2*RGNVOL*WIDTH(I)
         CENTER(I) = A(I) + WIDTH(I)
      END DO
      BASEST = 0
      RGNERT = 0
*
*     Compute basic rule and error
*
 10   RGNVAL = 0
      RGNERR = 0
      RGNCMP = 0
      RGNCPT = 0
      DO I = 1,LENRUL
         FSYMSM = FULSUM(NDIM,CENTER,WIDTH,Z,G(1,I),FUNCTN,TID)
*     Basic Rule
         RGNVAL = RGNVAL + W(I,1)*FSYMSM
*     First comparison rule
         RGNERR = RGNERR + W(I,2)*FSYMSM
*     Second comparison rule
         RGNCMP = RGNCMP + W(I,3)*FSYMSM
*     Third Comparison rule
         RGNCPT = RGNCPT + W(I,4)*FSYMSM
      END DO
*
*     Error estimation
*
      RGNERR = SQRT(RGNCMP**2 + RGNERR**2)
      RGNCMP = SQRT(RGNCPT**2 + RGNCMP**2)
      IF ( 4*RGNERR .LT. RGNCMP ) RGNERR = RGNERR/2
      IF ( 2*RGNERR .GT. RGNCMP ) RGNERR = MAX( RGNERR, RGNCMP )
      RGNERT = RGNERT +  RGNVOL*RGNERR
      BASEST = BASEST +  RGNVOL*RGNVAL
*
*     When subregion has more than one piece, determine next piece and
*      loop back to apply basic rule.
*
      DO I = 1,NDIM
         CENTER(I) = CENTER(I) + 2*WIDTH(I)
         IF ( CENTER(I) .LT. B(I) ) GO TO 10
         CENTER(I) = A(I) + WIDTH(I)
      END DO
      END
      DOUBLE PRECISION FUNCTION FULSUM(S,CENTER,HWIDTH,X,G,F,TID)
*
****  To compute fully symmetric basic rule sum
*
      EXTERNAL F
      INTEGER S, IXCHNG, LXCHNG, I, L, TID
      DOUBLE PRECISION CENTER(S), HWIDTH(S), X(S), G(S), F
      DOUBLE PRECISION INTSUM, GL, GI
      FULSUM = 0
*
*     Compute centrally symmetric sum for permutation of G
*
 10   INTSUM = 0
      DO I = 1,S
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
      END DO
 20   INTSUM = INTSUM + F(S,X,TID)
      DO I = 1,S
         G(I) = -G(I)
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
         IF ( G(I) .LT. 0 ) GO TO 20
      END DO
      FULSUM = FULSUM + INTSUM
*
*     Find next distinct permuation of G and loop back for next sum
*
      DO I = 2,S
         IF ( G(I-1) .GT. G(I) ) THEN
            GI = G(I)
            IXCHNG = I - 1
            DO L = 1,(I-1)/2
               GL = G(L)
               G(L) = G(I-L)
               G(I-L) = GL
               IF (  GL  .LE. GI ) IXCHNG = IXCHNG - 1
               IF ( G(L) .GT. GI ) LXCHNG = L
            END DO
            IF ( G(IXCHNG) .LE. GI ) IXCHNG = LXCHNG
            G(I) = G(IXCHNG)
            G(IXCHNG) = GI
            GO TO 10
         ENDIF
      END DO
*
*     End loop for permutations of G and associated sums
*
*     Restore original order to G's
*
      DO I = 1,S/2
         GI = G(I)
         G(I) = G(S+1-I)
         G(S+1-I) = GI
      END DO
      END
      SUBROUTINE DIFFER(NDIM, A, B, WIDTH, Z, DIF, FUNCTN,
     &     DIVAXN, DIFCLS, TID)
*
*     Compute fourth differences and subdivision axes
*
      EXTERNAL FUNCTN
      INTEGER I, NDIM, DIVAXN, DIFCLS, TID
      DOUBLE PRECISION
     &     A(NDIM), B(NDIM), WIDTH(NDIM), Z(NDIM), DIF(NDIM), FUNCTN
      DOUBLE PRECISION FRTHDF, FUNCEN, WIDTHI
      DIFCLS = 0
      DIVAXN = MOD( DIVAXN, NDIM ) + 1
      IF ( NDIM .GT. 1 ) THEN
         DO I = 1,NDIM
            DIF(I) = 0
            Z(I) = A(I) + WIDTH(I)
         END DO
 10      FUNCEN = FUNCTN(NDIM, Z, TID)
         DO I = 1,NDIM
            WIDTHI = WIDTH(I)/5
            FRTHDF = 6*FUNCEN
            Z(I) = Z(I) - 4*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z,TID)
            Z(I) = Z(I) + 2*WIDTHI
            FRTHDF = FRTHDF - 4*FUNCTN(NDIM,Z,TID)
            Z(I) = Z(I) + 4*WIDTHI
            FRTHDF = FRTHDF - 4*FUNCTN(NDIM,Z,TID)
            Z(I) = Z(I) + 2*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z,TID)
*     Do not include differences below roundoff
            IF ( FUNCEN + FRTHDF/8 .NE. FUNCEN )
     &           DIF(I) = DIF(I) + ABS(FRTHDF)*WIDTH(I)
            Z(I) = Z(I) - 4*WIDTHI
         END DO
         DIFCLS = DIFCLS + 4*NDIM + 1
         DO I = 1,NDIM
            Z(I) = Z(I) + 2*WIDTH(I)
            IF ( Z(I) .LT. B(I) ) GO TO 10
            Z(I) = A(I) + WIDTH(I)
         END DO
         DO I = 1,NDIM
            IF ( DIF(DIVAXN) .LT. DIF(I) ) DIVAXN = I
         END DO
      ENDIF
      END
      SUBROUTINE TRESTR(POINTR, SBRGNS, PONTRS, RGNERS)
****BEGIN PROLOGUE TRESTR
****PURPOSE TRESTR maintains a heap for subregions.
****DESCRIPTION TRESTR maintains a heap for subregions.
*            The subregions are ordered according to the size of the
*            greatest error estimates of each subregion (RGNERS).
*
*   PARAMETERS
*
*     POINTR Integer.
*            The index for the subregion to be inserted in the heap.
*     SBRGNS Integer.
*            Number of subregions in the heap.
*     PONTRS Real array of dimension SBRGNS.
*            Used to store the indices for the greatest estimated errors
*            for each subregion.
*     RGNERS Real array of dimension SBRGNS.
*            Used to store the greatest estimated errors for each
*            subregion.
*
****ROUTINES CALLED NONE
****END PROLOGUE TRESTR
*
*   Global variables.
*
      INTEGER POINTR, SBRGNS
      DOUBLE PRECISION PONTRS(*), RGNERS(*)
*
*   Local variables.
*
*   RGNERR Intermediate storage for the greatest error of a subregion.
*   SUBRGN Position of child/parent subregion in the heap.
*   SUBTMP Position of parent/child subregion in the heap.
*
      INTEGER SUBRGN, SUBTMP
      DOUBLE PRECISION RGNERR
*
****FIRST PROCESSING STATEMENT TRESTR
*
      RGNERR = RGNERS(POINTR)
      IF ( POINTR .EQ. PONTRS(1)) THEN
*
*        Move the new subregion inserted at the top of the heap
*        to its correct position in the heap.
*
         SUBRGN = 1
 10      SUBTMP = 2*SUBRGN
         IF ( SUBTMP .LE. SBRGNS ) THEN
            IF ( SUBTMP .NE. SBRGNS ) THEN
*
*              Find maximum of left and right child.
*
          IF ( RGNERS(INT(PONTRS(SUBTMP))) .LT.
     $         RGNERS(INT(PONTRS(SUBTMP+1))) ) SUBTMP = SUBTMP + 1
            ENDIF
*
*           Compare maximum child with parent.
*           If parent is maximum, then done.
*
            IF ( RGNERR .LT. RGNERS(INT(PONTRS(SUBTMP))) ) THEN
*
*              Move the pointer at position subtmp up the heap.
*
               PONTRS(SUBRGN) = INT(PONTRS(SUBTMP))
               SUBRGN = SUBTMP
               GO TO 10
            ENDIF
         ENDIF
      ELSE
*
*        Insert new subregion in the heap.
*
         SUBRGN = SBRGNS
 20      SUBTMP = SUBRGN/2
         IF ( SUBTMP .GE. 1 ) THEN
*
*           Compare child with parent. If parent is maximum, then done.
*
            IF ( RGNERR .GT. RGNERS(INT(PONTRS(SUBTMP))) ) THEN
*
*              Move the pointer at position subtmp down the heap.
*
               PONTRS(SUBRGN) = PONTRS(SUBTMP)
               SUBRGN = SUBTMP
               GO TO 20
            ENDIF
         ENDIF
      ENDIF
      PONTRS(SUBRGN) = POINTR
*
****END TRESTR
*
      END
