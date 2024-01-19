C     GL0502       19 MAY 90                                         MRP
C?IBM/GLD/GBR/VAX/UNX
      SUBROUTINE UHFOPN
C??
C?CDC
C     PROGRAM UHFOPN
C??
C=SCF OPTIONS
C1INTRODUCTION
C     LINK 0502
C*
C     --------------
C     U OF T VERSION
C     MAY 1990
C     --------------
C*
C     SOLUTION OF THE POPLE-NESBET EQUATIONS FOR SPIN-UNRESTRICTED
C     (UHF) OPEN SHELL SYSTEMS.
C     J.CHEM.PHYS. 22, 571 (1954)
C*
C1OPTIONS
C     ******************************************************************
C     OPTIONS ... IOP( )
C     ******************************************************************
C     SEE SUBROUTINE CLOSED, LINK 0501.
C     ******************************************************************
C==
C*
C/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C##
      PARAMETER (NA= 36, NB=200)
C###
      PARAMETER (NBB=NB*(NB+1)/2)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      COMMON /A/ IOP(99)
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3)
      COMMON /A/ CDUM(4),ICDUM(401)
C
      COMMON/SCF1/ACURCY,DA(NBB),DB(NBB)
      COMMON/SCF2/F(NB,NB),FF(NB)
      COMMON/SCF3/VV(NB),V(NB,NB)
C
      COMMON/GEN/E1,E2,E3,EDUM(3),DCONV,S(2),FX,FY,FZ,FE,EMOL,ESOL,
     1 DPOL(4),VCM
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,IODUM1(5),IFORMT,IODUM2(206)
C
      DATA TEN/10.0D0/
C
 1000 FORMAT('0AT TERMINATION TOTAL ENERGY IS',F17.6,'  HARTREES')
 1010 FORMAT('1UHF OPEN SHELL SCF',15X,'NUCLEAR REPULSION ENERGY IS',
     $ F17.9,'  HARTREES'/'0CONVERGENCE ON DENSITY MATRIX REQUIRED TO ',
     $ 'EXIT IS',1PD12.4/'0 CYCLE   ELECTRONIC ENERGY',6X,
     $ 'TOTAL ENERGY    CONVERGENCE   EXTRAPOLATION'/)
 1020 FORMAT('0THIS CALCULATION WAS PERFORMED WITH AN APPLIED',
     $ ' EXTERNAL ELECTRIC FIELD')
 1030 FORMAT(/'0THIS CALCULATION WAS DONE ACCORDING TO SELF-CONSISTENT',
     1 ' REACTION FIELD THEORY. ENERGY VALUES IN HARTREES ARE'//
     2 10X,'E(Q)',7X,'=',6X,'EMOL',6X,'+',6X,'ESOL'/1X,3F17.6)
C*
      IF(IOP(10).EQ.0)RETURN
C     FOR ED/BSSE, SEE IF THIS LINK IS REQUIRED FOR THIS PASS.
      IF (IOP(4).EQ.1 .OR. IOP(4).EQ.5 .OR. IOP(4).GE.7) RETURN
C     FOR ED/BSSE, SET CORRECT 2-ELECTRON INTEGRAL UNIT NUMBER.
      IF (IOP(4) .EQ. 6) THEN
         ISAVE = ITWOEL
         ITWOEL = IFORMT
      END IF
C     CLEAR ENERGY ORIGIN FLAG.
      IOP(94) = 0
      NTT=(NBASIS*(NBASIS+1))/2
C     FORMV READS IN THE OVERLAP MATRIX FROM DISK ... FILE 8
C     COMPUTES TRANSFORMATION MATRIX AND WRITES IT ON DISK FILE 4
      CALL FORMV(NBASIS)
C     DENSITY MATRICES READ IN FROM DISK ... FILES 19 AND 21
      CALL TREAD(19,DA,NBB,1,NTT,1,0)
      CALL TREAD(21,DB,NBB,1,NTT,1,0)
      ACURCY=TEN**(-(5+IOP(13)))
      MAXCYC=30+20*(IOP(14)+1)
C     CALCULATE NUCLEAR REPULSION ENERGY
      E1=ZERO
      NM1=NATOMS-1
      IF(NM1.EQ.0) GO TO 81
      DO 80 I=1,NM1
      IANI=IAN(I)
      IF(IANI.LE.0)GO TO 80
      IP1=I+1
      DO 70 J=IP1,NATOMS
      IANJ=IAN(J)
      IF(IANJ.LE.0)GO TO 70
      ANIJ=DFLOAT(IANI*IANJ)
      RIJ=DSQRT((C(I,1)-C(J,1))**2+(C(I,2)-C(J,2))**2+
     $ (C(I,3)-C(J,3))**2)
      E1=E1+ANIJ/RIJ
   70 CONTINUE
   80 CONTINUE
C     INSERT CALCULATED ENERGIES IN COMMON GEN
C     E1 IS NUCLEAR REPULSION ENERGY.
C     E2 IS ELECTRONIC ENERGY.
C     E3 IS TOTAL ENERGY.
   81 IF(IOP(3).NE.1)GO TO 90
C     ADD THE FIELD CONTRIBUTION TO THE NUCLEAR ENERGY, FOR EXTERNAL
C     UNIFORM ELECTRIC FIELD ONLY.
      DO 85 I=1,NATOMS
      IANI=IAN(I)
      IF(IANI.LE.0)GO TO 85
      E1=E1-DFLOAT(IANI)*(FX*C(I,1)+FY*C(I,2)+FZ*C(I,3))
   85 CONTINUE
C
   90 IF(IOP(21).EQ.0)WRITE(IOUT,1010)E1,ACURCY
C     ENTER ITERATION ROUTINE CYCOPN - RETURN WHEN THE RMS DIFFERENCE
C     BETWEEN 2 SUCCESSIVE DENSITY MATRICES IS LESS THAN ACURCY
      CALL CYCOPN(MAXCYC)
C     FOR ED/BSSE, RESET 2-ELECTRON INTEGRAL UNIT NUMBER.
      IF (IOP(4) .EQ. 6) THEN
         ITWOEL = ISAVE
      END IF
C     CONVERGED IF IOP(1) IS 0.
      IF(IOP(1).EQ.0.AND.IOP(21).EQ.0)WRITE(IOUT,1000)E3
      DCONV=ACURCY
C     CALCULATE S, S(S+1).
      T=-TEN
      CALL SPIN (T, ONE)
      IF (IOP(21) .NE. 0) RETURN
      IF(IOP(3).EQ.0)RETURN
      IF(IOP(3).EQ.2)GO TO 100
      WRITE(IOUT,1020)
      RETURN
  100 WRITE(IOUT,1030)E3,EMOL,ESOL
      RETURN
      END
      SUBROUTINE CYCOPN(MAXCYC)
C*
C     --------------
C     U OF T VERSION
C     MAY 1990
C     --------------
C*
C     ITERATION ROUTINE
C     PARAMETERS AS FOLLOWS
C     ACURCY ... CONVERGENCE REQUIRED ON DENSITY MATRIX
C     RETURNS WITH ACTUAL CONVERGENCE ACHEIVED
C     MAXCYC ... MAXIMUM NUMBER OF ITERATIONS ALLOWED
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C##
      PARAMETER (NA= 36, NB=200)
C###
      PARAMETER (NBB=NB*(NB+1)/2, NBB1=NBB+1)
      PARAMETER (ZERO=0.0D0)
C
      COMMON /A/ IOP(99)
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3)
      COMMON /A/ DUMC(4),IDUMC(401)
C
      COMMON/SCF1/ACURCY,D(NB,NB),DD(NB)
      COMMON/SCF2/F(NB,NB),FF(NB)
      COMMON/SCF3/VV(NB),V(NB,NB)
C
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,IODUM(12),NFILE(100,2)
      COMMON/GEN/E1,E2,E3,EDUM(3),DCONV,SPIN(2),FIELD(3),FE,EMOL,ESOL,
     1 DPOL(4),VCM
C
      CHARACTER*8 LABEL
C
      DIMENSION DA(NBB),DB(NBB),FA(NBB),FB(NBB),
     1 VA(NBB),DEQ(NBB1),FEQ(NBB1)
C
      EQUIVALENCE (DA(1),D(1,1),DEQ(1)), (FA(1),F(1,1),FEQ(1)),
     1 (VA(1),V(1,1))
      EQUIVALENCE (DB(1),DEQ(NBB1)), (FB(1),FEQ(NBB1))
C
      DATA PT5/0.5D0/, TEN/10.0D0/
C
 1000 FORMAT(1X,I4,1X,2F20.9,1PD14.5,6X,A8)
 1010 FORMAT(///)
 1020 FORMAT('1COEFFICIENTS (ALPHA SPIN)')
 1030 FORMAT('1COEFFICIENTS (BETA SPIN)')
 1040 FORMAT('1DENSITY MATRIX (ALPHA SPIN)')
 1050 FORMAT('0*** SCF DID NOT CONVERGE AFTER',I4,' ITERATIONS ***')
 1060 FORMAT('1DENSITY MATRIX (BETA SPIN)')
 1070 FORMAT(' E =',F20.9,10X,'DELTA E =',F15.9,10X,'CONVERGENCE =',
     1 1PD13.5,10X,'RUN ABORTED')
 1080 FORMAT(' E =',F20.9,10X,'DELTA E =',F15.9,10X,'CONVERGENCE =',
     1 1PD13.5,10X,'RUN ALLOWED TO CONTINUE')
C*
C     PRESET CONVERGENCE ROUTINES
C     SETOPN IS AN ENTRY POINT IN CONOPN
      CALL SETOPN
      JCYCLE=1
      NTT=(NBASIS*(NBASIS+1))/2
      IEXT=1
      IPRINT=IOP(15)+IOP(16)
      IOP21=IOP(21)
C     PRESET 2-ELECTRON INTEGRAL INPUT ROUTINE.
      CALL E2SETO
C     JUMP INTO LOOP FOR FIRST TIME ... CALCULATE ENERGY.
      E2=ZERO
      GO TO 30
   10 CALL CONOPN(KEY,IEXT,LABEL)
C     PRINT ITERATION RESULTS.
      IF(IOP21.NE.0)GO TO 15
      IF(IPRINT.NE.0)WRITE(IOUT,1010)
      IF(IEXT.EQ.0)WRITE(IOUT,1000)JCYCLE,E2,E3,ACURCY,LABEL
      IF(IEXT.NE.0)WRITE(IOUT,1000)JCYCLE,E2,E3
   15 JCYCLE=JCYCLE+1
C     HAS CONVERGENCE BEEN MET ... IF SO ... EXIT
C     IF DENSITY MATRIX HAS JUST BEEN EXTRAPOLATED, CONTINUE FOR 2
C     MORE ITERATIONS TO ALLOW POSSIBLE CONVERGENCE
      IF(KEY)30,120,20
C     HAVE WE EXCEEDED THE MAXIMUM ALLOWED NUMBER OF CYCLES ... IF
C     SO ... EXIT
   20 IF(MAXCYC-JCYCLE)121,30,30
   30 IF(IEXT.EQ.0)EPREV=E2
C     FORM TWO ELECTRON PART OF THE FOCK MATRIX
      CALL FOFOPN (NTT, IOP, NBASIS)
      IF (IOP(1) .NE. 0) RETURN
C     THE CORE HAMILTONIAN IS READ FROM DISK ... FILE 13 ... INTO VA.
      CALL TREAD(13,VA,NBB,1,NTT,1,0)
C     ADD THE CORE HAMILTONIAN TO THE FOCK MATRICES AND CALCULATE
C     THE ENERGY.
      E2=ZERO
      K=0
      DO 90 I=1,NBASIS
      DO 80 J=1,I
      K=K+1
      FA(K)=FA(K)+VA(K)
      FB(K)=FB(K)+VA(K)
   80 E2=E2+(FA(K)+VA(K))*DA(K)+(FB(K)+VA(K))*DB(K)
   90 E2=E2-((FA(K)+VA(K))*DA(K)+(FB(K)+VA(K))*DB(K))*PT5
      E3=E1+E2
C
C     FORM ALPHA AND BETA DENSITY DENSITY MATRICES INDEPENDENTLY
C
C     FIRST ALPHA
C     WRITE FA ON DISK ... FILE 15 ... IF REQUESTED
      IF(IOP(19).NE.0)CALL TWRITE(15,FA,NBB,1,NTT,1,0)
C     PLACE FB ON DISK ... FILE 17 ... FOR A WHILE
      CALL TWRITE(17,FB,NBB,1,NTT,1,0)
C     EXPAND FA INTO A SQUARE ARRAY F ... F IS ACTUALLY EQUIVALENCED
C     TO FA
      CALL SQUARE(FA,F,NBASIS)
C
C     CHECK FOR THE PRESENCE OF A FIELD (EXTERNAL OR SCRFT).
C
      IF (IOP(3) .NE. 0) CALL SCRFT (15, 19, 1, F, D)
C     READ IN THE V MATRIX ... DISK FILE 4
      CALL TREAD(4,V,NB,NB,NBASIS,NBASIS,0)
C     FA  (IN F)  IS TRANSFORMED BY V ... TWO STEPS
      CALL MATPAC(V,F,D,NBASIS,1)
      CALL MATPAC(D,V,F,NBASIS,3)
C     THEN DIAGONALIZED
      CALL TRED12(NBASIS,NB,F,D,FF,DD,2)
C     EIGENVALUES WILL BE SAVED ... TO BE COMBINED WITH BETA EIGENVALUES
C     AND EVENTUALLY WRITTEN ON THE DISK ... FILE 5
C     THE EIGENVECTORS ARE TRANSFORMED BY THE V MATRIX
C     AND THE MOLECULAR ORBITAL COEFFICIENT MATRIX IS FORMED
      CALL MATPAC(V,D,F,NBASIS,2)
C     UPDATE DISK WITH M.O. COEFFICIENTS ... FILE 25
      CALL TWRITE(25,F,NB,NB,NBASIS,NBASIS,0)
      IF(IOP(15).NE.1)GO TO 98
      WRITE(IOUT,1020)
      CALL GBSOUT(F,FF,NB,NB,NBASIS,1)
C     THE DENSITY MATRIX IS FORMED FROM THE COEFFICIENTS
   98 DO 100 I=1,NBASIS
      DO 100 J=1,I
      SUM=ZERO
      DO 101 K=1,NAE
  101 SUM=SUM+F(I,K)*F(J,K)
      D(I,J)=SUM
  100 D(J,I)=SUM
C     WRITE OUT D  (DA)  ON DISK ... FILE 19
      CALL TWRITE(19,D,NB,NB,NBASIS,NBASIS,1)
      IF(IOP(16).NE.1)GO TO 104
      WRITE(IOUT,1040)
      CALL GBSOUT(D,FF,NB,NB,NBASIS,0)
C
C     BETA SPIN.
C
C     BRING IN FB ... DISK FILE 17
C     ENTER INTO SLOT MARKED FA ... JUST FOR CONVENIENCE
  104 CALL TREAD(17,FA,NBB,1,NTT,1,0)
C     NOW REPEAT THE WHOLE DAMN PROCESS TO FORM ... THIS TIME ... DB
      CALL SQUARE(FA,F,NBASIS)
C
C     CHECK FOR THE PRESENCE OF A FIELD (EXTERNAL OR SCRFT).
C
      IF (IOP(3) .NE. 0) THEN
	 CALL SCRFT (17, 21, 0, F, D)
C        READ IN THE V MATRIX ... DISK FILE 4
         CALL TREAD(4,V,NB,NB,NBASIS,NBASIS,0)
      END IF
C
C     FOR ED/BSSE, SAVE THE FIRST CYCLE ENERGY IN /GEN/ IF NECESSARY.
C
      IF (JCYCLE .EQ. 1) THEN
         IF (IOP(4).EQ.4 .OR. IOP(4).EQ.6) EDUM(1) = E3
      END IF
C
      CALL MATPAC(V,F,D,NBASIS,1)
      CALL MATPAC(D,V,F,NBASIS,3)
C     THIS TIME PUT EIGENVALUES IN VV SO AS NOT TO DESTROY ALPHA
C     EIGENVALUES
      CALL TRED12(NBASIS,NB,F,D,VV,DD,2)
C     COMBINE ALPHA AND BETA EIGENVALUES IN FA ... AND SHIP OUT ON
C     DISK ... FILE 5
      DO 105 I=1,NBASIS
      FA(I)=FF(I)
  105 FA(I+NBASIS)=VV(I)
      CALL TWRITE(5,FA,NBB,1,NBASIS+NBASIS,1,0)
      CALL MATPAC(V,D,F,NBASIS,2)
C     BETA M.O. COEFFICIENTS TO DISK ... FILE 27
      CALL TWRITE(27,F,NB,NB,NBASIS,NBASIS,0)
      IF(IOP(15).NE.1)GO TO 108
      WRITE(IOUT,1030)
      CALL GBSOUT(F,VV,NB,NB,NBASIS,1)
  108 DO 111 I=1,NBASIS
      DO 111 J=1,I
      SUM=ZERO
      IF(NBE.EQ.0)GO TO 110
      DO 109 K=1,NBE
  109 SUM=SUM+F(I,K)*F(J,K)
  110 D(I,J)=SUM
  111 D(J,I)=SUM
C     BETA DENSITY MATRIX TO DISK ... FILE 21
      CALL TWRITE(21,D,NB,NB,NBASIS,NBASIS,1)
      IF(IOP(16).NE.1)GO TO 118
      WRITE(IOUT,1060)
      CALL GBSOUT(D,FF,NB,NB,NBASIS,0)
C     READ BOTH DENSITY MATRICES BACK IN ... INTO LINEAR FORM
  118 CALL TREAD(19,DA,NBB,1,NTT,1,0)
      CALL TREAD(21,DB,NBB,1,NTT,1,0)
      GO TO 10
C     MAXIMUM ITERATION COUNT EXCEEDED WITHOUT CONVERGENCE
  121 WRITE(IOUT,1050)MAXCYC
      IOP(1)=-1
      IF(IOP(22).EQ.0)RETURN
C     ALLOW RUN TO CONTINUE IF ENERGY CHANGE IS SMALL ENOUGH.
      EPREV=E2-EPREV
      IF(DABS(EPREV).LE.TEN**(-8+IOP(22)))GO TO 122
      WRITE(IOUT,1070)E3,EPREV,ACURCY
      RETURN
C     ALLOW RUN TO CONTINUE.
  122 WRITE(IOUT,1080)E3,EPREV,ACURCY
      IOP(1)=0
  120 RETURN
      END
      SUBROUTINE FOFOPN (NTT, IOP, NBASIS)
C*
C     --------------
C     U OF T VERSION
C     DECEMBER 1988
C     --------------
C*
C     SUBROUTINE TO FORM THE TWO ELECTRON PART OF THE ALPHA AND BETA
C     FOCK MATRICES.
C     GIVEN THE ALPHA AND BETA DENSITY MATRICES AND ANY OR ALL TWO
C     ELECTRON INTEGRALS.
C     BASED ON GAUSSIAN 76.
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C##
      PARAMETER (NA= 36, NB=200)
C###
      PARAMETER (NBB=NB*(NB+1)/2, NB1=NB*(NB+1),
     1 NBP1=NB+1)
      PARAMETER (ZERO=0.0D0)
C
      INTEGER P,Q,R,SINDX
      LOGICAL USEOFF
C
      COMMON/SCF1/ACURCY,DA(NBB),DB(NBB)
      COMMON/SCF2/FA(NBB),FB(NBB)
      COMMON/SCF3/FILL(NB1)
C
      COMMON/BSINFO/NVO(NA),IAOFF(NA)
      COMMON/PACKED/I,J,K,L,JA
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,IODUM(12),NFILE(100,2)
C
      DIMENSION IJ(NBP1), IX(3072,2), IOP(99), IV(2)
C
      EQUIVALENCE (IX(1,1),FILL(1)), (VALINT,IV(1))
      EQUIVALENCE (P,I),(Q,J),(R,MINDX,K),(SINDX,L)
C
      SAVE IJ
C
C?IBM/GLD/GBR/VAX/UNX
      DATA NWORDS/3/, MAXBUF/3070/
C??
C?CDC NO OFFSETS ARE USED.
C     DATA NWORDS/2/, MAXBUF/3071/, USEOFF/.FALSE./
C??
C?IBM/GLD/GBR/VAX/UNX
      USEOFF = NBASIS .GT. 127
C??
C     SETUP FOR DOUBLE BUFFERING (ASYNCHRONOUS I/O).
C     LBUFF IS THE CURRENT BUFFER NUMBER (BEING UNPACKED).
C     NBUFF IS THE NEXT BUFFER NUMBER (I/O IN PROGRESS).
C
      NBUFF = 1
      CALL NWIORD (ITWOEL, IX(1,NBUFF), 3072, .FALSE., 'FOFOPN',
     1 IOP(1), *9000)
C
C     ZERO THE ARRAYS FA AND FB, WHERE THE F MATRICES ARE TO BE PLACED.
C
      DO 10 M=1,NTT
      FA(M)=ZERO
   10 FB(M)=ZERO
C
C     THE FILLED BUFFER LOOKS AS FOLLOWS FOR 32 BIT MACHINES:
C     LABEL-INTEGRAL-INTEGRAL-LABEL-INTEGRAL-INTEGRAL-LABEL-----------
C     THE LABEL BEING A PACKED WORD CONTAINING THE NUMBERS OF THE FOUR
C     ATOMIC ORBITALS ASSOCIATED WITH THE INTEGRAL - THE SIGN
C     SERVES TO CLASSIFY THE INTEGRAL TYPE (1 OR 2) - SEE SHLOUT.
C
C     WAIT FOR PREVIOUS BUFFER TO BE READ, THEN START READING THE NEXT
C     ONE.
C
   20 LBUFF = NBUFF
      CALL NWIOWT (ITWOEL, 3072, .FALSE., 'FOFOPN', IOP(1), *9000)
      NBUFF = 3 - NBUFF
      CALL NWIORD (ITWOEL, IX(1,NBUFF), 3072, .FALSE., 'FOFOPN',
     1 IOP(1), *9000)
C
C     THIS IS THE LOOP OVER THE BUFFER OF INTEGRALS AND LABELS
C     TREATING EACH LABEL-INTEGRAL COMBINATION SEPARATELY.
C     THE ROUTINE USES A LABEL OF 0 TO SIGNIFY THAT A NEW SET OF ATOM
C     OFFSETS IS TO BE APPLIED.
C     FOR 32 BIT MACHINES, EACH LABEL/INTEGRAL OCCUPIES 3 WORDS.
C     FOR 60 BIT MACHINES, EACH LABEL/INTEGRAL OCCUPIES 2 WORDS.
C
      DO 300 M=1,MAXBUF,NWORDS
      JA=IX(M,LBUFF)
      IV(1)=IX(M+1,LBUFF)
C?IBM/GLD/GBR/VAX/UNX
      IV(2)=IX(M+2,LBUFF)
C??
C     UNPACK THE LABELS PREVIOUSLY STORED AS ONE WORD.
C     THE FOLLOWING CODING FINDS ALL THE UNIQUE PERMUTATIONS OF THE
C     INTEGRAL LABELS AND CORRESPONDINGLY ADDS THE VARIOUS
C     CONTRIBUTIONS OF THE INTEGRAL TO THE F MATRICES.
C     CAUTION ** ONLY A SINGLE PERMUTATION OF THE INTEGRAL LABELS MAY
C     APPEAR IN THE INPUT - OTHERWISE ALL THE PERMUTATIONS ON THAT
C     INTEGRAL WILL BE ADDED TO THE F MATRIX MORE THAN ONCE.
C     I, J, K AND L ARE THE VARIOUS LABELS.
C     SEE PROGRAM SHLOUT (LINK 0306) FOR MORE DETAILS.
C
      IF(JA)30,200,40
C
C     TYPE 2 INTEGRALS ARE CLASSIFIED BY SINDX.
C
   30 JA=-JA
C?IBM/GLD/GBR/VAX/CDC
C     CALL UNPACK
C??
C?UNX
      CALL UNPACK (JA, I, J, K, L)
C??
      IF (USEOFF) GO TO 35
C     DO NOT USE OFFSETS.
      GO TO (101,111,121,131,141,151,161,171),SINDX
C     USE OFFSETS.
   35 GO TO (100,110,120,130,140,150,160,170),SINDX
C
C     ( AB I CD ) - TYPE 1 INTEGRALS.
C
C?IBM/GLD/GBR/VAX/CDC
C  40 CALL UNPACK
C??
C?UNX
   40 CALL UNPACK (JA, I, J, K, L)
C??
      IF (USEOFF) THEN
         I=I+IOFF
         J=J+JOFF
         K=K+KOFF
         L=L+LOFF
      END IF
      JS=IJ(I)+J
      JT=IJ(K)+L
      TERINT=VALINT+VALINT
      SECINT=(DA(JT)+DB(JT))*TERINT
      FA(JS)=FA(JS)+SECINT
      FB(JS)=FB(JS)+SECINT
      SECINT=(DA(JS)+DB(JS))*TERINT
      FA(JT)=FA(JT)+SECINT
      FB(JT)=FB(JT)+SECINT
      JS=IJ(I)+K
      IF(J.LT.L)GO TO 50
      JT=IJ(J)+L
      GO TO 60
   50 JT=IJ(L)+J
   60 FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      JS=IJ(I)+L
      IF(J.LT.K)GO TO 80
      JT=IJ(J)+K
      GO TO 90
   80 JT=IJ(K)+J
   90 FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      GO TO 300
C
C     ( AB I AC )
C
  100 P=P+IOFF
      Q=Q+JOFF
      R=R+LOFF
  101 JS=IJ(P)+Q
      JT=IJ(P)+R
      TERINT=VALINT+VALINT
      SECINT=(DA(JT)+DB(JT))*TERINT
      FA(JS)=FA(JS)+SECINT
      FB(JS)=FB(JS)+SECINT
      SECINT=(DA(JS)+DB(JS))*TERINT
      FA(JT)=FA(JT)+SECINT
      FB(JT)=FB(JT)+SECINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      JS=IJ(P+1)
      JT=IJ(Q)+R
      FA(JS)=FA(JS)-DA(JT)*TERINT
      FB(JS)=FB(JS)-DB(JT)*TERINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      GO TO 300
C
C     ( AB I BC )
C
  110 P=P+IOFF
      Q=Q+JOFF
      R=R+LOFF
  111 JS=IJ(P)+Q
      JT=IJ(Q)+R
      TERINT=VALINT+VALINT
      SECINT=(DA(JT)+DB(JT))*TERINT
      FA(JS)=FA(JS)+SECINT
      FB(JS)=FB(JS)+SECINT
      SECINT=(DA(JS)+DB(JS))*TERINT
      FA(JT)=FA(JT)+SECINT
      FB(JT)=FB(JT)+SECINT
      FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      JS=IJ(P)+R
      JT=IJ(Q+1)
      FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*TERINT
      FB(JT)=FB(JT)-DB(JS)*TERINT
      GO TO 300
C
C     ( AC I BC )
C
  120 P=P+IOFF
      Q=Q+JOFF
      R=R+KOFF
  121 JS=IJ(P)+Q
      JT=IJ(R)+Q
      TERINT=VALINT+VALINT
      SECINT=(DA(JT)+DB(JT))*TERINT
      FA(JS)=FA(JS)+SECINT
      FB(JS)=FB(JS)+SECINT
      SECINT=(DA(JS)+DB(JS))*TERINT
      FA(JT)=FA(JT)+SECINT
      FB(JT)=FB(JT)+SECINT
      FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      JS=IJ(P)+R
      JT=IJ(Q+1)
      FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*TERINT
      FB(JT)=FB(JT)-DB(JS)*TERINT
      GO TO 300
C
C     ( AB I AB )
C
  130 P=P+IOFF
      Q=Q+JOFF
  131 JS=IJ(P)+Q
      SECINT=DA(JS)+DB(JS)
      FA(JS)=FA(JS)+(SECINT+DB(JS))*VALINT
      FB(JS)=FB(JS)+(SECINT+DA(JS))*VALINT
      JS=IJ(P+1)
      JT=IJ(Q+1)
      FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      GO TO 300
C
C     ( AA I BC )
C
  140 P=P+IOFF
      Q=Q+KOFF
      R=R+LOFF
  141 JS=IJ(P+1)
      JT=IJ(Q)+R
      SECINT=(DA(JT)+DB(JT))*VALINT
      SECINT=SECINT+SECINT
      FA(JS)=FA(JS)+SECINT
      FB(JS)=FB(JS)+SECINT
      SECINT=(DA(JS)+DB(JS))*VALINT
      FA(JT)=FA(JT)+SECINT
      FB(JT)=FB(JT)+SECINT
      JS=IJ(P)+Q
      JT=IJ(P)+R
      FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      GO TO 300
C
C     ( AB I CC )
C
  150 P=P+IOFF
      Q=Q+JOFF
      R=R+KOFF
  151 JS=IJ(P)+Q
      JT=IJ(R+1)
      SECINT=(DA(JT)+DB(JT))*VALINT
      FA(JS)=FA(JS)+SECINT
      FB(JS)=FB(JS)+SECINT
      SECINT=(DA(JS)+DB(JS))*VALINT
      SECINT=SECINT+SECINT
      FA(JT)=FA(JT)+SECINT
      FB(JT)=FB(JT)+SECINT
      JS=IJ(P)+R
      IF(Q.LT.R)GO TO 155
      JT=IJ(Q)+R
      GO TO 156
  155 JT=IJ(R)+Q
  156 FA(JS)=FA(JS)-DA(JT)*VALINT
      FB(JS)=FB(JS)-DB(JT)*VALINT
      FA(JT)=FA(JT)-DA(JS)*VALINT
      FB(JT)=FB(JT)-DB(JS)*VALINT
      GO TO 300
C
C     SINDX=6 ... 3 SUBCASES.
C
  160 P=P+IOFF
      Q=Q+LOFF
  161 IF(MINDX-2)165,166,167
C
C     ( AB I BB )
C
  165 JS=IJ(P)+Q
      JT=IJ(Q+1)
      FA(JS)=FA(JS)+DB(JT)*VALINT
      FB(JS)=FB(JS)+DA(JT)*VALINT
      SECINT=VALINT+VALINT
      FA(JT)=FA(JT)+DB(JS)*SECINT
      FB(JT)=FB(JT)+DA(JS)*SECINT
      GO TO 300
C
C     ( AA I BB )
C
  166 JS=IJ(P+1)
      JT=IJ(Q+1)
      SECINT=(DA(JT)+DB(JT))*VALINT
      FA(JS)=FA(JS)+SECINT
      FB(JS)=FB(JS)+SECINT
      SECINT=(DA(JS)+DB(JS))*VALINT
      FA(JT)=FA(JT)+SECINT
      FB(JT)=FB(JT)+SECINT
      JS=IJ(P)+Q
      FA(JS)=FA(JS)-DA(JS)*VALINT
      FB(JS)=FB(JS)-DB(JS)*VALINT
      GO TO 300
C
C     ( AA I AB )
C
  167 JS=IJ(P+1)
      JT=IJ(P)+Q
      SECINT=VALINT+VALINT
      FA(JS)=FA(JS)+DB(JT)*SECINT
      FB(JS)=FB(JS)+DA(JT)*SECINT
      FA(JT)=FA(JT)+DB(JS)*VALINT
      FB(JT)=FB(JT)+DA(JS)*VALINT
      GO TO 300
C
C     ( AA I AA )
C
  170 P=P+IOFF
  171 JS=IJ(P+1)
      FA(JS)=FA(JS)+DB(JS)*VALINT
      FB(JS)=FB(JS)+DA(JS)*VALINT
      GO TO 300
C
C     GET NEW ATOM OFFSET VALUES.
C
  200 JA = IV(1)
      IF (JA .EQ. 0) GO TO 310
C?IBM/GLD/GBR/VAX/CDC
C     CALL UNPACK
C??
C?UNX
      CALL UNPACK (JA, I, J, K, L)
C??
      IOFF=IAOFF(I)
      JOFF=IAOFF(J)
      KOFF=IAOFF(K)
      LOFF=IAOFF(L)
  300 CONTINUE
      GO TO 20
C
C     WAIT FOR FINAL BUFFER TO BE READ - IT SHOULD BE AT END OF FILE.
C
  310 CALL NWIOWT (ITWOEL, 3072, .TRUE., 'FOFOPN', IOP(1), *9000)
      REWIND ITWOEL
      RETURN
C
C     ERROR EXIT.
C
 9000 RETURN
C
      ENTRY E2SETO
C*
C     FORM INDEX ARRAY IJ.
C     THIS IS USED IN CONJUNCTION WITH THE CONVERSION OF DOUBLE
C     TO SINGLE SUBSCRIPTS.
C
      DO 330 M=1,NBP1
  330 IJ(M)=(M*(M-1))/2
      RETURN
      END
      SUBROUTINE CONOPN(KEY,IEXT,LABEL)
C*
C     --------------
C     U OF T VERSION
C     DECEMBER 1988
C     --------------
C*
C     UHF OPEN SHELL SCF CONVERGENCE ROUTINES
C     USES DISK FILES 31 TO 33 FOR SCRATCH
C     EACH MAXIMUM OF #NB1 DOUBLE WORDS LONG
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C##
      PARAMETER (NA= 36, NB=200)
C###
      PARAMETER (NBB=NB*(NB+1)/2, NBB1=NBB+1, NB1=NB*(NB+1))
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
      COMMON /A/ IOP(99),IC1(NA),IC2(6),N,CD1(NA,3),CD2(4),IC3(401)
      COMMON/SCF1/ACURCY,A3(NB1)
      COMMON/SCF2/A2(NB1)
      COMMON/SCF3/A1(NB1)
C
      CHARACTER*8 LABEL
C
      SAVE
C
      DATA PT99/0.99D0/, PT995/0.995D0/, ONEPT9/1.9D0/
      DATA FOUR/4.0D0/, PT06/0.06D0/
C
      KEY=-1
      IEXT=1
      LABEL = ' '
      IOP17=IOP(17)+1
      ICOUNT=ICOUNT+1
C     SWITCH SCRATCH FILE DEFINITIONS.
      LOC1=LOC2
      LOC2=65-LOC1
C     COMPRESS ALPHA AND BETA DENSITY MATRICES INTO ONE ARRAY
      DO 60 I=1,NTT
   60 A3(NTT+I)=A3(NBB+I)
C     SKIP TO END IN FIRST CYCLE OR IMMEDIATELY AFTER EXTRAPOLATION
      IF(ICOUNT.EQ.1)GO TO 260
C     P(N) NOW IN A3.  READ P(N-1) INTO A1
      CALL TREAD(31,A1,NB1,1,NTT2,1,0)
C     FORM P(N)-P(N-1) IN A2
      DO 80 I=1,NTT2
   80 A2(I)=A3(I)-A1(I)
      IEXT=0
      KEY=1
C     FIND LENGTH DP1
      SP11=TRAOPN(A2,A2,N)
      DP1=DSQRT(SP11/TWO)
C     TEST FOR CONVERGENCE BY FINDING IF ROOT MEAN SQUARE DP IS LESS
C     THAN ACURCY - GET SAME CONVERGENCE AS CLOSED SHELL.
      ACURCY=(DP1+DP1)/DFLOAT(N)
      IF(ACURCY.GE.CRIT)GO TO 90
C     CONVERGENCE CRITERION MET ... EXIT
C     WATCH OUT FOR USE OF EXTRAPOLATED DENSITY MATRIX IN DEWAR METHOD.
      IF(IOP17.NE.2.OR.ICOUNT.GT.2)KEY=0
      RETURN
C     WHAT TYPE OF EXTRAPOLATION?
   90 GO TO (100,300,100,260),IOP17
C     POPLE'S GAUSSIAN 70 EXTRAPOLATION.
  100 IF(ICOUNT.LE.2)GO TO 240
      IF(ICOUNT.EQ.3)GO TO 110
      CALL TREAD(LOC1,A1,NB1,1,NTT2,1,0)
      SP23=SP12
      SP33=SP22
      SP13=TRAOPN(A1,A2,N)
C     FIND LENGTH DP3
      DP3=DSQRT(SP33/TWO)
C     READ P(N-1)-P(N-2) INTO A1
  110 CALL TREAD(LOC2,A1,NB1,1,NTT2,1,0)
      SP12=TRAOPN(A1,A2,N)
      SP22=TRAOPN(A1,A1,N)
C     FIND LENGTH DP2
      DP2=DSQRT(SP22/TWO)
      IF(ICOUNT.EQ.3)GO TO 240
C     FIND COSINE OF ANGLE BETWEEN SUCCESSIVE DISPLACEMENTS
      COSPHI=SP12/(TWO*DP1*DP2)
C     FIND COSINE OF ANGLE BETWEEN DP3 AND PLANE OF DP1 AND DP2
      X=(SP13*SP22-SP12*SP23)/(SP11*SP22-SP12*SP12)
      Y=(SP23*SP11-SP12*SP13)/(SP11*SP22-SP12*SP12)
      COSPSI=DSQRT((X*X*SP11+Y*Y*SP22+TWO*X*Y*SP12)/TWO)/DP3
C     DO NOT EXTRAPOLATE UNLESS 4 CONSECUTIVE POINTS ARE NEARLY COPLANAR
      IF(COSPSI.LE.PT99)GO TO 240
C     EXPRESS VECTOR DP1 AS X*DP3(PROJECTED)+Y*DP2
      Y=-Y/X
      X=ONE/X
C     TEST IF 2*2 MATRIX HAS REAL EIGENVALUES BETWEEN -.95 AND +.95
      XY=Y*Y+FOUR*X
      IF(XY.LT.ZERO)GO TO 190
      XY=DABS(Y)+DSQRT(XY)
      IF(XY.LE.ONEPT9)GO TO 220
  190 IF(IOP(18).EQ.1)GO TO 240
C     IF 4-POINT EXTRAPOLATION IS NOT POSSIBLE TRY 3-POINT
      IF(DABS(COSPHI).LE.PT995)GO TO 240
      X=DP1/(DP2*COSPHI-DP1)
      DO 210 I=1,NTT2
  210 A3(I)=A3(I)+X*A2(I)
      LABEL = '3-POINT'
      ICOUNT=0
      KEY=-1
      GO TO 270
  220 XXX=X/(ONE-X-Y)
      YYY=(X+Y)/(ONE-X-Y)
      DO 230 I=1,NTT2
  230 A3(I)=A3(I)+XXX*A1(I)+YYY*A2(I)
      LABEL = '4-POINT'
      ICOUNT=0
      KEY=-1
      GO TO 270
  240 CALL TWRITE(LOC1,A2,NB1,1,NTT2,1,0)
  260 CALL TWRITE(31,A3,NB1,1,NTT2,1,0)
C     SEPARATE ALPHA AND BETA DENSITY MATRICES
  270 DO 280 I=1,NTT
      I1=NBB1+NTT-I
      I2=NTT2+1-I
  280 A3(I1)=A3(I2)
      RETURN
C     DEWAR EXTRAPOLATION: MJS DEWAR AND PK WEINER, COMP & CHEM, 2,
C     31 (1978). CALCULATE A NEW VALUE FOR LAMBDA AND SAVE PREVIOUS
C     VALUE IN Y2.
  300 IDEWIT=IDEWIT+1
C     IDEWIT CONTROLS THE EXTRAPOLATION PROCEDURE:
C     1 ... JUST EXTRAPOLATED THE DENSITY MATRIX,
C     2 ... CALCULATE LAMBDA2,
C     3 ... CALCULATE LAMBDA3 AND EXTRAPOLATE.
      IF(IDEWIT-2)240,305,303
  303 Y2=Y
  305 CALL TREAD(LOC2,A1,NB1,1,NTT2,1,0)
      SP22=TRAOPN(A1,A1,N)
      SP12=TRAOPN(A1,A2,N)
      Y=SP12/SP22
C     CHECK FOR DIVERGENCE.
      IF(DABS(Y).GE.ONE)GO TO 330
C     EXTRAPOLATE IF 2 SUCCESSIVE LAMBDA'S AGREE WITHIN DCRIT.
      IF(IDEWIT.EQ.2)GO TO 240
      IF(DABS(Y-Y2).GT.DCRIT)GO TO 240
      LABEL = ' DEWAR'
  310 X=Y/(ONE-Y)
      DO 320 I=1,NTT2
  320 A3(I)=A3(I)+X*A2(I)
      ICOUNT=1
      KEY=-1
      IDEWIT=0
C     SAVE THIS EXTRAPOLATED DENSITY MATRIX EVEN THOUGH NOT IDEMPOTENT.
      GO TO 260
C     SPECIAL VALUE OF LAMBDA FOR DIVERGENT CASES.
  330 Y=SP12/SP11
      LABEL = ' DEWAR'''
      GO TO 310
C*
C     ENTRY TO PRESET CONVERGENCE ROUTINES
C*
      ENTRY SETOPN
      ICOUNT=0
      NTT=(N*(N+1))/2
      NTT2=NTT+NTT
      LOC2=32
C     PRESET REQUIRED CONVERGENCE ON THE DENSITY MATRIX
      CRIT=ACURCY
      IDEWIT=0
      DCRIT=PT06*DFLOAT(IOP(20)+1)
      IF(IOP(20).GE.7)DCRIT=TWO
      RETURN
      END
      SUBROUTINE SPIN (AOPEN, C0)
C*
C     --------------
C     U OF T VERSION
C     JANUARY 1987
C     --------------
C*
C     COMPUTES S AND S(S+1) BY TRACE(DB*S*DA*S).
C
C     REFERENCE: F.L. PILAR, "ELEMENTARY QUANTUM CHEMISTRY",
C     MCGRAW-HILL, TORONTO, 1968, PP 287-293.
C
C     AOPEN GIVES THE TOTAL NUMBER OF ALPHA-BETA OPEN SHELL SINGLET
C     PAIRS IN THE OPEN SHELLS. A NEGATIVE VALUE INDICATES THAT
C     THIS ROUTINE IS TO CALCULATE S(S+1), THEN SUBTRACT 1.0 IF THE
C     VALUE DEVIATES GREATLY FROM ZERO - THIS WILL WORK FOR UHF
C     CALCULATIONS ON CLOSED SHELL SINGLETS (WHICH SHOULD GIVE A VALUE
C     OF 0.0 FOR S(S+1)), AND FOR OPEN SHELL SINGLETS WITH ONLY 2
C     OPEN SHELL ORBITALS. AOPEN IS IGNORED EXCEPT FOR SINGLET STATES,
C     IF AOPEN IS LESS THAN ZERO.
C
C     C0 GIVES THE COEFFICIENT OF THE REFERENCE CONFIGURATION FOR
C     THIS WAVEFUNCTION, AND WOULD BE 1.0D0 FOR UHF OR RHF OPEN
C     SHELL SCFS AS THIS VALUE IS ONLY NEEDED FOR MC-SCF OR GVB
C     WAVEFUNCTIONS.
C     THE TOTAL SPIN IS USUALLY INCORRECT FOR MULTI-CONFIGURATION
C     SCF OR GVB SCF CALCULATIONS - THE EXPRESSION TO CALCULATE THE
C     SPIN FOR SUCH WAVEFUNCTIONS IS UNDER INVESTIGATION, BUT FOR NOW
C     THE SPIN IS CALCULATED USING THE STANDARD FORMULA, THEN
C     REPLACED BY THE EXACT VALUE DERIVED FROM THE MULTIPLICITY.
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C##
      PARAMETER (NA= 36, NB=200)
C###
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
      COMMON /A/ IOP(99)
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3)
      COMMON /A/ CDUM(4),ICDUM(401)
C
      COMMON/SCF1/ACURCY,A(NB,NB),AA(NB)
      COMMON/SCF2/B(NB,NB),BB(NB)
      COMMON/SCF3/S(NB,NB),SS(NB)
C
      COMMON/GEN/E(6),DCONV,XS,SP,FIELD(3),FE,EMOL,ESOL,DPOL(4),VCM
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,IODUM(12),NFILE(100,2)
C
      LOGICAL MULTIC
C
      DATA PT25/0.25D0/, PT5/0.5D0/, FOUR/4.0D0/,
     1 PT01/0.01D0/, PT05/0.05D0/
C
 1000 FORMAT(/'0S =',F6.3,' (',F6.3,' EXPECTED),   S(S+1) =',F8.4,
     1 ' (',F8.4,' EXPECTED)')
 1010 FORMAT('0*** WARNING: SPIN ',
     1 'CONTAMINATION - GEOMETRIES/ENERGIES MAY BE AFFECTED ***')
 1020 FORMAT('0*** WARNING: SUBSTANTIAL ',
     1 'SPIN CONTAMINATION - GEOMETRIES/ENERGIES ARE SIGNIFICANTLY ',
     2 'AFFECTED ***')
 1030 FORMAT ('0*** THE S AND S(S+1) VALUES ABOVE HAVE BEEN REPLACED ',
     1 'BY THE EXACT VALUES DETERMINED FROM THE MULTIPLICITY'/
     2 ' AS THE FORMULA TO CALCULATE THE SPIN FOR MC-SCF/GVB ',
     3 'WAVEFUNCTIONS IS NOT YET AVAILABLE. ***')
C
C     READ OVERLAP MATRIX FROM DISK ... INTO S
      CALL TREAD(8,S,NB,NB,NBASIS,NBASIS,1)
C     AND ALPHA DENSITY MATRIX ... INTO A
      CALL TREAD(19,A,NB,NB,NBASIS,NBASIS,1)
C     FORM MATRIX PRODUCT S*A*S.
      CALL MATPAC(S,A,B,NBASIS,1)
      CALL MATPAC(B,S,A,NBASIS,1)
C     READ IN BETA DENSITY MATRIX ... INTO S
      CALL TREAD(21,S,NB,NB,NBASIS,NBASIS,1)
C     FORM TRACE OF THE PRODUCT  BETA*(S*ALPHA*S).
      TRACE=ZERO
      DO 60 I=1,NBASIS
      DO 50 J=1,I
      XY=A(I,J)*S(J,I)
   50 TRACE=TRACE+XY+XY
   60 TRACE=TRACE-XY
      APB=NAE+NBE
      AMB=NAE-NBE
      SP=PT5*APB+PT25*AMB**2-TRACE
C
C     WATCH FOR OPEN SHELL SINGLETS.
C
      IF(AOPEN.LT.ZERO)GO TO 70
C     REDUCE SP BY 1.0 FOR EACH BETA SPIN OPEN SHELL ELECTRON IN
C     OPEN SINGLET PAIRS.
      SP=SP-AOPEN
      GO TO 100
C     UHF CASE - IF SP IS CLOSE TO ZERO, ASSUME CLOSED SHELL; ELSE
C     SUBTRACT 1.0 FOR 1 OPEN SHELL SINGLET PAIR OF ELECTRONS.
   70 IF(AMB.NE.ZERO)GO TO 100
      SP=DABS(SP)
      IF(SP.GT.PT25)SP=SP-ONE
C
  100 SP=DABS(SP)
      MULTIC = .FALSE.
      IF (DABS(C0) .EQ. ONE) GO TO 200
C
C     APPLY EMPIRICAL CORRECTION FOR MC-SCF AND GVB WAVEFUNCTIONS:
C     THE TOTAL SPIN SEEMS TO BE TOO LARGE BY A FACTOR OF
C     2.0*(1.0-C0**2), SO REDUCE S**2 BY THIS AMOUNT.
C
      MULTIC = .TRUE.
      SP = DABS(SP-TWO*(ONE-C0*C0))
C
  200 XS=DABS(-PT5+PT5*DSQRT(ONE+FOUR*SP))
      XSE=PT5*DFLOAT(MULTIP-1)
      SPE=XSE*(XSE+ONE)
      IF (IOP(1).EQ.0 .AND. IOP(21).NE.0) GO TO 210
      WRITE (IOUT,1000) XS,XSE,SP,SPE
C
  210 IF (MULTIC) GO TO 400
      IF (IOP(1).EQ.0 .AND. IOP(21).NE.0) RETURN
      CONTAM=DABS(XS-XSE)
      IF(XSE.NE.ZERO)CONTAM=CONTAM/XSE
      IF(CONTAM.GE.PT05)GO TO 300
      IF(CONTAM.GT.PT01)WRITE(IOUT,1010)
      RETURN
C
  300 WRITE(IOUT,1020)
      RETURN
C
C     REPLACE CALCULATED SPIN WITH EXACT SPIN FOR MC-SCF/GVB.
C
  400 XS = XSE
      SP = SPE
      IF (IOP(1).EQ.0 .AND. IOP(21).NE.0) RETURN
      WRITE (IOUT,1030)
      RETURN
      END
      FUNCTION TRAOPN(A,B,N)
C*
C     --------------
C     QCPE VERSION
C     JANUARY 1987
C     --------------
C*
C     TRACE OF A PRODUCT OF DOUBLE SYMMETRIC MATRICES A AND B STORED
C     AS LINEAR VECTORS
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NB=#NB)
C##
      PARAMETER (NB=200)
C###
      PARAMETER (NB1=NB*(NB+1))
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
C
      DIMENSION A(NB1),B(NB1)
C
      TRAOPN=ZERO
      K=0
      DO 3 L=1,2
      DO 2 J=1,N
      DO 1 I=1,J
      K=K+1
    1 TRAOPN=TRAOPN+TWO*A(K)*B(K)
    2 TRAOPN=TRAOPN-A(K)*B(K)
    3 CONTINUE
      RETURN
      END
