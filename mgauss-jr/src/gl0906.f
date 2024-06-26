C     GL0906       20 APR 87                                         MRP
C?IBM/GLD/GBR/VAX/UNX
      SUBROUTINE CIROOT
C??
C?CDC
C     PROGRAM CIROOT
C??
C1INTRODUCTION
C     LINK 0906
C
C     AUTHORS: N.C. HANDY, J.D. GODDARD AND H.F. SCHAEFER, BERKELEY.
C     MODIFIED BY M.R. PETERSON AND R.A. POIRIER,
C     UNIVERSITY OF TORONTO CHEMISTRY DEPT., TORONTO, ONTARIO, CANADA.
C     OPEN SHELL SINGLET ADDED BY P.S. MARTIN AND M.R. PETERSON,
C     UNIVERSITY OF TORONTO CHEMISTRY DEPT., TORONTO, ONTARIO, CANADA.
C     VERSION: APRIL 1987.
C
C     ITERATIVE SOLUTION FOR THE LOWEST EIGENVALUE AND EIGENVECTOR BY
C     DAVIDSON'S SCHEME: E.R. DAVIDSON, J. COMP. PHYS., 17, 87 (1975).
C
C     THE SORTED INTEGRALS (ON UNIT IFINT) ARE READ IN CONJUNCTION WITH
C     THE FORMULA TAPE (UNIT IFORMT). UNIT ISCRCH IS USED FOR SCRATCH.
C
C     THE CURRENT PACKING OF THE FORMULA TAPE GIVES AN ABSOLUTE
C     MAXIMUM OF 32000 CONFIGURATIONS.
C*
C/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA)
C     PARAMETER (NT=#NT, NCON=#NCON )
C##
      PARAMETER (NA= 36)
      PARAMETER (NT= 70, NCON= 32000)
C###
      PARAMETER (NT1=NT*(NT+1)/2)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      COMMON /A/ IOP(99),NATOMS,ICHARG,MULTIP,IC1(NA),NAE,NBE,NE,NBASIS,
     1 CD1(NA,3),CD2(4),IC3(401)
C
      COMMON/C906A/NOCC,NOCC1,NHALF,NMO,NDF,NSF,NSTAT,NDOUB,NBF,NS,
     1 ITY(NT,2),ITV(NT,2),INN(NT,NT),
     2 A(15,15),B(15),C(15),AA(15,15),FILL1,IFILL(2,1680),Y(1680)
      COMMON/C906B/VC(NCON)
      COMMON/C906C/SIG(NCON)
C
      COMMON/GEN/E1,E2,E3,EDUM(3),DCONV,SPIN(2),FIELD(3),FE,EMOL,ESOL,
     1 DOPLE(4),VCM
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,
     1 IMOINT,ISCRCH,IFINT,ILABEL,ICOUPL,IFORMT,IODUM(6),NFILE(100,2)
C
      DIMENSION VEC(15),FF(NT1),ICONF(6)
C
      EQUIVALENCE (ANOCC,NOCC),(FF(1),Y(1))
C
      DATA PT9/0.9D0/, TEN/10.0D0/
      DATA A5TM7/5.0D-7/
      DATA MAXIT/15/
C?IBM/GLD/GBR/VAX/UNX
      DATA IINT/1/
C??
C?CDC
C     DATA IINT/2/
C??
C     MAXIT IS THE MAXIMUM NUMBER OF ITERATIONS.
C     A, AA, B, C AND VEC ARE DIMENSIONED MAXIT.
C*
 1000 FORMAT('1',20('*')/' *CI DIAGONALIZATION*'/1X,20('*')///)
 1010 FORMAT('0MAXIMUM NUMBER OF ITERATIONS EXCEEDED IN CIROOT')
 1020 FORMAT(/'0TOTAL SCF ENERGY',19X,F18.7/' CORRELATION ENERGY',17X,
     1 F18.7/'0TOTAL CI ENERGY',20X,F18.7/'0TOTAL CI ENERGY + ',
     2 '(1-C0*C0)*ECORR',F20.7/'0WHERE (1-C0*C0)*ECORR IS AN ',
     3 'APPROXIMATE CORRECTION TO THE CI ENERGY FOR QUADRUPLE ',
     4 'EXCITATIONS')
 1030 FORMAT(' ITERATION',5X,'CORRELATION ENERGY'/)
 1040 FORMAT(3X,I4,F22.7)
 1050 FORMAT('0C0 IS THE COEFFICIENT OF THE ROOT CONFIGURATION'/
     1 ' C0 = ',F12.7)
 1060 FORMAT('0WARNING: C0 LESS THAN 0.9, USE A MULTI-REFERENCE CI FOR',
     1 ' RELIABLE RESULTS')
 1070 FORMAT(/'0')
 1080 FORMAT (/'0LIST OF CONFIGURATIONS WITH ABS(COEFFICIENT) AT ',
     1 'LEAST 10.0**(-',I1,')')
 1090 FORMAT (/'0LIST OF ALL CONFIGURATIONS')
 1100 FORMAT ('0',6('   CONF   COEFFICIENT ')/)
 1110 FORMAT (1X,6(I7,1PD15.6))
 1120 FORMAT ((1X,6(I7,1PD15.6)))
C
      IF (IOP(10).EQ.0) RETURN
      REWIND IFINT
      REWIND IFORMT
      REWIND ISCRCH
      IPRINT=IOP(19)
      IF (IPRINT.EQ.0) WRITE (IOUT,1000)
C
C     CALCULATE NUMBER OF INTEGERS TO READ FROM UNIT IMAT - BE CAREFUL
C     OF AN ODD NUMBER OF INTEGERS ON MACHINES THAT PACK 2 INTEGERS
C     INTO A DOUBLE WORD (WE NEED NT*(NT+4)+10 INTEGERS).
C
      L9 = ((NT*(NT+4)+12-IINT)*IINT) / 2
      CALL TREAD (9,ANOCC,L9,1,L9,1,0)
      LZ = INN(NOCC1,NOCC1)
      ICASE = MULTIP
      IF (IOP(11) .EQ. 1) ICASE = 5
C
C     INITIALIZE GUESS VECTOR VC AND CORRECTION VECTOR SIG.
C
      DO 30 I=1,NSTAT
      SIG(I)=ZERO
   30 VC(I)=ZERO
      VC(3)=ONE
      CALL SFOPEN(VC,SIG,NCON,IOP,LZ,ICASE)
      IF (IOP(1) .NE. 0) RETURN
C     WRITE VC AND SIG ALTERNATELY ON THE SCRATCH FILE.
      WRITE (ISCRCH) (VC(I),I=1,NSTAT)
      WRITE (ISCRCH) (SIG(I),I=1,NSTAT)
      REWIND ISCRCH
      ZZ=ZERO
      VEC(1)=ONE
      A(1,1)=ZZ
      MM=1
      IF (IPRINT.EQ.0) WRITE (IOUT,1030)
C
   80 CALL OPQUE(MM,ZZ,VEC,VC,SIG,NCON)
      VC(1)=ZERO
      VC(2)=ZERO
      CALL SFOPEN(VC,SIG,NCON,IOP,LZ,ICASE)
      IF (IOP(1) .NE. 0) RETURN
      MM=MM+1
      WRITE (ISCRCH) (VC(I),I=1,NSTAT)
      WRITE (ISCRCH) (SIG(I),I=1,NSTAT)
      REWIND ISCRCH
C
      DO 15 I=1,MM
C     SKIP VC.
      READ (ISCRCH)
C     RETRIEVE SIG.
      READ (ISCRCH) (SIG(J),J=1,NSTAT)
      Z=ZERO
      DO 10 J=3,NSTAT
   10 Z=Z+SIG(J)*VC(J)
      A(MM,I)=Z
      A(I,MM)=Z
   15 CONTINUE
C
      REWIND ISCRCH
      ZZOLD=ZZ
      CALL TRED12(MM,MAXIT,A,AA,B,C,2)
      ZZ=B(1)
      DO 90 I=1,MM
   90 VEC(I)=AA(I,1)
      IF (IPRINT.EQ.0) WRITE (IOUT,1040)MM,ZZ
      IF (DABS(ZZ-ZZOLD).LE.A5TM7) GO TO 95
C
      IF (MM.LT.MAXIT) GO TO 80
      IOP(1)=-2
      WRITE (IOUT,1010)
C
   95 DO 21 I=1,NSTAT
   21 SIG(I)=ZERO
C
      DO 100 I=1,MM
C     RETRIEVE VC.
      READ (ISCRCH) (VC(J),J=1,NSTAT)
C     SKIP SIG.
      READ (ISCRCH)
      DO 22 J=3,NSTAT
   22 SIG(J)=SIG(J)+VEC(I)*VC(J)
  100 CONTINUE
C
      REWIND ISCRCH
C     NORMALIZE THE CI EXPANSION.
      Z=ZERO
      DO 23 I=3,NSTAT
   23 Z=Z+SIG(I)*SIG(I)
      Z=ONE/DSQRT(Z)
      DO 24 I=3,NSTAT
   24 SIG(I)=SIG(I)*Z
      E4=E3+ZZ
      E5=E4+(ONE-SIG(3)*SIG(3))*ZZ
C
      IF (IPRINT-1) 120,110,130
  110 WRITE (IOUT,1070)
  120 WRITE (IOUT,1020) E3,ZZ,E4,E5
      WRITE (IOUT,1050)SIG(3)
C     SAVE CI ENERGY IF IOP(18) IS SET.
  130 IF (IOP(18).EQ.1) E3=E4
      IF (IOP(18).EQ.2) E3=E5
C     SET ENERGY/PROPERTY ORIGIN FLAG TO THE CORRECT TYPE OF CI.
      IOP(94) = IOP(18)
C     CHECK FOR A LARGE LEADING TERM IN THE CI EXPANSION.
      IF (DABS(SIG(3)) .LT. PT9) WRITE (IOUT,1060)
C
C     PRINT THE CI COEFFICIENTS IF NECESSARY.
C
      IF (IOP(21) .EQ. 0) RETURN
      IF (IOP(21) .LT. 7) GO TO 300
C
C     PRINT ALL CI COEFFICIENTS.
C
      WRITE (IOUT,1090)
      WRITE (IOUT,1100)
      WRITE (IOUT,1120) (I,SIG(I),I=1,NSTAT)
      RETURN
C
C     PRINT ONLY CI COEFFICIENTS AT LEAST 10.0**(-IOP(21)).
C
  300 WRITE (IOUT,1080) IOP(21)
      WRITE (IOUT,1100)
      IQ = 1
      CUTOFF = TEN ** (-IOP(21))
C
      DO 310 I=1,NSTAT
      IF (DABS(SIG(I)) .LT. CUTOFF) GO TO 310
      ICONF(IQ) = I
      VEC(IQ) = SIG(I)
      IQ = IQ + 1
      IF (IQ .GT. 6) THEN
         WRITE (IOUT,1110) (ICONF(J),VEC(J),J=1,6)
         IQ = 1
      END IF
  310 CONTINUE
C
C     PRINT ANY REMAINING COEFFICIENTS IN THE LINE BUFFER.
C
      IQ = IQ - 1
      IF (IQ .GT. 0) WRITE (IOUT,1110) (ICONF(J),VEC(J),J=1,IQ)
      RETURN
      END
      SUBROUTINE SFOPEN(VC,SIG,MXCONF,IOP,LZ,ICASE)
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NT=#NT)
C##
      PARAMETER (NT= 70)
C###
      PARAMETER (NT1=NT*(NT+1)/2)
C
      COMMON/C906A/NOCC,NOCC1,NHALF,NMO,NDF,NSF,NSTAT,NDOUB,NBF,NS,
     1 ITY(NT,2),ITV(NT,2),INN(NT,NT),
     2 A(15,15),X(128),V(128),I3(2,1680),Y(1680)
C
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,
     1 IMOINT,ISCRCH,IFINT,ILABEL,ICOUPL,IFORMT,IODUM(6),NFILE(100,2)
C
      DIMENSION FF(NT1),VC(MXCONF),SIG(MXCONF),IOP(99)
C?GLD/GBR/IBM/VAX/UNX
      DIMENSION IX(256)
C??
C?CDC
C     DIMENSION IX(128)
C??
      EQUIVALENCE (X(1),IX(1))
      EQUIVALENCE (FF(1),Y(1))
C*
C     READ THE SORTED MO INTEGRALS UNDER CONTROL OF THE FORMULA TAPE.
      NEXT = 1680
C     SKIP RECORD OF EIGENVALUE DIFFERENCES.
      READ (IFINT)
C?GLD/GBR/IBM/VAX/UNX
      NEED = 256
C??
C?CDC
C     NEED = 128
C??
      CALL NWIORD (IFORMT, IX, NEED, .TRUE., 'SFOPEN',
     1 IOP(1), *9000)
      NEED = 3360
C
   19 CALL NWIORD (IFORMT, I3, NEED, .TRUE., 'SFOPEN',
     1 IOP(1), *9000)
      MT=0
   20 MT=MT+1
      CALL CIUNPK(I3(1,MT),N1,M1)
      IF (M1.GT.32000) GO TO 30
      M2=I3(2,MT)
      SIG(M1)=SIG(M1)+V(N1)*VC(M2)
      IF (M1.EQ.M2) GO TO 21
      SIG(M2)=SIG(M2)+V(N1)*VC(M1)
   21 IF (MT.LT.1680) GO TO 20
      GO TO 19
C
   30 IF (M1.EQ.32002) GO TO 32
      NEXT = NEXT + 1
      IF (NEXT .LE. 1680) GO TO 31
      READ (IFINT) Y
      NEXT = 1
C
   31 U=Y(NEXT)
      DO 70 I=1,64
      V(I)=U*X(I)
   70 V(I+64)=-V(I)
      GO TO 21
C
C     READ THE FOCK MATRICES.
C
   32 READ (IFINT)FF
   40 MT=MT+1
      IF (MT.LE.1680) GO TO 60
      CALL NWIORD (IFORMT, I3, NEED, .TRUE., 'SFOPEN',
     1 IOP(1), *9000)
      MT=1
   60 CALL CIUNPK(I3(1,MT),N1,M1)
      IF (M1.GT.32000) GO TO 50
      CALL CIUNPK(I3(2,MT),N2,M2)
      SIG(M1)=SIG(M1)+FF(N2)*X(N1)*VC(M2)
      IF (M1.EQ.M2) GO TO 40
      SIG(M2)=SIG(M2)+FF(N2)*X(N1)*VC(M1)
      GO TO 40
C
   50 IF (M1.NE.32003) GO TO 32
      IF (ICASE .NE. 5) GO TO 100
C
C     OPEN SHELL SINGLET ONLY.
C
      Z = FF(LZ)
C
   59 MT=MT+1
      IF (MT.LE.1680) GO TO 58
      CALL NWIORD (IFORMT, I3, NEED, .TRUE., 'SFOPEN',
     1 IOP(1), *9000)
      MT=1
C
   58 CALL CIUNPK (I3(1,MT),N1,M1)
      IF (M1 .GT. 32000) GO TO 100
      SIG(M1) = SIG(M1) - Z*VC(M1)*X(N1)
      GO TO 59
C
  100 REWIND IFINT
      REWIND IFORMT
 9000 RETURN
      END
      SUBROUTINE CIUNPK(IPK,N,M)
C?IBM/GLD/GBR/UNX USE HALF-WORD INTEGER UNPACK
      INTEGER*2 IH1,IH2
C*
      COMMON/PACKED/IFILL(4),IH1,IH2
C*
      EQUIVALENCE (IH,IH1)
C*
C*** JR *** ON i386 packing order is the inverse (use VAX)
C      IH=IPK
C      N=IH1
C      M=IH2
C??
C?VAX USE HALF-WORD INTEGER UNPACK
C     INTEGER*2 IH1,IH2
C*
C     COMMON/PACKED/IFILL(4),IH1,IH2
C*
C     EQUIVALENCE (IH,IH1)
C*
C*** JR *** On i386 this code should be used instead of UNX
      IH=IPK
      N=IH2
      M=IH1
C??
C?CDC SHIFT AND AND TO UNPACK.
C     DATA IMASK/65535/
C*
C     N=SHIFT(IPK,-16)
C     M=AND(IPK,IMASK)
C??
      RETURN
      END
      SUBROUTINE OPQUE(MM,ZZ,VEC,VC,SIG,MXCONF)
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NT=#NT)
C##
      PARAMETER (NT= 70)
C###
      PARAMETER (NT1=NT*(NT+1)/2)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      COMMON/C906A/NOCC,NOCC1,NHALF,NMO,NDF,NSF,NSTAT,NDOUB,NBF,NS,
     1 ITY(NT,2),ITV(NT,2),INN(NT,NT),
     2 A(15,15),B(15),C(15),AA(15,15),FILL1,IFILL(2,1680),Y(1680)
C
      COMMON/IO/IODUM1(5),IMOINT,ISCRCH,IFINT,ILABEL,ICOUPL,IFORMT,
     1 IODUM2(206)
C
      DIMENSION VC(MXCONF),SIG(MXCONF),VEC(15),FF(NT1)
C
      EQUIVALENCE (FF(1),Y(1))
C
      DATA TENM3,PT1/1.0D-3,0.1D0/, THRES/1.0D-25/
C*
      DO 10 I=3,NSTAT
   10 SIG(I)=ZERO
      DO 25 I1=1,MM
C     RETRIEVE VC (PREVIOUS VECTORS).
      READ (ISCRCH) (VC(I),I=1,NSTAT)
      Z=ZZ*VEC(I1)
      DO 15 I=3,NSTAT
   15 SIG(I)=SIG(I)-Z*VC(I)
C     RETRIEVE SIG (PREVIOUS SIGMAS).
      READ (ISCRCH) (VC(I),I=1,NSTAT)
      DO 20 I=3,NSTAT
   20 SIG(I)=SIG(I)+VEC(I1)*VC(I)
   25 CONTINUE
      REWIND ISCRCH
      LT=3
      IF (MM.GT.1) GO TO 30
      LT=4
      SIG(3)=ZERO
C     RETRIEVE APPROXIMATE DIAGONAL ELEMENTS FORMED IN LINK 0903.
   30 READ (IFINT) (VC(I),I=1,NSTAT)
      REWIND IFINT
      DO 35 I=LT,NSTAT
      ZZ1=ZZ-VC(I)
C     WATCH FOR ACCIDENTAL ZEROES IN THE DENOMINATOR - THIS AFFECTS
C     ONLY THE RATE OF CONVERGENCE, NOT THE FINAL ENERGY.
      IF (DABS(ZZ1).LT.TENM3) ZZ1=PT1
   35 SIG(I)=SIG(I)/ZZ1
      DO 50 I1=1,MM
C     RETRIEVE VC
      READ (ISCRCH) (VC(I),I=1,NSTAT)
C     SKIP SIG
      READ (ISCRCH)
      Z=ZERO
      DO 40 I=3,NSTAT
   40 Z=Z+SIG(I)*VC(I)
      DO 45 I=3,NSTAT
   45 SIG(I)=SIG(I)-VC(I)*Z
   50 CONTINUE
      Z=ZERO
      DO 55 I=3,NSTAT
   55 Z=Z+SIG(I)*SIG(I)
C
C     WATCH FOR SIG VECTORS THAT ARE ALL ESSENTIALLY ZERO.
C
      IF (Z .GE. THRES) THEN
         Z = ONE / DSQRT(Z)
      ELSE
         Z = ZERO
      END IF
      DO 60 I=3,NSTAT
      VC(I)=SIG(I)*Z
   60 SIG(I)=ZERO
      RETURN
      END
