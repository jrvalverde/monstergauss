C     GL0902A      04 MAY 89                                         MRP
C?IBM/GLD/GBR/VAX/UNX
      SUBROUTINE COUPLE
C??
C?CDC
C     PROGRAM COUPLE
C??
C1INTRODUCTION
C     LINK 0902
C
C     THIS DIRECT CI PROGRAM IS DESCRIBED IN N.C. HANDY, J.D. GODDARD
C     AND H.F. SCHAEFER III, J. CHEM. PHYS., 71, 426 (1979).
C
C     CALCULATE CI COUPLING COEFFICIENTS FOR THIS CASE.
C
C     AUTHORS: N.C. HANDY, J.D. GODDARD AND H.F. SCHAEFER, BERKELEY.
C     MODIFIED: MIKE PETERSON, U OF T CHEMISTRY DEPT, TORONTO, CANADA.
C     OPEN SHELL SINGLET ADDED BY PETER MARTIN AND MIKE PETERSON,
C     UNIVERSITY OF TORONTO CHEMISTRY DEPT, TORONTO, CANADA.
C
C     THE COUPLING COEFFICIENTS ARE WRITTEN ON UNIT ICOUPL.
C*
C/
      IMPLICIT INTEGER(X), DOUBLE PRECISION (A-H,O-W,Y,Z)
C#
C     PARAMETER (NA=#NA)
C##
      PARAMETER (NA= 36)
C###
      REAL C,C1,C2
C
      COMMON /A/ IOP(99),IC1(2),MULTIP,IC2(NA),IC3(4),CD1(NA,3),CD2(4),
     1 IC4(401)
C
      COMMON/C902A/XJ1(24,24,2),XJ2(24,24,2),XJ3(24,24,2),XJ4(24,24,2),
     1 XJ5(24,24,2),XJ6(24,24,2),XK1(24,24,2),XK2(24,24,2),XK3(24,24,2),
     2 XK4(24,24,2),XK5(24,24,2),XK6(24,24,2)
      COMMON/C902B/XM1(6,6),XM2(6,6),XN1(6,6,3),XN2(6,6,3)
     1             ,XR1(24),XR2(24),XS1(6,3)
     2            ,XG1(24,6),XH1(24,6),XG2(24,6),XH2(24,6),XG3(24,6),
     3 XH3(24,6),XG4(24,6),XH4(24,6),XI1(24,6,3),XI2(24,6,3),XI3(24,6,3)
     4 ,XI4(24,6,3)
     5                                      ,XL1(24,24,2),XL1P(24,24,2),
     6 XL2(24,24,2),XL2P(24,24,2),XL3(24,24,2),XL3P(24,24,2),
     7 XL4(24,24,2),XL4P(24,24,2),XL1PP(24,24,2),XL2PP(24,24,2),
     8 XL3PP(24,24,2),XL4PP(24,24,2)
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2        ,NCD(24),NCS(6)
     3        ,Y(128),J7
      COMMON/C902D/Z1(12),Z2(12),C1(12,12),C2(12,12),ICUP(24,24,4)
C     PAD IE123 BY 3 WORDS TO ALIGN THE DOUBLE PRECISION VARIABLES.
      COMMON/C902E/INT123(30,10),Z1E123(120),IE123(128)
      COMMON/IO/IN,IOUT,IODUM(215)
C*
      DIMENSION XC902A(13824),XC902B(17058)
C
      EQUIVALENCE (XC902A(1),XJ1(1,1,1)),(XC902B(1),XM1(1,1))
C*
 1000 FORMAT('0CASE',I4,' CANNOT BE TREATED IN CI PACKAGE')
C*
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,THREE/3.0D0/,FOUR/4.0D0/,
     1 FIVE/5.0D0/,SIX/6.0D0/,F12/12.0D0/,PT25/0.25D0/,
     2 PT3/0.3D0/,PT5/0.5D0/,PT75/0.75D0/,ONEP5/1.5D0/,SEVNP5/7.5D0/,
     3 PT075/0.075D0/,PT1/0.1D0/,PT15/0.15D0/,PT2/0.2D0/,PT4/0.4D0/,
     4 PT45/0.45D0/,PT6/0.6D0/,F15/15.0D0/,F20/20.0D0/,F30/30.0D0/,
     5 F90/90.0D0/
C*
      ICASE=MULTIP
      IF(ICASE.LT.1.OR.ICASE.GT.4)GO TO 100
      IF(IOP(10).EQ.0)RETURN
C     SET ICASE TO 5 FOR THE OPEN SHELL SINGLET CASE.
      IF (IOP(11) .EQ. 1) ICASE = 5
C
C     ZERO COMMONS C902A AND C902B.
C
      DO 200 I=1,13824
  200 XC902A(I)=0
      DO 210 I=1,17058
  210 XC902B(I)=0
C
C     CLEAR COMMON C902C.
C
      DO 220 I=1,12
  220 C(I) = 0.0
      DO 230 I=1,26
  230 NDED(I) = 0
      DO 240 I=1,6
      NDES(I) = 0
  240 NCS(I) = 0
      DO 250 I=1,4
  250 NDER(I) = 0
      DO 260 I=1,10
      MCT(I) = 0
  260 MCTS(I) = 0
      DO 270 I=1,24
  270 NCD(I) = 0
C
C     INITIALIZE COUPLING CONSTANTS ARRAY.
C
      DO 10 I=51,64
   10 Y(I) = ZERO
      J7 = 50
      Y(1) = ZERO
      Y(2) = ONE
      Y(3) = TWO
      Y(4) = FOUR
      Y(5) = DSQRT(TWO)
      Y(6) = DSQRT(THREE)
      Y(7) = Y(5) + Y(5)
      Y(8) = Y(5) / THREE
      Y(9) = Y(5) / SIX
      Y(10) = ONE / Y(3)
      Y(11) = ONE / Y(4)
      Y(12) = ONE / Y(5)
      Y(13) = ONE / Y(6)
      Y(14) = ONE / Y(7)
      Y(15) = ONE / DSQRT(SIX)
      Y(16) = ONE / DSQRT(F12)
      Y(17) = ONEP5
      Y(18) = ONE / Y(17)
      Y(19) = DSQRT(ONEP5)
      Y(20) = ONE / Y(19)
      Y(21) = PT5 * Y(6)
      Y(22) = Y(5) / F12
      Y(23) = ONE / SIX
      Y(24) = ONE / THREE
      Y(25) = PT25 * Y(6)
      Y(26) = PT5 * Y(19)
      Y(27) = PT75 * Y(5)
      Y(28) = TWO * Y(19) / THREE
      Y(29) = Y(19) / SIX
      Y(30) = FIVE * Y(23)
      Y(31) = PT25 / Y(6)
      Y(32) = ONE / DSQRT(SEVNP5)
      Y(33) = DSQRT(PT3)
      Y(34) = DSQRT(PT2)
      Y(35) = Y(34) / THREE
      Y(36) = ONE / DSQRT(F15)
      Y(37) = ONE / SEVNP5
      Y(38) = Y(35) + Y(35)
      Y(39) = DSQRT(PT4) / THREE
      Y(40) = ONE / DSQRT(F90)
      Y(41) = ONE / DSQRT(F30)
      Y(42) = DSQRT(PT15)
      Y(43) = ONE / DSQRT(F20)
      Y(44) = DSQRT(PT45)
      Y(45) = DSQRT(PT1)
      Y(46) = PT2
      Y(47) = PT5 * Y(45)
      Y(48) = DSQRT(PT075)
      Y(49) = PT6
      Y(50) = PT3
C
      DO 20 J=1,64
   20 Y(J+64) = -Y(J)
C
      GO TO (1,2,3,4,5),ICASE
    1 CALL SETUP1
      GO TO 30
    2 CALL SETUP2
      GO TO 30
    3 CALL SETUP3
      GO TO 30
    4 CALL SETUP4
      GO TO 30
    5 CALL SETUP5
C
   30 CALL LHSRHS
      RETURN
C
  100 WRITE(IOUT,1000)ICASE
      IOP(1) = -2
      RETURN
      END
      SUBROUTINE SETUP1
C
C     SETUP FOR CLOSED SHELL SINGLET CASE
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2        ,NCD(24),NCS(6)
     3 ,FILL(128),J7
C
      NCS(1)=1
      NCD(1)=1
      NCD(2)=1
      NCD(3)=1
      NCD(4)=2
      NCT=1
      MCT(1)=4
      MCTS(1)=2
      NCTS=1
      NSC=1
      NDC=5
      MNE=8
      NSO = 0
      NO=4
      NDER(1)=1
      NDES(1)=2
      NDED(1)=1
      NDED(2)=2
      NDED(3)=2
      NDED(4)=6
      NDED(5)=4
      NDED(6)=1
C
C     TEST EXAMPLE ROOT CONFIGURATION
C     .1 STANDS FOR ALPHA SPIN, .2 FOR BETA SPIN.
C     VARIABLES CONTAINING THESE NUMBERS ARE ALL REAL*4 ON THE IBM,
C     GOULD AND VAX COMPUTERS.
C
      C(1)=1.1
      C(2)=1.2
      C(3)=2.1
      C(4)=2.2
      C(5)=3.1
      C(6)=3.2
      C(7)=4.1
      C(8)=4.2
C
C     VIRTUALS 5 THROUGH 9.
C
      MNX1=5
      MNX2=9
      RETURN
      END
      SUBROUTINE SETUP2
C
C     SETUP FOR OPEN SHELL DOUBLET CASE.
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2        ,NCD(24),NCS(6)
     3 ,FILL(128),J7
C
      NCS(1)=1
      NCS(2)=1
      NCS(3)=2
      NCD(1)=1
      NCD(2)=1
      NCD(3)=1
      NCD(4)=1
      NCD(5)=1
      NCD(6)=1
      NCD(7)=2
      NCD(9)=2
      NCD(11)=2
      NDC=12
      NSC=4
      MNE=9
      NDES(1)=1
      NDES(2)=1
      NDES(3)=2
      NDES(4)=2
      NSO=1
      NO=5
      NCTS=1
      MCTS(1)=3
      NCT=3
      MCT(1)=7
      MCT(2)=9
      MCT(3)=11
      NDED(1)=1
      NDED(2)=2
      NDED(3)=2
      NDED(4)=1
      NDED(5)=1
      NDED(6)=1
      NDED(7)=6
      NDED(8)=4
      NDED(9)=3
      NDED(10)=2
      NDED(11)=3
      NDED(12)=2
      NDED(13)=1
      C(1)=1.1
      C(2)=1.2
      C(3)=2.1
      C(4)=2.2
      C(5)=3.1
      C(6)=3.2
      C(7)=4.1
      C(8)=4.2
      C(9)=5.1
      MNX1=5
      MNX2=9
      RETURN
      END
      SUBROUTINE SETUP3
C
C     SETUP FOR OPEN SHELL TRIPLET CASE.
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2        ,NCD(24),NCS(6)
     3 ,FILL(128),J7
C
      NCS(1)=1
      NCS(2)=1
      NCS(3)=2
      NCD(1)=1
      NCD(2)=1
      NCD(3)=1
      NCD(4)=1
      NCD(5)=1
      NCD(6)=1
      NCD(7)=1
      NCD(8)=1
      NCD(9)=2
      NCD(11)=2
      NCD(13)=2
      NCD(15)=2
      NCD(17)=2
      NDC=18
      NSC=4
      MNE=10
      NDES(1)=1
      NDES(2)=1
      NDES(3)=2
      NDES(4)=2
      NSO=2
      NO=6
      NCTS=1
      MCTS(1)=3
      NCT=5
      MCT(1)=9
      MCT(2)=11
      MCT(3)=13
      MCT(4)=15
      MCT(5)=17
      NDED(1)=1
      NDED(2)=2
      NDED(3)=2
      NDED(4)=1
      NDED(5)=1
      NDED(6)=1
      NDED(7)=1
      NDED(8)=1
      NDED(9)=6
      NDED(10)=4
      NDED(11)=3
      NDED(12)=2
      NDED(13)=3
      NDED(14)=2
      NDED(15)=1
      NDED(16)=1
      NDED(17)=1
      NDED(18)=1
      NDED(19)=1
      C(1)=1.1
      C(2)=1.2
      C(3)=2.1
      C(4)=2.2
      C(5)=3.1
      C(6)=3.2
      C(7)=4.1
      C(8)=4.2
      C(9)=5.1
      C(10)=6.1
      MNX1=5
      MNX2=11
      RETURN
      END
      SUBROUTINE SETUP4
C
C     SETUP FOR OPEN SHELL QUARTET CASE.
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2        ,NCD(24),NCS(6)
     3 ,FILL(128),J7
C
      NCS(1)=1
      NCS(2)=1
      NCS(3)=2
      NCD(1)=1
      NCD(2)=1
      NCD(3)=1
      NCD(4)=1
      NCD(5)=1
      NCD(6)=1
      NCD(7)=1
      NCD(8)=1
      NCD(9)=2
      NCD(11)=2
      NCD(13)=2
      NCD(15)=2
      NCD(17)=3
      NCD(20)=3
      NDC=22
      NSC=4
      MNE=11
      NDES(1)=1
      NDES(2)=1
      NDES(3)=2
      NDES(4)=2
      NSO=3
      NO=7
      NCTS=1
      MCTS(1)=3
      NCT=8
      MCT(1)=9
      MCT(2)=11
      MCT(3)=13
      MCT(4)=15
      MCT(5)=17
      MCT(6)=18
      MCT(7)=20
      MCT(8)=21
      NDED(1)=1
      NDED(2)=2
      NDED(3)=2
      NDED(4)=1
      NDED(5)=1
      NDED(6)=1
      NDED(7)=1
      NDED(8)=1
      NDED(9)=6
      NDED(10)=4
      NDED(11)=3
      NDED(12)=2
      NDED(13)=3
      NDED(14)=2
      NDED(15)=1
      NDED(16)=1
      NDED(17)=1
      NDED(18)=1
      NDED(19)=1
      NDED(20)=1
      NDED(21)=1
      NDED(22)=1
      NDED(23)=1
      C(1)=1.1
      C(2)=1.2
      C(3)=2.1
      C(4)=2.2
      C(5)=3.1
      C(6)=3.2
      C(7)=4.1
      C(8)=4.2
      C(9)=5.1
      C(10)=6.1
      C(11)=7.1
      MNX1=5
      MNX2=11
      RETURN
      END
      SUBROUTINE SETUP5
C
C     SETUP FOR OPEN SHELL SINGLET CASE
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2        ,NCD(24),NCS(6)
     3 ,FILL(128),J7
C
      NDES(1)=2
      NDES(2)=2
      NDES(3)=2
      NDES(4)=2
      NDES(5)=4
      NDES(6)=4
      NCS(1)=1
      NCS(2)=1
      NCS(3)=1
      NCS(4)=1
      NCS(5)=2
      NCD(1)=1
      NCD(2)=1
      NCD(3)=1
      NCD(4)=1
      NCD(5)=1
      NCD(6)=1
      NCD(7)=1
      NCD(8)=1
      NCD(9)=1
      NCD(10)=1
      NCD(11)=1
      NCD(12)=2
      NCD(14)=2
      NCD(16)=2
      NCD(18)=2
      NCD(20)=2
      NCD(22)=1
      NCD(23)=1
      NDED(1)=2
      NDED(2)=4
      NDED(3)=4
      NDED(4)=2
      NDED(5)=2
      NDED(6)=2
      NDED(7)=2
      NDED(8)=2
      NDED(9)=2
      NDED(10)=2
      NDED(11)=2
      NDED(12)=6
      NDED(13)=4
      NDED(14)=6
      NDED(15)=4
      NDED(16)=6
      NDED(17)=4
      NDED(18)=6
      NDED(19)=4
      NDED(20)=12
      NDED(21)=8
      NDED(22)=2
      NDED(23)=2
      NDED(24)=2
      NDC=23
      NSC=6
      MNE=10
      NSO=2
      NO=6
      NCTS=1
      MCTS(1)=5
      NCT=5
      MCT(1)=20
      MCT(2)=12
      MCT(3)=14
      MCT(4)=16
      MCT(5)=18
      MNX1=5
      MNX2=10
      C(1)=1.1
      C(2)=1.2
      C(3)=2.1
      C(4)=2.2
      C(5)=3.1
      C(6)=3.2
      C(7)=4.1
      C(8)=4.2
      C(9)=5.1
      C(10)=6.2
      RETURN
      END
      SUBROUTINE LHSRHS
C#
C     PARAMETER (NA=#NA)
C##
      PARAMETER (NA= 36)
C###
      DOUBLE PRECISION FILL,Z1,Z2,CD1,CD2
C
      COMMON /A/ IOP(99),IC1(NA),IC2(7),CD1(NA,3),CD2(4),IC3(401)
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
      COMMON/C902D/Z1(12),Z2(12),C1(12,12),C2(12,12),IMATR(24,24,2),
     1 IMATR1(24,24,2)
C
      DIMENSION IA(2),IB(2),JA(2),JB(2),KA(2),KB(2),LY(5)
C
      DO 140 J8=1,19
      IF(J8.LT.11)J=J8
      IF(J8.GT.10)J=J8-10
      DO 25 L=1,2
      DO 25 K=1,NDC
      DO 25 M=1,NDC
   25 IMATR(M,K,L)=0
      DO 220 I1=1,NO
      DO 220 I2=I1,NO
      DO 220 I3=MNX1,MNX2
      DO 220 I4=I3,MNX2
      IA(1)=I1
      IA(2)=I2
      IB(1)=I3
      IB(2)=I4
C
      GO TO (1001,1002,1003,1004,1005),ICASE
 1001 CALL TYPDD1(IA,IB,MU,IP,JA,JB,I1,KA,KB)
      GO TO 1006
C
 1002 CALL TYPDD2(IA,IB,MU,IP,JA,JB,I1,KA,KB)
      GO TO 1006
C
 1003 CALL TYPDD3(IA,IB,MU,IP,JA,JB,I1,KA,KB)
      GO TO 1006
C
 1004 CALL TYPDD4(IA,IB,MU,IP,JA,JB,I1,KA,KB)
      GO TO 1006
C
 1005 CALL TYPDD5(IA,IB,MU,IP,JA,JB,I1,KA,KB)
C
 1006 IF(MU.EQ.61)GO TO 220
   30 NDETL=NDED(MU)
C
      GO TO (2001,2002,2003,2004,2005),ICASE
 2001 CALL CODED1(IA,IB,MU,C1,NDETL,Z1)
      GO TO 2006
C
 2002 CALL CODED2(IA,IB,MU,C1,NDETL,Z1)
      GO TO 2006
C
 2003 CALL CODED3(IA,IB,MU,C1,NDETL,Z1)
      GO TO 2006
C
 2004 CALL CODED4(IA,IB,MU,C1,NDETL,Z1)
      GO TO 2006
C
 2005 CALL CODED5(IA,IB,MU,C1,NDETL,Z1)
C
 2006 GO TO(35,45,55,65,75,85,95,105,115,125,150,160,170,180,190,195,
     1 200,205,198),J8
   35 DO 40 K=1,NO
      DO 40 L=K,NO
      JA(1)=K
      JA(2)=L
      JB(1)=IB(1)
      JB(2)=IB(2)
   40 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
   45 DO 50 K=MNX1,MNX2
      DO 50 L=K,MNX2
      JA(1)=IA(1)
      JA(2)=IA(2)
      JB(1)=K
      JB(2)=L
   50 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
   55 DO 60 K=1,NO
      DO 60 L=MNX1,MNX2
      JA(1)=K
      JA(2)=IA(2)
      JB(1)=L
      JB(2)=IB(2)
   60 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
   65 DO 70 K=1,NO
      DO 70 L=MNX1,MNX2
      JA(1)=IA(1)
      JA(2)=K
      JB(1)=IB(1)
      JB(2)=L
   70 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
   75 DO 80 K=1,NO
      DO 80 L=MNX1,MNX2
      JA(1)=K
      JA(2)=IA(2)
      JB(1)=IB(1)
      JB(2)=L
   80 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
   85 DO 90 K=1,NO
      DO 90 L=MNX1,MNX2
      JA(1)=IA(1)
      JA(2)=K
      JB(1)=L
      JB(2)=IB(2)
   90 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
   95 DO 100 K=MNX1,MNX2
      JA(1)=IA(1)
      JA(2)=IA(2)
      JB(1)=IB(1)
      JB(2)=K
  100 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
  105 DO 110 K=MNX1,MNX2
      JA(1)=IA(1)
      JA(2)=IA(2)
      JB(1)=K
      JB(2)=IB(2)
  110 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
  115 DO 120 K=1,NO
      JA(1)=IA(1)
      JA(2)=K
      JB(1)=IB(1)
      JB(2)=IB(2)
  120 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 135
C
  125 DO 130 K=1,NO
      JA(1)=K
      JA(2)=IA(2)
      JB(1)=IB(1)
      JB(2)=IB(2)
  130 CALL SUBR1(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
  135 CONTINUE
      GO TO 210
C
  150 DO 155 K=MNX1,MNX2
      JA(1)=IA(1)
      JB(1)=K
  155 CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 210
C
  160 DO 165 K=MNX1,MNX2
      JA(1)=IA(2)
      JB(1)=K
  165 CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 210
C
  170 DO 175 K=1,NO
      JA(1)=K
      JB(1)=IB(2)
  175 CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 210
C
  180 DO 185 K=1,NO
      JA(1)=K
      JB(1)=IB(1)
  185 CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 210
C
  190 JA(1)=IA(1)
      JB(1)=IB(1)
      CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 210
C
  195 JA(1)=IA(1)
      JB(1)=IB(2)
      CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 210
C
  200 JA(1)=IA(2)
      JB(1)=IB(1)
      CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 210
C
  205 JA(1)=IA(2)
      JB(1)=IB(2)
      CALL SUBR2(MU,C1,NDETL,Z1,J,JA,JB,IA,IB,IMATR)
  210 CONTINUE
      GO TO 197
C
  198 IP=300
      J=1
      NU=NDC+1
      GO TO (3001,3002,3003,3004,3005),ICASE
C
 3001 CALL CODED1(IA,IB,NU,C2,NDETR,Z2)
      GO TO 3006
C
 3002 CALL CODED2(IA,IB,NU,C2,NDETR,Z2)
      GO TO 3006
C
 3003 CALL CODED3(IA,IB,NU,C2,NDETR,Z2)
      GO TO 3006
C
 3004 CALL CODED4(IA,IB,NU,C2,NDETR,Z2)
      GO TO 3006
C
 3005 CALL CODED5(IA,IB,NU,C2,NDETR,Z2)
C
 3006 IP4=0
      IF (ICASE.EQ.5) LY(1)=0
      CALL MATRG(C1,C2,Z1,Z2,NDETL,NDETR,J,MU,NU,IP,IA,IB,JA,JB,IP4,LY)
  197 DO 1313 I=1,NCT
      IF(MU.EQ.MCT(I))GO TO 215
 1313 CONTINUE
      GO TO 220
C
  215 MU=MU+1
      GO TO 30
C
  220 CONTINUE
  140 CONTINUE
C
      DO 259 J8=1,4
      IF(J8.LT.4)J=J8
      DO 15 K=1,NSC
      DO 15 M=1,NSC
   15 IMATR(M,K,1)=0
      DO 265 I1=1,NO
      DO 265 I2=MNX1,MNX2
      IA(1)=I1
      IB(1)=I2
      GO TO (4001,4002,4003,4004,4005),ICASE
C
 4001 CALL TYPES1(IA,IB,MU)
      GO TO 4006
C
 4002 CALL TYPES2(IA,IB,MU)
      GO TO 4006
C
 4003 CALL TYPES3(IA,IB,MU)
      GO TO 4006
C
 4004 CALL TYPES4(IA,IB,MU)
      GO TO 4006
C
 4005 CALL TYPES5(IA,IB,MU)
C
 4006 IF(MU.EQ.62)GO TO 265
  513 GO TO (5001,5002,5003,5004,5005),ICASE
C
 5001 CALL CODES1(IA,IB,MU,C1,NDETSL,Z1)
      GO TO 5006
C
 5002 CALL CODES2(IA,IB,MU,C1,NDETSL,Z1)
      GO TO 5006
C
 5003 CALL CODES3(IA,IB,MU,C1,NDETSL,Z1)
      GO TO 5006
C
 5004 CALL CODES4(IA,IB,MU,C1,NDETSL,Z1)
      GO TO 5006
C
 5005 CALL CODES5(IA,IB,MU,C1,NDETSL,Z1)
C
 5006 GO TO(230,240,250,196),J8
C
  230 DO 235 K=1,NO
      DO 235 L=MNX1,MNX2
      JA(1)=K
      JB(1)=L
  235 CALL SUBR3(MU,C1,NDETSL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 260
C
  240 DO 245 K=MNX1,MNX2
      JA(1)=IA(1)
      JB(1)=K
  245 CALL SUBR3(MU,C1,NDETSL,Z1,J,JA,JB,IA,IB,IMATR)
      GO TO 260
C
  250 DO 255 K=1,NO
      JA(1)=K
      JB(1)=IB(1)
  255 CALL SUBR3(MU,C1,NDETSL,Z1,J,JA,JB,IA,IB,IMATR)
  260 CONTINUE
      GO TO 194
C
  196 IP=300
      J=2
      NU=NDC+1
      GO TO (6001,6002,6003,6004,6005),ICASE
C
 6001 CALL CODED1(IA,IB,NU,C2,NDETR,Z2)
      GO TO 6006
C
 6002 CALL CODED2(IA,IB,NU,C2,NDETR,Z2)
      GO TO 6006
C
 6003 CALL CODED3(IA,IB,NU,C2,NDETR,Z2)
      GO TO 6006
C
 6004 CALL CODED4(IA,IB,NU,C2,NDETR,Z2)
      GO TO 6006
C
 6005 CALL CODED5(IA,IB,NU,C2,NDETR,Z2)
C
 6006 IP4=0
      IF (ICASE.EQ.5) THEN
         LY(1)=1
         LY(2)=IA(1)
         LY(3)=IB(1)
         LY(4)=5
         LY(5)=6
      END IF
      CALL MATRG(C1,C2,Z1,Z2,NDETSL,NDETR,J,MU,NU,IP,IA,IB,JA,JB,IP4,LY)
C
  194 DO 1414 I=1,NCTS
      IF(MU.EQ.MCTS(I))GO TO 1515
 1414 CONTINUE
      GO TO 265
C
 1515 MU=MU+1
      GO TO 513
C
  265 CONTINUE
  259 CONTINUE
C
      CALL OUTPUT (IOP(15))
      CALL WRTOUT
      RETURN
      END
      SUBROUTINE TYPDD1(JA,JB,NU,IP,KA,KB,J,IA,IB)
C
C     CLOSED SHELL SINGLET CASE
C     CONFIGURATION TYPE FOR DOUBLE REPLACEMENT STATES
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION IA(2),IB(2),JA(2),JB(2),KA(2),KB(2)
C
      IP=1
      KA(1)=JA(1)
      KA(2)=JA(2)
      KB(1)=JB(1)
      KB(2)=JB(2)
      IF(JA(2).GE.JA(1))GO TO 10
      IP=-IP
      KA(2)=JA(1)
      KA(1)=JA(2)
   10 IF(JB(2).GE.JB(1))GO TO 15
      IP=-IP
      KB(2)=JB(1)
      KB(1)=JB(2)
   15 IF(J.EQ.2)GO TO 25
      IF(J.EQ.1)GO TO 20
      GO TO 30
C
   20 IF(IA(1).EQ.IA(2).OR.JA(1).EQ.JA(2))GO TO 30
      IF(IA(1).EQ.JA(2))IP=-1
      IF(IA(2).EQ.JA(1))IP=-1
      GO TO 30
C
   25 IF(IB(1).EQ.IB(2).OR.JB(1).EQ.JB(2))GO TO 30
      IF(IB(1).EQ.JB(2))IP=-1
      IF(IB(2).EQ.JB(1))IP=-1
   30 IF(KB(1).NE.KB(2))GO TO 45
      IF(KA(1).NE.KA(2))GO TO 35
C
C     II TO AA
C
      NU=1
      RETURN
C
C     IJ TO AA
C
   35 NU=2
      RETURN
C
   45 IF(KA(1).NE.KA(2))GO TO 55
C
C     II TO AB
C
      NU=3
      RETURN
C
C     IJ TO AB
C
   55 NU=4
      RETURN
      END
      SUBROUTINE TYPDD2(JA,JB,NU,IP,KA,KB,I1,IA,IB)
C
C     OPEN SHELL DOUBLET CASE
C     CONFIGURATION TYPE FOR DOUBLE REPLACEMENT STATES
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION IA(2),IB(2),JA(2),JB(2),KA(2),KB(2)
C
      NDO=NO-NSO
      IP=1
      KA(1)=JA(1)
      KA(2)=JA(2)
      KB(1)=JB(1)
      KB(2)=JB(2)
      IF(JA(2).GE.JA(1))GO TO 10
      IP=-IP
      KA(2)=JA(1)
      KA(1)=JA(2)
   10 IF(JB(2).GE.JB(1))GO TO 15
      IP=-IP
      KB(2)=JB(1)
      KB(1)=JB(2)
   15 IF(I1.EQ.2)GO TO 25
      IF(I1.EQ.1)GO TO 20
      GO TO 30
C
   20 IF(IA(1).EQ.IA(2).OR.JA(1).EQ.JA(2))GO TO 30
      IF(IA(1).EQ.JA(2))IP=-1
      IF(IA(2).EQ.JA(1))IP=-1
      GO TO 30
C
   25 IF(IB(1).EQ.IB(2).OR.JB(1).EQ.JB(2))GO TO 30
      IF(IB(1).EQ.JB(2))IP=-1
      IF(IB(2).EQ.JB(1))IP=-1
   30 IF(JB(1).LE.NO.OR.JB(2).LE.NO)GO TO 70
      IF(KB(1).NE.KB(2))GO TO 45
      IF(KA(1).NE.KA(2))GO TO 35
      IF(KA(1).GT.NDO)GO TO 65
      NU=1
      RETURN
C
   35 IF(KA(1).GT.NDO)GO TO 65
      IF(KA(2).GT.NDO)GO TO 40
      NU=2
      RETURN
C
   40 NU=5
      RETURN
C
   45 IF(KA(1).EQ.KA(2).AND.KA(1).GT.NDO)GO TO 65
      IF(KA(1).EQ.KA(2))GO TO 60
      IF(KA(1).GT.NDO)GO TO 65
      IF(KA(2).LE.NDO)GO TO 50
      NU=11
      RETURN
C
   50 NU=7
      RETURN
C
   60 NU=3
      RETURN
C
   65 NU=61
      RETURN
C
   70 NU=61
      IF(KA(1).EQ.KA(2).AND.KA(1).LE.NDO.AND.KB(2).GT.NO)NU=4
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).GT.NO)NU=9
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.KB(1).AND.KA(2).EQ.
     1 5)NU=6
      RETURN
      END
      SUBROUTINE TYPDD3(JA,JB,NU,IP,KA,KB,I1,IA,IB)
C
C     OPEN SHELL TRIPLET CASE
C     CONFIGURATION TYPE FOR DOUBLE REPLACEMENT STATES
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION IA(2),IB(2),JA(2),JB(2),KA(2),KB(2)
C
      NDO=NO-NSO
      IP=1
      KA(1)=JA(1)
      KA(2)=JA(2)
      KB(1)=JB(1)
      KB(2)=JB(2)
      IF(JA(2).GE.JA(1))GO TO 10
      IP=-IP
      KA(2)=JA(1)
      KA(1)=JA(2)
   10 IF(JB(2).GE.JB(1))GO TO 15
      IP=-IP
      KB(2)=JB(1)
      KB(1)=JB(2)
   15 IF(I1.EQ.2)GO TO 25
      IF(I1.EQ.1)GO TO 20
      GO TO 30
C
   20 IF(IA(1).EQ.IA(2).OR.JA(1).EQ.JA(2))GO TO 30
      IF(IA(1).EQ.JA(2))IP=-1
      IF(IA(2).EQ.JA(1))IP=-1
      GO TO 30
C
   25 IF(IB(1).EQ.IB(2).OR.JB(1).EQ.JB(2))GO TO 30
      IF(IB(1).EQ.JB(2))IP=-1
      IF(IB(2).EQ.JB(1))IP=-1
   30 IF(JB(1).LE.NO.OR.JB(2).LE.NO)GO TO 70
      IF(KB(1).NE.KB(2))GO TO 45
      IF(KA(1).NE.KA(2))GO TO 35
      IF(KA(1).GT.NDO)GO TO 65
      NU=1
      RETURN
C
C     II TO AA
C
   35 IF(KA(1).GT.NDO)GO TO 65
      IF(KA(2).GT.NDO)GO TO 40
      NU=2
      RETURN
C
C     II TO AB
C
   40 NU=5
      RETURN
C
C     II TO AR
C
   45 IF(KA(1).EQ.KA(2).AND.KA(1).GT.NDO)GO TO 65
      IF(KA(1).EQ.KA(2))GO TO 60
      IF(KA(1).GT.NDO)GO TO 55
      IF(KA(2).LE.NDO)GO TO 50
      NU=13
      RETURN
C
C     IJ TO RS
   50 NU=9
      RETURN
C
C     RI TO AA
C
   55 NU=8
      RETURN
C
C     IJ TO AR
C
   60 NU=3
      RETURN
C
C     IJ TO AA
C
   65 NU=61
      RETURN
C
   70 NU=61
      IF(KA(1).EQ.KA(2).AND.KA(1).LE.NDO.AND.KB(2).GT.NO)NU=4
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).GT.NO)NU=11
      IF(KA(1).LE.NDO.AND.KA(2).GT.NDO.AND.KB(2).GT.NO.AND.
     1 KA(2).NE.KB(1))NU=6
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).LE.NO.AND.KB(1).
     1 NE.KB(2))NU=7
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.KB(1).AND.KA(2).EQ.
     1 5)NU=15
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.KB(1).AND.KA(2).EQ.
     1 6)NU=17
      RETURN
      END
      SUBROUTINE TYPDD4(JA,JB,NU,IP,KA,KB,I1,IA,IB)
C
C     OPEN SHELL QUARTET CASE
C     CONFIGURATION TYPE FOR DOUBLE REPLACEMENT STATES
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION IA(2),IB(2),JA(2),JB(2),KA(2),KB(2)
C
      NDO=NO-NSO
      IP=1
      KA(1)=JA(1)
      KA(2)=JA(2)
      KB(1)=JB(1)
      KB(2)=JB(2)
      IF(JA(2).GE.JA(1))GO TO 10
      IP=-IP
      KA(2)=JA(1)
      KA(1)=JA(2)
   10 IF(JB(2).GE.JB(1))GO TO 15
      IP=-IP
      KB(2)=JB(1)
      KB(1)=JB(2)
   15 IF(I1.EQ.2)GO TO 25
      IF(I1.EQ.1)GO TO 20
      GO TO 30
   20 IF(IA(1).EQ.IA(2).OR.JA(1).EQ.JA(2))GO TO 30
      IF(IA(1).EQ.JA(2))IP=-1
      IF(IA(2).EQ.JA(1))IP=-1
      GO TO 30
C
   25 IF(IB(1).EQ.IB(2).OR.JB(1).EQ.JB(2))GO TO 30
      IF(IB(1).EQ.JB(2))IP=-1
      IF(IB(2).EQ.JB(1))IP=-1
   30 IF(JB(1).LE.NO.OR.JB(2).LE.NO)GO TO 70
      IF(KB(1).NE.KB(2))GO TO 45
      IF(KA(1).NE.KA(2))GO TO 35
      IF(KA(1).GT.NDO)GO TO 65
      NU=1
      RETURN
C
   35 IF(KA(1).GT.NDO)GO TO 65
      IF(KA(2).GT.NDO)GO TO 40
      NU=2
      RETURN
C
   40 NU=5
      RETURN
C
   45 IF(KA(1).EQ.KA(2).AND.KA(1).GT.NDO)GO TO 65
      IF(KA(1).EQ.KA(2))GO TO 60
      IF(KA(1).GT.NDO)GO TO 55
      IF(KA(2).LE.NDO)GO TO 50
      NU=13
      RETURN
C
   50 NU=9
      RETURN
C
   55 NU=8
      RETURN
C
   60 NU=3
      RETURN
C
   65 NU=61
      RETURN
C
   70 NU=61
      IF(KA(1).EQ.KA(2).AND.KA(1).LE.NDO.AND.KB(2).GT.NO)NU=4
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).GT.NO)NU=11
      IF(KA(1).LE.NDO.AND.KA(2).GT.NDO.AND.KB(2).GT.NO.AND.
     1 KA(2).NE.KB(1))NU=6
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).LE.NO.AND.KB(1).
     1 NE.KB(2))NU=7
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.KB(1).AND.KA(2).EQ.
     1 5)NU=15
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.KB(1).AND.KA(2).EQ.
     1 6)NU=17
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.KB(1).AND.KA(2).EQ.
     1 7)NU=20
      RETURN
      END
      SUBROUTINE TYPDD5(JA,JB,NU,IP,KA,KB,I1,IA,IB)
C
C     OPEN SHELL SINGLET CASE
C     CONFIGURATION TYPE FOR DOUBLE REPLACEMENT STATES
C
      DOUBLE PRECISION FILL
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION IA(2),IB(2),JA(2),JB(2),KA(2),KB(2)
C
      NDO = NO - NSO
      IP=1
      KA(1)=JA(1)
      KA(2)=JA(2)
      KB(1)=JB(1)
      KB(2)=JB(2)
      IF(JA(2).GE.JA(1))GO TO 10
      IP=-IP
      KA(2)=JA(1)
      KA(1)=JA(2)
   10 IF(JB(2).GE.JB(1))GO TO 15
      IP=-IP
      KB(2)=JB(1)
      KB(1)=JB(2)
   15 IF(I1.EQ.2)GO TO 25
      IF(I1.EQ.1)GO TO 20
      GO TO 30
C
   20 IF(IA(1).EQ.IA(2).OR.JA(1).EQ.JA(2))GO TO 30
      IF(IA(1).EQ.JA(2))IP=-1
      IF(IA(2).EQ.JA(1))IP=-1
      GO TO 30
C
   25 IF(IB(1).EQ.IB(2).OR.JB(1).EQ.JB(2))GO TO 30
      IF(IB(1).EQ.JB(2))IP=-1
      IF(IB(2).EQ.JB(1))IP=-1
C
   30 IF(JB(1).LE.NO.OR.JB(2).LE.NO)GO TO 70
      IF(KB(1).NE.KB(2))GO TO 45
      IF(KA(1).NE.KA(2))GO TO 35
      IF(KA(1).GT.NDO)GO TO 65
      NU=1
      RETURN
C
   35 IF(KA(1).GT.NDO)GO TO 65
      IF(KA(2).GT.NDO)GO TO 40
      NU=2
      RETURN
C
   40 IF(KA(2).EQ.5)NU=6
      IF(KA(2).EQ.6)NU=7
      RETURN
C
   45 IF(KA(1).EQ.KA(2).AND.KA(1).GT.NDO)GO TO 65
      IF(KA(1).EQ.KA(2))GO TO 60
      IF(KA(1).GT.NDO)GO TO 55
      IF(KA(2).LE.NDO)GO TO 50
      IF(KA(2).EQ.5)NU=12
      IF(KA(2).EQ.6)NU=14
      RETURN
C
   50 NU=20
      RETURN
C
   55 NU=9
      RETURN
C
   60 NU=3
      RETURN
C
   65 NU=61
      RETURN
C
   70 NU=61
      IF(KA(1).EQ.KA(2).AND.KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KB(1).EQ.
     1 5)NU=4
      IF(KA(1).EQ.KA(2).AND.KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KB(1).EQ.
     1 6)NU=5
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).GT.NO.AND.KB(1).EQ.
     1 5)NU=16
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).GT.NO.AND.KB(1).EQ.
     1 6)NU=18
      IF(KA(1).LE.NDO.AND.KA(2).EQ.5.AND.KB(2).GT.NO.AND.KB(1).EQ.6)
     1 NU=10
      IF(KA(1).LE.NDO.AND.KA(2).EQ.6.AND.KB(2).GT.NO.AND.KB(1).EQ.5)
     1 NU=11
      IF(KA(1).NE.KA(2).AND.KA(2).LE.NDO.AND.KB(2).EQ.6.AND.KB(1).EQ.
     1 5)NU=8
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.5.AND.KB(1).EQ.5)
     1 NU=22
      IF(KA(1).LE.NDO.AND.KB(2).GT.NO.AND.KA(2).EQ.6.AND.KB(1).EQ.6)
     1 NU=23
      RETURN
      END
      SUBROUTINE CODED1(IA,IB,NST,C2,NDET,Z2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL C,C2
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION C2(12,12),Z2(12),IA(2),IB(2)
C
      DATA PT5/0.5D0/,ONE/1.0D0/,TWO/2.0D0/,THREE/3.0D0/
C
      NDET=NDED(NST)
      DO 10 I=1,MNE
      DO 10 J=1,NDET
   10 C2(J,I)=C(I)
      GO TO(20,30,25,35,40,15),NST
   15 Z2(1)=ONE
      RETURN
C
   20 NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NI+1)=IB(2)+0.2
      Z2(1)=ONE
      RETURN
C
   25 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI-1)=IB(1)+0.1
      C2(1,NI)=IB(2)+0.2
      C2(2,NI-1)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.2
      RETURN
C
   30 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI)=IB(1)+0.2
      C2(2,NI-1)=IB(2)+0.1
      NJ=2*IA(2)
      C2(1,NJ-1)=IB(1)+0.1
      C2(2,NJ)=IB(2)+0.2
      RETURN
C
   35 Z2(1)=PT5/DSQRT(THREE)
      Z2(2)=Z2(1)
      Z2(3)=-Z2(1)
      Z2(4)=-Z2(1)
      Z2(5)=Z2(1)+Z2(1)
      Z2(6)=Z2(5)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.1
      C2(2,NJ+1)=IB(2)+0.2
      C2(3,NI+1)=IB(1)+0.2
      C2(3,NJ+1)=IB(2)+0.2
      C2(4,NI)=IB(1)+0.1
      C2(4,NJ)=IB(2)+0.1
      C2(5,NI)=IB(2)+0.1
      C2(5,NJ+1)=IB(1)+0.2
      C2(6,NI+1)=IB(2)+0.2
      C2(6,NJ)=IB(1)+0.1
      RETURN
C
   40 Z2(1)=PT5
      Z2(2)=PT5
      Z2(3)=PT5
      Z2(4)=PT5
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ+1)=IB(2)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(1)+0.1
      C2(3,NJ)=IB(2)+0.1
      C2(4,NI+1)=IB(1)+0.2
      C2(4,NJ+1)=IB(2)+0.2
      RETURN
      END
      SUBROUTINE CODED2(IA,IB,NSI,C2,NDET,Z2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL C,C2
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION C2(12,12),Z2(12),IA(2),IB(2)
C
      DATA PT5/0.5D0/,ONE/1.0D0/,TWO/2.0D0/,THREE/3.0D0/,SIX/6.0D0/
C
      NDO=NO-NSO
      NDET=NDED(NSI)
      DO 10 I=1,MNE
      DO 10 J=1,NDET
   10 C2(J,I)=C(I)
      GO TO(20,25,30,35,40,45,50,55,60,65,70,75,15),NSI
   15 Z2(1)=ONE
      RETURN
C
   20 NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NI+1)=IB(2)+0.2
      Z2(1)=ONE
      RETURN
   25 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI)=IB(1)+0.2
      C2(2,NI-1)=IB(2)+0.1
      NJ=2*IA(2)
      C2(1,NJ-1)=IB(1)+0.1
      C2(2,NJ)=IB(2)+0.2
      RETURN
C
   30 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI-1)=IB(1)+0.1
      C2(1,NI)=IB(2)+0.2
      C2(2,NI-1)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.2
      RETURN
C
   50 Z2(1)=PT5/DSQRT(THREE)
      Z2(2)=Z2(1)
      Z2(3)=-Z2(1)
      Z2(4)=-Z2(1)
      Z2(5)=Z2(1)+Z2(1)
      Z2(6)=Z2(5)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.1
      C2(2,NJ+1)=IB(2)+0.2
      C2(3,NI+1)=IB(1)+0.2
      C2(3,NJ+1)=IB(2)+0.2
      C2(4,NI)=IB(1)+0.1
      C2(4,NJ)=IB(2)+0.1
      C2(5,NI)=IB(2)+0.1
      C2(5,NJ+1)=IB(1)+0.2
      C2(6,NI+1)=IB(2)+0.2
      C2(6,NJ)=IB(1)+0.1
      RETURN
C
   55 Z2(1)=PT5
      Z2(2)=PT5
      Z2(3)=PT5
      Z2(4)=PT5
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ+1)=IB(2)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(1)+0.1
      C2(3,NJ)=IB(2)+0.1
      C2(4,NI+1)=IB(1)+0.2
      C2(4,NJ+1)=IB(2)+0.2
      RETURN
C
   40 Z2(1)=ONE
      NI=2*IA(1)
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      RETURN
C
   70 Z2(2)=ONE/DSQRT(SIX)
      Z2(1)=Z2(2)+Z2(2)
      Z2(3)=-Z2(2)
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI+1)=IB(2)+0.2
      C2(3,NJ)=IB(1)+0.1
      RETURN
C
   75 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI+1)=IB(2)+0.2
      C2(2,NJ)=IB(1)+0.1
      RETURN
C
   35 Z2(1)=ONE
      NI=2*IA(1)-1
      C2(1,NI)=IB(2)+0.1
      C2(1,NI+1)=IB(1)+0.2
      RETURN
C
   60 Z2(2)=ONE/DSQRT(SIX)
      Z2(1)=Z2(2)+Z2(2)
      Z2(3)=-Z2(2)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ+1)=IB(2)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(2)+0.1
      C2(3,NJ+1)=IB(1)+0.2
      RETURN
C
   65 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(2)+0.1
      C2(2,NJ+1)=IB(1)+0.2
      RETURN
C
   45 Z2(1)=TWO/DSQRT(SIX)
      NI=2*IA(1)-1
      C2(1,NI+1)=5.2
      C2(1,9)=IB(2)+0.1
      RETURN
      END
      SUBROUTINE CODED3(IA,IB,NSI,C2,NDET,Z2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL C,C2
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION C2(12,12),Z2(12),IA(2),IB(2)
C
      DATA PT5/0.5D0/,ONE/1.0D0/,TWO/2.0D0/,THREE/3.0D0/,SIX/6.0D0/
C
      NDO=NO-NSO
      NDET=NDED(NSI)
      DO 10 I=1,MNE
      DO 10 J=1,NDET
   10 C2(J,I)=C(I)
      GO TO(20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,15),
     1 NSI
   15 Z2(1)=ONE
      RETURN
C
   20 NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NI+1)=IB(2)+0.2
      Z2(1)=ONE
      RETURN
C
   25 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI)=IB(1)+0.2
      C2(2,NI-1)=IB(2)+0.1
      NJ=2*IA(2)
      C2(1,NJ-1)=IB(1)+0.1
      C2(2,NJ)=IB(2)+0.2
      RETURN
C
   30 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI-1)=IB(1)+0.1
      C2(1,NI)=IB(2)+0.2
      C2(2,NI-1)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.2
      RETURN
C
   60 Z2(1)=PT5/DSQRT(THREE)
      Z2(2)=Z2(1)
      Z2(3)=-Z2(1)
      Z2(4)=-Z2(1)
      Z2(5)=Z2(1)+Z2(1)
      Z2(6)=Z2(5)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.1
      C2(2,NJ+1)=IB(2)+0.2
      C2(3,NI+1)=IB(1)+0.2
      C2(3,NJ+1)=IB(2)+0.2
      C2(4,NI)=IB(1)+0.1
      C2(4,NJ)=IB(2)+0.1
      C2(5,NI)=IB(2)+0.1
      C2(5,NJ+1)=IB(1)+0.2
      C2(6,NI+1)=IB(2)+0.2
      C2(6,NJ)=IB(1)+0.1
      RETURN
C
   65 Z2(1)=PT5
      Z2(2)=PT5
      Z2(3)=PT5
      Z2(4)=PT5
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ+1)=IB(2)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(1)+0.1
      C2(3,NJ)=IB(2)+0.1
      C2(4,NI+1)=IB(1)+0.2
      C2(4,NJ+1)=IB(2)+0.2
      RETURN
C
   40 Z2(1)=ONE
      NI=2*IA(1)
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      RETURN
C
   80 Z2(2)=ONE/DSQRT(SIX)
      Z2(1)=Z2(2)+Z2(2)
      Z2(3)=-Z2(2)
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI+1)=IB(2)+0.2
      C2(3,NJ)=IB(1)+0.1
      RETURN
C
   85 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI+1)=IB(2)+0.2
      C2(2,NJ)=IB(1)+0.1
      RETURN
C
   55 Z2(1)=ONE
      NI=NDO+IA(1)
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ)=IB(2)+0.1
      RETURN
C
   35 Z2(1)=ONE
      NI=2*IA(1)-1
      C2(1,NI)=IB(2)+0.1
      C2(1,NI+1)=IB(1)+0.2
      RETURN
C
   70 Z2(2)=ONE/DSQRT(SIX)
      Z2(1)=Z2(2)+Z2(2)
      Z2(3)=-Z2(2)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ+1)=IB(2)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(2)+0.1
      C2(3,NJ+1)=IB(1)+0.2
      RETURN
C
   75 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(2)+0.1
      C2(2,NJ+1)=IB(1)+0.2
      RETURN
C
   45 Z2(1)=ONE
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      RETURN
C
   50 NI=2*IA(1)-1
      NJ=2*IA(2)-1
      Z2(1)=ONE
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ+1)=IB(2)+0.2
      RETURN
C
   90 Z2(1)=PT5
   26 NI=2*IA(1)-1
      C2(1,NI+1)=5.2
      C2(1,9)=IB(2)+0.1
      RETURN
C
   95 Z2(1)=-ONE/DSQRT(TWO)
      GO TO 26
  100 Z2(1)=PT5
   28 NI=2*IA(1)-1
      C2(1,NI+1)=6.2
      C2(1,10)=IB(2)+0.1
      RETURN
C
  105 Z2(1)=ONE/DSQRT(TWO)
      GO TO 28
      END
      SUBROUTINE CODED4(IA,IB,NSI,C2,NDET,Z2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL C,C2
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION C2(12,12),Z2(12),IA(2),IB(2)
C
      DATA PT5/0.5D0/,ONE/1.0D0/,TWO/2.0D0/,THREE/3.0D0/,SIX/6.0D0/,
     1 F30/30.0D0/
C
      NDO=NO-NSO
      NDET=NDED(NSI)
      DO 10 I=1,MNE
      DO 10 J=1,NDET
   10 C2(J,I)=C(I)
      GO TO(20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,
     1 110,115,120,125,15),NSI
   15 Z2(1)=ONE
      RETURN
C
   20 NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NI+1)=IB(2)+0.2
      Z2(1)=ONE
      RETURN
C
   25 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI)=IB(1)+0.2
      C2(2,NI-1)=IB(2)+0.1
      NJ=2*IA(2)
      C2(1,NJ-1)=IB(1)+0.1
      C2(2,NJ)=IB(2)+0.2
      RETURN
C
   30 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)
      C2(1,NI-1)=IB(1)+0.1
      C2(1,NI)=IB(2)+0.2
      C2(2,NI-1)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.2
      RETURN
C
   60 Z2(1)=PT5/DSQRT(THREE)
      Z2(2)=Z2(1)
      Z2(3)=-Z2(1)
      Z2(4)=-Z2(1)
      Z2(5)=Z2(1)+Z2(1)
      Z2(6)=Z2(5)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(1)+0.1
      C2(2,NJ+1)=IB(2)+0.2
      C2(3,NI+1)=IB(1)+0.2
      C2(3,NJ+1)=IB(2)+0.2
      C2(4,NI)=IB(1)+0.1
      C2(4,NJ)=IB(2)+0.1
      C2(5,NI)=IB(2)+0.1
      C2(5,NJ+1)=IB(1)+0.2
      C2(6,NI+1)=IB(2)+0.2
      C2(6,NJ)=IB(1)+0.1
      RETURN
C
   65 Z2(1)=PT5
      Z2(2)=PT5
      Z2(3)=PT5
      Z2(4)=PT5
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ+1)=IB(2)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(1)+0.1
      C2(3,NJ)=IB(2)+0.1
      C2(4,NI+1)=IB(1)+0.2
      C2(4,NJ+1)=IB(2)+0.2
      RETURN
C
   40 Z2(1)=ONE
      NI=2*IA(1)
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      RETURN
C
   80 Z2(2)=ONE/DSQRT(SIX)
      Z2(1)=Z2(2)+Z2(2)
      Z2(3)=-Z2(2)
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI+1)=IB(2)+0.2
      C2(3,NJ)=IB(1)+0.1
      RETURN
C
   85 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI+1)=IB(2)+0.2
      C2(2,NJ)=IB(1)+0.1
      RETURN
C
   35 Z2(1)=ONE
      NI=2*IA(1)-1
      C2(1,NI)=IB(2)+0.1
      C2(1,NI+1)=IB(1)+0.2
      RETURN
C
   70 Z2(2)=ONE/DSQRT(SIX)
      Z2(1)=Z2(2)+Z2(2)
      Z2(3)=-Z2(2)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ+1)=IB(2)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(2)+0.1
      C2(3,NJ+1)=IB(1)+0.2
      RETURN
C
   75 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(2)+0.1
      C2(2,NJ+1)=IB(1)+0.2
      RETURN
C
   45 Z2(1)=ONE
      NI=2*IA(1)-1
      NJ=NDO+IA(2)
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ)=IB(2)+0.1
      RETURN
C
   50 NI=2*IA(1)-1
      NJ=2*IA(2)-1
      Z2(1)=ONE
      C2(1,NI+1)=IB(1)+0.2
      C2(1,NJ+1)=IB(2)+0.2
      RETURN
C
   55 Z2(1)=ONE
      NI=NDO+IA(1)
      NJ=NDO+IA(2)
      C2(1,NI)=IB(1)+0.1
      C2(1,NJ)=IB(2)+0.1
      RETURN
C
   90 Z2(1)=TWO/DSQRT(F30)
      NI=2*IA(1)-1
      C2(1,NI+1)=5.2
      C2(1,9)=IB(2)+0.1
      RETURN
C
   95 Z2(1)=TWO/DSQRT(SIX)
      NI=2*IA(1)-1
      C2(1,NI+1)=5.2
      C2(1,9)=IB(2)+0.1
      RETURN
C
  100 Z2(1)=TWO/DSQRT(F30)
      NI=2*IA(1)-1
      C2(1,NI+1)=6.2
      C2(1,10)=IB(2)+0.1
      RETURN
C
  105 Z2(1)=-ONE/DSQRT(SIX)
      NI=2*IA(1)-1
      C2(1,NI+1)=6.2
      C2(1,10)=IB(2)+0.1
      RETURN
C
  110 Z2(1)=ONE/DSQRT(TWO)
      NI=2*IA(1)-1
      C2(1,NI+1)=6.2
      C2(1,10)=IB(2)+0.1
      RETURN
C
  115 Z2(1)=TWO/DSQRT(F30)
      NI=2*IA(1)-1
      C2(1,NI+1)=7.2
      C2(1,11)=IB(2)+0.1
      RETURN
C
  120 Z2(1)=-ONE/DSQRT(SIX)
      NI=2*IA(1)-1
      C2(1,NI+1)=7.2
      C2(1,11)=IB(2)+0.1
      RETURN
C
  125 Z2(1)=-ONE/DSQRT(TWO)
      NI=2*IA(1)-1
      C2(1,NI+1)=7.2
      C2(1,11)=IB(2)+0.1
      RETURN
      END
      SUBROUTINE CODED5(IA,IB,NSI,C2,NDET,Z2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL C,C2
C
      COMMON/C902C/NSC,C(12),NDED(26),NO,NSO,MNE,MNX1,MNX2,
     1 NDES(6),NDC,NDER(4),ICASE,NCT,MCT(10),NCTS,MCTS(10)
     2 ,IFILL(30),FILL(128),J7
C
      DIMENSION C2(12,12),Z2(12),IA(2),IB(2)
C
      DATA PT5/0.5D0/,ONE/1.0D0/,TWO/2.0D0/,THREE/3.0D0/,EIGHT/8.0D0/,
     1 TWELVE/12.0D0/,TWENT4/24.0D0/
C
      NDET=NDED(NSI)
      NDET1=NDET/2
      DO 100 J=1,NDET1
      DO 100 I=1,MNE
  100 C2(J,I)=C(I)
      DO 102 J=1,NDET1
      DO 101 I=1,MNE
  101 C2(J+NDET1,I)=C(I)
      C2(J+NDET1,9)=6.1
  102 C2(J+NDET1,10)=5.2
      GO TO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     1 23,24)NSI
   24 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      RETURN
C
    1 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,NI+1)=IB(2)+0.2
      C2(2,NI)=C2(1,NI)
      C2(2,NI+1)=C2(1,NI+1)
      RETURN
C
    2 K=-2
  200 K=K+2
      IF(K.EQ.4)RETURN
      Z2(1+K)=PT5
      Z2(2+K)=Z2(1+K)
      NI=2*IA(1)
      C2(1+K,NI)=IB(1)+0.2
      C2(2+K,NI-1)=IB(2)+0.1
      NJ=2*IA(2)
      C2(1+K,NJ-1)=IB(1)+0.1
      C2(2+K,NJ)=IB(2)+0.2
      GO TO 200
C
    3 K=-2
   30 K=K+2
      IF(K.EQ.4)RETURN
      Z2(1+K)=PT5
      Z2(2+K)=Z2(1+K)
      NI=2*IA(1)
      C2(K+1,NI-1)=IB(1)+0.1
      C2(K+1,NI)=IB(2)+0.2
      C2(K+2,NI-1)=IB(2)+0.1
      C2(2+K,NI)=IB(1)+0.2
      GO TO 30
C
   20 K=-6
   40 K=K+6
      IF(K.EQ.12)RETURN
      Z2(1+K)=ONE/DSQRT(TWENT4)
      Z2(2+K)=Z2(1+K)
      Z2(3+K)=-Z2(1+K)
      Z2(4+K)=-Z2(1+K)
      Z2(5+K)=TWO*Z2(1+K)
      Z2(6+K)=Z2(5+K)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1+K,NI+1)=IB(1)+0.2
      C2(3+K,NI+1)=IB(1)+0.2
      C2(1+K,NJ)=IB(2)+0.1
      C2(3+K,NJ+1)=IB(2)+0.2
      C2(2+K,NI)=IB(1)+0.1
      C2(4+K,NI)=IB(1)+0.1
      C2(2+K,NJ+1)=IB(2)+0.2
      C2(4+K,NJ)=IB(2)+0.1
      C2(5+K,NI)=IB(2)+0.1
      C2(5+K,NJ+1)=IB(1)+0.2
      C2(6+K,NI+1)=IB(2)+0.2
      C2(6+K,NJ)=IB(1)+0.1
      GO TO 40
C
   21 K=-4
   50 K=K+4
      IF(K.EQ.8)RETURN
      Z2(1+K)=ONE/DSQRT(EIGHT)
      Z2(2+K)=Z2(1+K)
      Z2(3+K)=Z2(1+K)
      Z2(4+K)=Z2(3+K)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1+K,NI)=IB(1)+0.1
      C2(1+K,NJ+1)=IB(2)+0.2
      C2(2+K,NI+1)=IB(1)+0.2
      C2(2+K,NJ)=IB(2)+0.1
      C2(3+K,NI)=IB(1)+0.1
      C2(3+K,NJ)=IB(2)+0.1
      C2(4+K,NI+1)=IB(1)+0.2
      C2(4+K,NJ+1)=IB(2)+0.2
      GO TO 50
C
    6 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(2,10)=IB(1)+0.2
      C2(1,9)=IB(1)+0.1
      C2(2,NI)=IB(1)+0.1
      RETURN
C
    7 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,10)=IB(1)+0.2
      C2(2,NI+1)=IB(1)+0.2
      C2(2,9)=IB(1)+0.1
      RETURN
C
    8 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=5.2
      C2(1,NJ)=6.1
      C2(2,NI)=5.1
      C2(2,NJ+1)=6.2
      RETURN
C
    9 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      C2(1,9)=IB(1)+0.1
      C2(1,10)=IB(2)+0.2
      C2(2,9)=IB(2)+0.1
      C2(2,10)=IB(1)+0.2
      RETURN
C
   10 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      C2(1,NI)=IB(2)+0.1
      C2(1,9)=6.1
      C2(2,NI+1)=IB(2)+0.2
      C2(2,10)=6.2
      RETURN
C
   11 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      C2(1,NI+1)=IB(2)+0.2
      C2(1,10)=5.2
      C2(2,NI)=IB(2)+0.1
      C2(2,9)=5.1
      RETURN
C
   12 Z2(2)=-ONE/DSQRT(TWELVE)
      Z2(3)=-Z2(2)
      Z2(1)=TWO*Z2(2)
      Z2(4)=Z2(1)
      Z2(5)=Z2(2)
      Z2(6)=Z2(3)
      NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,9)=IB(2)+0.1
      C2(2,NI+1)=IB(1)+0.2
      C2(2,9)=IB(2)+0.1
      C2(3,NI+1)=IB(2)+0.2
      C2(3,9)=IB(1)+0.1
      C2(4,NI+1)=IB(1)+0.2
      C2(4,10)=IB(2)+0.2
      C2(5,NI)=IB(1)+0.1
      C2(5,10)=IB(2)+0.2
      C2(6,NI)=IB(2)+0.1
      C2(6,10)=IB(1)+0.2
      RETURN
C
   13 Z2(1)=PT5
      Z2(2)=Z2(1)
      Z2(3)=Z2(2)
      Z2(4)=Z2(3)
      NI=2*IA(1)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,9)=IB(2)+0.1
      C2(2,NI+1)=IB(2)+0.2
      C2(2,9)=IB(1)+0.1
      C2(3,NI)=IB(1)+0.1
      C2(3,10)=IB(2)+0.2
      C2(4,NI)=IB(2)+0.1
      C2(4,10)=IB(1)+0.2
      RETURN
C
   14 Z2(2)=-ONE/DSQRT(TWELVE)
      Z2(3)=-Z2(2)
      Z2(1)=TWO*Z2(2)
      Z2(4)=Z2(1)
      Z2(5)=Z2(2)
      Z2(6)=Z2(3)
      NI=2*IA(1)-1
      C2(1,NI+1)=IB(1)+0.2
      C2(1,10)=IB(2)+0.2
      C2(2,NI)=IB(1)+0.1
      C2(2,10)=IB(2)+0.2
      C2(3,NI)=IB(2)+0.1
      C2(3,10)=IB(1)+0.2
      C2(4,NI)=IB(1)+0.1
      C2(4,9)=IB(2)+0.1
      C2(5,NI+1)=IB(1)+0.2
      C2(5,9)=IB(2)+0.1
      C2(6,NI+1)=IB(2)+0.2
      C2(6,9)=IB(1)+0.1
      RETURN
C
   15 Z2(1)=PT5
      Z2(2)=Z2(1)
      Z2(3)=Z2(2)
      Z2(4)=Z2(3)
      NI=2*IA(1)-1
      C2(1,NI)=IB(1)+0.1
      C2(1,10)=IB(2)+0.2
      C2(2,NI)=IB(2)+0.1
      C2(2,10)=IB(1)+0.2
      C2(3,NI+1)=IB(1)+0.2
      C2(3,9)=IB(2)+0.1
      C2(4,NI+1)=IB(2)+0.2
      C2(4,9)=IB(1)+0.1
      RETURN
C
   16 Z2(2)=ONE/DSQRT(TWELVE)
      Z2(1)=TWO*Z2(2)
      Z2(3)=-Z2(2)
      Z2(4)=Z2(1)
      Z2(5)=-Z2(2)
      Z2(6)=-Z2(3)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=5.2
      C2(1,NJ+1)=IB(2)+0.2
      C2(4,NI)=5.1
      C2(4,NJ)=IB(2)+0.1
      C2(2,NI+1)=5.2
      C2(2,NJ)=IB(2)+0.1
      C2(3,NI)=IB(2)+0.1
      C2(3,NJ+1)=5.2
      C2(5,NI+1)=IB(2)+0.2
      C2(5,NJ)=5.1
      C2(6,NJ+1)=IB(2)+0.2
      C2(6,NI)=5.1
      RETURN
C
   17 Z2(1)=PT5
      Z2(2)=Z2(1)
      Z2(3)=Z2(1)
      Z2(4)=Z2(1)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI+1)=5.2
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=IB(2)+0.1
      C2(2,NJ+1)=5.2
      C2(3,NI+1)=IB(2)+0.2
      C2(3,NJ)=5.1
      C2(4,NJ+1)=IB(2)+0.2
      C2(4,NI)=5.1
      RETURN
C
   18 Z2(2)=ONE/DSQRT(TWELVE)
      Z2(3)=-Z2(2)
      Z2(1)=TWO*Z2(2)
      Z2(4)=Z2(1)
      Z2(5)=-Z2(2)
      Z2(6)=-Z2(3)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI)=6.1
      C2(1,NJ)=IB(2)+0.1
      C2(2,NI)=6.1
      C2(2,NJ+1)=IB(2)+0.2
      C2(3,NI+1)=IB(2)+0.2
      C2(3,NJ)=6.1
      C2(4,NJ+1)=IB(2)+0.2
      C2(4,NI+1)=6.2
      C2(5,NI)=IB(2)+0.1
      C2(5,NJ+1)=6.2
      C2(6,NJ)=IB(2)+0.1
      C2(6,NI+1)=6.2
      RETURN
C
   19 Z2(1)=PT5
      Z2(2)=Z2(1)
      Z2(3)=Z2(1)
      Z2(4)=Z2(1)
      NI=2*IA(1)-1
      NJ=2*IA(2)-1
      C2(1,NI)=6.1
      C2(1,NJ+1)=IB(2)+0.2
      C2(3,NI)=IB(2)+0.1
      C2(3,NJ+1)=6.2
      C2(2,NI+1)=IB(2)+0.2
      C2(2,NJ)=6.1
      C2(4,NJ)=IB(2)+0.1
      C2(4,NI+1)=6.2
      RETURN
C
    4 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      C2(1,NI)=IB(2)+0.1
      C2(1,NI+1)=5.2
      C2(2,NI+1)=IB(2)+0.2
      C2(2,NI)=5.1
      RETURN
C
    5 Z2(1)=ONE/DSQRT(TWO)
      Z2(2)=Z2(1)
      NI=2*IA(1)-1
      C2(1,NI+1)=IB(2)+0.2
      C2(1,NI)=6.1
      C2(2,NI)=IB(2)+0.1
      C2(2,NI+1)=6.2
      RETURN
C
   22 Z2(1)=-TWO/DSQRT(TWELVE)
      NI=2*IA(1)-1
      C2(1,NI+1)=5.2
      C2(1,9)=IB(2)+0.1
      Z2(2)=Z2(1)
      C2(2,NI)=5.1
      C2(2,10)=IB(2)+0.2
      RETURN
C
   23 Z2(1)=TWO/DSQRT(TWELVE)
      NI=2*IA(1)-1
      C2(1,NI)=6.1
      C2(1,10)=IB(2)+0.2
      Z2(2)=Z2(1)
      C2(2,NI+1)=6.2
      C2(2,9)=IB(2)+0.1
      RETURN
      END
