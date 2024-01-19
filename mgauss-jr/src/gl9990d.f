C/    GL9990D      24 FEB 87                                         MRP
      SUBROUTINE F02CTR(NZERO)
C***********************************************************************
C     ROUTINE TO COMPUTE TWO-DIMENSIONAL 2-CENTER INTEGRALS
C     GIVEN THE G'S AND A'S.
C     THIS ROUTINE FILLS UP TO NZERO*LPMAX*LQMAX LOCATIONS IN EACH
C     OF GX, GY AND GZ.
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/FST2D0/GX(45),GY(45),GZ(45),
     $              X2(125),Y2(125),Z2(125),
     $              X3(225),Y3(225),Z3(225),
     $              X4(405),Y4(405),Z4(405)
C
C     (LPMAX,LQMAX) INDEPENDENT INITIALIZATION.
      N2=LPMAX*LQMAX
      I2=-N2
      IG=-LPQMAX
C
C     BRANCH TO APPROPRIATE CODE.
      GO TO(10,20,30,40,50),LPMAX
   10 GO TO(100,101,102,103,104),LQMAX
   20 GO TO(110,111,112,113,114),LQMAX
   30 GO TO(120,121,122,123,124),LQMAX
   40 GO TO(130,131,132,133,134),LQMAX
   50 GO TO(140,141,142,143,144),LQMAX
C
C.....  CASE (0,0).
  100 DO 210 I=1,NZERO
      X2(I)=GX(I)
      Y2(I)=GY(I)
  210 Z2(I)=GZ(I)
      GO TO 500
C
C.....  CASE (0,1)
CDIR$ IVDEP
101   DO 220 I=1,NZERO
      I2=I2+N2
      IG=IG+LPQMAX
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Z2(1+I2)=      GZ(1+IG)
  220 Z2(2+I2)=A( 1)*GZ(2+IG)
      GO TO 500
C
C.....  CASE (0,2)
CDIR$ IVDEP
102   I2=-3
      IG=-3
      DO 230 I=1,NZERO
      I2=I2+3
      IG=IG+3
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
  230 Z2(3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      GO TO 500
C
C.....  CASE (0,3)
CDIR$ IVDEP
103   I2=-4
      IG=-4
      DO 240 I=1,NZERO
      I2=I2+4
      IG=IG+4
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2(4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2(4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
  240 Z2(4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      GO TO 500
C
C.....  CASE (0,4)
CDIR$ IVDEP
104   I2=-5
      IG=-5
      DO 250 I=1,NZERO
      I2=I2+5
      IG=IG+5
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2(4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2(5+I2)=A( 6)*GX(1+IG)+A( 7)*GX(3+IG)+A( 8)*GX(5+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2(4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2(5+I2)=A( 6)*GY(1+IG)+A( 7)*GY(3+IG)+A( 8)*GY(5+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2(4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
  250 Z2(5+I2)=A( 6)*GZ(1+IG)+A( 7)*GZ(3+IG)+A( 8)*GZ(5+IG)
      GO TO 500
C
C.....  CASE (1,0)
CDIR$ IVDEP
110   I2=-2
      IG=-2
      DO 260 I=1,NZERO
      I2=I2+2
      IG=IG+2
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A(16)*GX(2+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A(16)*GY(2+IG)
      Z2(1+I2)=      GZ(1+IG)
  260 Z2(2+I2)=A(16)*GZ(2+IG)
      GO TO 500
C
C.....  CASE (1,1).
CDIR$ IVDEP
111   I2=-4
      IG=-3
      DO 270 I=1,NZERO
      I2=I2+4
      IG=IG+3
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A(16)*GX(2+IG)
      X2(4+I2)=A(17)*GX(3+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A(16)*GY(2+IG)
      Y2(4+I2)=A(17)*GY(3+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A(16)*GZ(2+IG)
  270 Z2(4+I2)=A(17)*GZ(3+IG)
      GO TO 500
C
C.....  CASE (1,2).
CDIR$ IVDEP
112   I2=-6
      IG=-4
      DO 280 I=1,NZERO
      I2=I2+6
      IG=IG+4
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2(4+I2)=A(16)*GX(2+IG)
      X2(5+I2)=A(17)*GX(3+IG)
      X2(6+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2(4+I2)=A(16)*GY(2+IG)
      Y2(5+I2)=A(17)*GY(3+IG)
      Y2(6+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2(4+I2)=A(16)*GZ(2+IG)
      Z2(5+I2)=A(17)*GZ(3+IG)
  280 Z2(6+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      GO TO 500
C
C.....  CASE (1,3).
CDIR$ IVDEP
113   I2=-8
      IG=-5
      DO 290 I=1,NZERO
      I2=I2+8
      IG=IG+5
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2(4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2(5+I2)=A(16)*GX(2+IG)
      X2(6+I2)=A(17)*GX(3+IG)
      X2(7+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2(8+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2(4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2(5+I2)=A(16)*GY(2+IG)
      Y2(6+I2)=A(17)*GY(3+IG)
      Y2(7+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2(8+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2(4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2(5+I2)=A(16)*GZ(2+IG)
      Z2(6+I2)=A(17)*GZ(3+IG)
      Z2(7+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
  290 Z2(8+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
      GO TO 500
C
C.....  CASE (1,4).
CDIR$ IVDEP
114   I2=-10
      IG=-6
      DO 300 I=1,NZERO
      I2=I2+10
      IG=IG+6
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2( 5+I2)=A( 6)*GX(1+IG)+A( 7)*GX(3+IG)+A( 8)*GX(5+IG)
      X2( 6+I2)=A(16)*GX(2+IG)
      X2( 7+I2)=A(17)*GX(3+IG)
      X2( 8+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 9+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      X2(10+I2)=A(22)*GX(2+IG)+A(23)*GX(4+IG)+A(24)*GX(6+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2( 5+I2)=A( 6)*GY(1+IG)+A( 7)*GY(3+IG)+A( 8)*GY(5+IG)
      Y2( 6+I2)=A(16)*GY(2+IG)
      Y2( 7+I2)=A(17)*GY(3+IG)
      Y2( 8+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 9+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Y2(10+I2)=A(22)*GY(2+IG)+A(23)*GY(4+IG)+A(24)*GY(6+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2( 5+I2)=A( 6)*GZ(1+IG)+A( 7)*GZ(3+IG)+A( 8)*GZ(5+IG)
      Z2( 6+I2)=A(16)*GZ(2+IG)
      Z2( 7+I2)=A(17)*GZ(3+IG)
      Z2( 8+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 9+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
  300 Z2(10+I2)=A(22)*GZ(2+IG)+A(23)*GZ(4+IG)+A(24)*GZ(6+IG)
      GO TO 500
C
C.....  CASE (2,0)
CDIR$ IVDEP
120   I2=-3
      IG=-3
      DO 310 I=1,NZERO
      I2=I2+3
      IG=IG+3
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A(16)*GX(2+IG)
      X2(3+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A(16)*GY(2+IG)
      Y2(3+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A(16)*GZ(2+IG)
  310 Z2(3+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      GO TO 500
C
C.....  CASE (2,1).
CDIR$ IVDEP
121   I2=-6
      IG=-4
      DO 320 I=1,NZERO
      I2=I2+6
      IG=IG+4
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A(16)*GX(2+IG)
      X2(4+I2)=A(17)*GX(3+IG)
      X2(5+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(6+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A(16)*GY(2+IG)
      Y2(4+I2)=A(17)*GY(3+IG)
      Y2(5+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(6+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A(16)*GZ(2+IG)
      Z2(4+I2)=A(17)*GZ(3+IG)
      Z2(5+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
  320 Z2(6+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      GO TO 500
C
C.....  CASE (2,2).
CDIR$ IVDEP
122   I2=-9
      IG=-5
      DO 330 I=1,NZERO
      I2=I2+9
      IG=IG+5
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2(4+I2)=A(16)*GX(2+IG)
      X2(5+I2)=A(17)*GX(3+IG)
      X2(6+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2(7+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(8+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(9+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2(4+I2)=A(16)*GY(2+IG)
      Y2(5+I2)=A(17)*GY(3+IG)
      Y2(6+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2(7+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(8+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(9+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2(4+I2)=A(16)*GZ(2+IG)
      Z2(5+I2)=A(17)*GZ(3+IG)
      Z2(6+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2(7+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(8+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
  330 Z2(9+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      GO TO 500
C
C.....  CASE (2,3).
CDIR$ IVDEP
123   I2=-12
      IG=-6
      DO 340 I=1,NZERO
      I2=I2+12
      IG=IG+6
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2( 5+I2)=A(16)*GX(2+IG)
      X2( 6+I2)=A(17)*GX(3+IG)
      X2( 7+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 8+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      X2( 9+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(10+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(11+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(12+I2)=A(39)*GX(2+IG)+A(40)*GX(4+IG)+A(41)*GX(6+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2( 5+I2)=A(16)*GY(2+IG)
      Y2( 6+I2)=A(17)*GY(3+IG)
      Y2( 7+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 8+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Y2( 9+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(10+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(11+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(12+I2)=A(39)*GY(2+IG)+A(40)*GY(4+IG)+A(41)*GY(6+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2( 5+I2)=A(16)*GZ(2+IG)
      Z2( 6+I2)=A(17)*GZ(3+IG)
      Z2( 7+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 8+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
      Z2( 9+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(10+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2(11+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
  340 Z2(12+I2)=A(39)*GZ(2+IG)+A(40)*GZ(4+IG)+A(41)*GZ(6+IG)
      GO TO 500
C
C.....  CASE (2,4).
CDIR$ IVDEP
124   I2=-15
      IG=-7
      DO 350 I=1,NZERO
      I2=I2+15
      IG=IG+7
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2( 5+I2)=A( 6)*GX(1+IG)+A( 7)*GX(3+IG)+A( 8)*GX(5+IG)
      X2( 6+I2)=A(16)*GX(2+IG)
      X2( 7+I2)=A(17)*GX(3+IG)
      X2( 8+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 9+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      X2(10+I2)=A(22)*GX(2+IG)+A(23)*GX(4+IG)+A(24)*GX(6+IG)
      X2(11+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(12+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(13+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(14+I2)=A(39)*GX(2+IG)+A(40)*GX(4+IG)+A(41)*GX(6+IG)
      X2(15+I2)=A(42)*GX(1+IG)+A(43)*GX(3+IG)+A(44)*GX(5+IG)+
     $          A(45)*GX(7+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2( 5+I2)=A( 6)*GY(1+IG)+A( 7)*GY(3+IG)+A( 8)*GY(5+IG)
      Y2( 6+I2)=A(16)*GY(2+IG)
      Y2( 7+I2)=A(17)*GY(3+IG)
      Y2( 8+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 9+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Y2(10+I2)=A(22)*GY(2+IG)+A(23)*GY(4+IG)+A(24)*GY(6+IG)
      Y2(11+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(12+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(13+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(14+I2)=A(39)*GY(2+IG)+A(40)*GY(4+IG)+A(41)*GY(6+IG)
      Y2(15+I2)=A(42)*GY(1+IG)+A(43)*GY(3+IG)+A(44)*GY(5+IG)+
     $          A(45)*GY(7+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2( 5+I2)=A( 6)*GZ(1+IG)+A( 7)*GZ(3+IG)+A( 8)*GZ(5+IG)
      Z2( 6+I2)=A(16)*GZ(2+IG)
      Z2( 7+I2)=A(17)*GZ(3+IG)
      Z2( 8+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 9+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
      Z2(10+I2)=A(22)*GZ(2+IG)+A(23)*GZ(4+IG)+A(24)*GZ(6+IG)
      Z2(11+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(12+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2(13+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      Z2(14+I2)=A(39)*GZ(2+IG)+A(40)*GZ(4+IG)+A(41)*GZ(6+IG)
  350 Z2(15+I2)=A(42)*GZ(1+IG)+A(43)*GZ(3+IG)+A(44)*GZ(5+IG)+
     $          A(45)*GZ(7+IG)
      GO TO 500
C
C.....  CASE (3,0)
CDIR$ IVDEP
130   I2=-4
      IG=-4
      DO 360 I=1,NZERO
      I2=I2+4
      IG=IG+4
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A(16)*GX(2+IG)
      X2(3+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(4+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A(16)*GY(2+IG)
      Y2(3+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(4+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A(16)*GZ(2+IG)
      Z2(3+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
  360 Z2(4+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      GO TO 500
C
C.....  CASE (3,1)
CDIR$ IVDEP
131   I2=-8
      IG=-5
      DO 370 I=1,NZERO
      I2=I2+8
      IG=IG+5
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A( 1)*GX(2+IG)
      X2(3+I2)=A(16)*GX(2+IG)
      X2(4+I2)=A(17)*GX(3+IG)
      X2(5+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(6+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(7+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(8+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A( 1)*GY(2+IG)
      Y2(3+I2)=A(16)*GY(2+IG)
      Y2(4+I2)=A(17)*GY(3+IG)
      Y2(5+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(6+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(7+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(8+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A( 1)*GZ(2+IG)
      Z2(3+I2)=A(16)*GZ(2+IG)
      Z2(4+I2)=A(17)*GZ(3+IG)
      Z2(5+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(6+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2(7+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
  370 Z2(8+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
      GO TO 500
C
C.....  CASE (3,2)
CDIR$ IVDEP
132   I2=-12
      IG=-6
      DO 380 I=1,NZERO
      I2=I2+12
      IG=IG+6
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A(16)*GX(2+IG)
      X2( 5+I2)=A(17)*GX(3+IG)
      X2( 6+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 7+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2( 8+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2( 9+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(10+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(11+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      X2(12+I2)=A(59)*GX(2+IG)+A(60)*GX(4+IG)+A(61)*GX(6+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A(16)*GY(2+IG)
      Y2( 5+I2)=A(17)*GY(3+IG)
      Y2( 6+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 7+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2( 8+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2( 9+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(10+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(11+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Y2(12+I2)=A(59)*GY(2+IG)+A(60)*GY(4+IG)+A(61)*GY(6+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A(16)*GZ(2+IG)
      Z2( 5+I2)=A(17)*GZ(3+IG)
      Z2( 6+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 7+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2( 8+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2( 9+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      Z2(10+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      Z2(11+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
  380 Z2(12+I2)=A(59)*GZ(2+IG)+A(60)*GZ(4+IG)+A(61)*GZ(6+IG)
      GO TO 500
C
C.....  CASE (3,3)
CDIR$ IVDEP
133   I2=-16
      IG=-7
      DO 390 I=1,NZERO
      I2=I2+16
      IG=IG+7
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2( 5+I2)=A(16)*GX(2+IG)
      X2( 6+I2)=A(17)*GX(3+IG)
      X2( 7+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 8+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      X2( 9+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(10+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(11+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(12+I2)=A(39)*GX(2+IG)+A(40)*GX(4+IG)+A(41)*GX(6+IG)
      X2(13+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(14+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      X2(15+I2)=A(59)*GX(2+IG)+A(60)*GX(4+IG)+A(61)*GX(6+IG)
      X2(16+I2)=A(62)*GX(3+IG)+A(63)*GX(5+IG)+A(64)*GX(7+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2( 5+I2)=A(16)*GY(2+IG)
      Y2( 6+I2)=A(17)*GY(3+IG)
      Y2( 7+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 8+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Y2( 9+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(10+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(11+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(12+I2)=A(39)*GY(2+IG)+A(40)*GY(4+IG)+A(41)*GY(6+IG)
      Y2(13+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(14+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Y2(15+I2)=A(59)*GY(2+IG)+A(60)*GY(4+IG)+A(61)*GY(6+IG)
      Y2(16+I2)=A(62)*GY(3+IG)+A(63)*GY(5+IG)+A(64)*GY(7+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2( 5+I2)=A(16)*GZ(2+IG)
      Z2( 6+I2)=A(17)*GZ(3+IG)
      Z2( 7+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 8+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
      Z2( 9+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(10+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2(11+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      Z2(12+I2)=A(39)*GZ(2+IG)+A(40)*GZ(4+IG)+A(41)*GZ(6+IG)
      Z2(13+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      Z2(14+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
      Z2(15+I2)=A(59)*GZ(2+IG)+A(60)*GZ(4+IG)+A(61)*GZ(6+IG)
  390 Z2(16+I2)=A(62)*GZ(3+IG)+A(63)*GZ(5+IG)+A(64)*GZ(7+IG)
      GO TO 500
C
C.....  CASE (3,4)
CDIR$ IVDEP
134   I2=-20
      IG=-8
      DO 400 I=1,NZERO
      I2=I2+20
      IG=IG+8
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2( 5+I2)=A(6)*GX(1+IG)+A( 7)*GX(3+IG)+A( 8)*GX(5+IG)
      X2( 6+I2)=A(16)*GX(2+IG)
      X2( 7+I2)=A(17)*GX(3+IG)
      X2( 8+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 9+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      X2(10+I2)=A(22)*GX(2+IG)+A(23)*GX(4+IG)+A(24)*GX(6+IG)
      X2(11+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(12+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(13+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(14+I2)=A(39)*GX(2+IG)+A(40)*GX(4+IG)+A(41)*GX(6+IG)
      X2(15+I2)=A(42)*GX(1+IG)+A(43)*GX(3+IG)+A(44)*GX(5+IG)+
     $          A(45)*GX(7+IG)
      X2(16+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(17+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      X2(18+I2)=A(59)*GX(2+IG)+A(60)*GX(4+IG)+A(61)*GX(6+IG)
      X2(19+I2)=A(62)*GX(3+IG)+A(63)*GX(5+IG)+A(64)*GX(7+IG)
      X2(20+I2)=A(65)*GX(2+IG)+A(66)*GX(4+IG)+A(67)*GX(6+IG)+
     $          A(68)*GX(8+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2( 5+I2)=A(6)*GY(1+IG)+A( 7)*GY(3+IG)+A( 8)*GY(5+IG)
      Y2( 6+I2)=A(16)*GY(2+IG)
      Y2( 7+I2)=A(17)*GY(3+IG)
      Y2( 8+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 9+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Y2(10+I2)=A(22)*GY(2+IG)+A(23)*GY(4+IG)+A(24)*GY(6+IG)
      Y2(11+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(12+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(13+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(14+I2)=A(39)*GY(2+IG)+A(40)*GY(4+IG)+A(41)*GY(6+IG)
      Y2(15+I2)=A(42)*GY(1+IG)+A(43)*GY(3+IG)+A(44)*GY(5+IG)+
     $          A(45)*GY(7+IG)
      Y2(16+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(17+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Y2(18+I2)=A(59)*GY(2+IG)+A(60)*GY(4+IG)+A(61)*GY(6+IG)
      Y2(19+I2)=A(62)*GY(3+IG)+A(63)*GY(5+IG)+A(64)*GY(7+IG)
      Y2(20+I2)=A(65)*GY(2+IG)+A(66)*GY(4+IG)+A(67)*GY(6+IG)+
     $          A(68)*GY(8+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2( 5+I2)=A(6)*GZ(1+IG)+A( 7)*GZ(3+IG)+A( 8)*GZ(5+IG)
      Z2( 6+I2)=A(16)*GZ(2+IG)
      Z2( 7+I2)=A(17)*GZ(3+IG)
      Z2( 8+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 9+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
      Z2(10+I2)=A(22)*GZ(2+IG)+A(23)*GZ(4+IG)+A(24)*GZ(6+IG)
      Z2(11+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(12+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2(13+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      Z2(14+I2)=A(39)*GZ(2+IG)+A(40)*GZ(4+IG)+A(41)*GZ(6+IG)
      Z2(15+I2)=A(42)*GZ(1+IG)+A(43)*GZ(3+IG)+A(44)*GZ(5+IG)+
     $          A(45)*GZ(7+IG)
      Z2(16+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      Z2(17+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
      Z2(18+I2)=A(59)*GZ(2+IG)+A(60)*GZ(4+IG)+A(61)*GZ(6+IG)
      Z2(19+I2)=A(62)*GZ(3+IG)+A(63)*GZ(5+IG)+A(64)*GZ(7+IG)
  400 Z2(20+I2)=A(65)*GZ(2+IG)+A(66)*GZ(4+IG)+A(67)*GZ(6+IG)+
     $          A(68)*GZ(8+IG)
      GO TO 500
C
C.....  CASE (4,0)
CDIR$ IVDEP
140   I2=-5
      IG=-5
      DO 410 I=1,NZERO
      I2=I2+5
      IG=IG+5
      X2(1+I2)=      GX(1+IG)
      X2(2+I2)=A(16)*GX(2+IG)
      X2(3+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(4+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(5+I2)=A(78)*GX(1+IG)+A(79)*GX(3+IG)+A(80)*GX(5+IG)
      Y2(1+I2)=      GY(1+IG)
      Y2(2+I2)=A(16)*GY(2+IG)
      Y2(3+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(4+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(5+I2)=A(78)*GY(1+IG)+A(79)*GY(3+IG)+A(80)*GY(5+IG)
      Z2(1+I2)=      GZ(1+IG)
      Z2(2+I2)=A(16)*GZ(2+IG)
      Z2(3+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(4+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
  410 Z2(5+I2)=A(78)*GZ(1+IG)+A(79)*GZ(3+IG)+A(80)*GZ(5+IG)
      GO TO 500
C
C.....  CASE (4,1)
CDIR$ IVDEP
141   I2=-10
      IG=-6
      DO 420 I=1,NZERO
      I2=I2+10
      IG=IG+6
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A(16)*GX(2+IG)
      X2( 4+I2)=A(17)*GX(3+IG)
      X2( 5+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2( 6+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2( 7+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2( 8+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      X2( 9+I2)=A(78)*GX(1+IG)+A(79)*GX(3+IG)+A(80)*GX(5+IG)
      X2(10+I2)=A(81)*GX(2+IG)+A(82)*GX(4+IG)+A(83)*GX(6+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A(16)*GY(2+IG)
      Y2( 4+I2)=A(17)*GY(3+IG)
      Y2( 5+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2( 6+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2( 7+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2( 8+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Y2( 9+I2)=A(78)*GY(1+IG)+A(79)*GY(3+IG)+A(80)*GY(5+IG)
      Y2(10+I2)=A(81)*GY(2+IG)+A(82)*GY(4+IG)+A(83)*GY(6+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A(16)*GZ(2+IG)
      Z2( 4+I2)=A(17)*GZ(3+IG)
      Z2( 5+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2( 6+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2( 7+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      Z2( 8+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
      Z2( 9+I2)=A(78)*GZ(1+IG)+A(79)*GZ(3+IG)+A(80)*GZ(5+IG)
  420 Z2(10+I2)=A(81)*GZ(2+IG)+A(82)*GZ(4+IG)+A(83)*GZ(6+IG)
      GO TO 500
C
C.....  CASE (4,2)
CDIR$ IVDEP
142   I2=-15
      IG=-7
      DO 430 I=1,NZERO
      I2=I2+15
      IG=IG+7
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A(16)*GX(2+IG)
      X2( 5+I2)=A(17)*GX(3+IG)
      X2( 6+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 7+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2( 8+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2( 9+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(10+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(11+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      X2(12+I2)=A(59)*GX(2+IG)+A(60)*GX(4+IG)+A(61)*GX(6+IG)
      X2(13+I2)=A(78)*GX(1+IG)+A(79)*GX(3+IG)+A(80)*GX(5+IG)
      X2(14+I2)=A(81)*GX(2+IG)+A(82)*GX(4+IG)+A(83)*GX(6+IG)
      X2(15+I2)=A(84)*GX(1+IG)+A(85)*GX(3+IG)+A(86)*GX(5+IG)+
     $          A(87)*GX(7+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A(16)*GY(2+IG)
      Y2( 5+I2)=A(17)*GY(3+IG)
      Y2( 6+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 7+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2( 8+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2( 9+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(10+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(11+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Y2(12+I2)=A(59)*GY(2+IG)+A(60)*GY(4+IG)+A(61)*GY(6+IG)
      Y2(13+I2)=A(78)*GY(1+IG)+A(79)*GY(3+IG)+A(80)*GY(5+IG)
      Y2(14+I2)=A(81)*GY(2+IG)+A(82)*GY(4+IG)+A(83)*GY(6+IG)
      Y2(15+I2)=A(84)*GY(1+IG)+A(85)*GY(3+IG)+A(86)*GY(5+IG)+
     $          A(87)*GY(7+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A(16)*GZ(2+IG)
      Z2( 5+I2)=A(17)*GZ(3+IG)
      Z2( 6+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 7+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2( 8+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2( 9+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      Z2(10+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      Z2(11+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
      Z2(12+I2)=A(59)*GZ(2+IG)+A(60)*GZ(4+IG)+A(61)*GZ(6+IG)
      Z2(13+I2)=A(78)*GZ(1+IG)+A(79)*GZ(3+IG)+A(80)*GZ(5+IG)
      Z2(14+I2)=A(81)*GZ(2+IG)+A(82)*GZ(4+IG)+A(83)*GZ(6+IG)
  430 Z2(15+I2)=A(84)*GZ(1+IG)+A(85)*GZ(3+IG)+A(86)*GZ(5+IG)+
     $          A(87)*GZ(7+IG)
      GO TO 500
C
C.....  CASE (4,3)
CDIR$ IVDEP
143   I2=-20
      IG=-8
      DO 440 I=1,NZERO
      I2=I2+20
      IG=IG+8
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2( 5+I2)=A(16)*GX(2+IG)
      X2( 6+I2)=A(17)*GX(3+IG)
      X2( 7+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 8+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      X2( 9+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(10+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(11+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(12+I2)=A(39)*GX(2+IG)+A(40)*GX(4+IG)+A(41)*GX(6+IG)
      X2(13+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(14+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      X2(15+I2)=A(59)*GX(2+IG)+A(60)*GX(4+IG)+A(61)*GX(6+IG)
      X2(16+I2)=A(62)*GX(3+IG)+A(63)*GX(5+IG)+A(64)*GX(7+IG)
      X2(17+I2)=A(78)*GX(1+IG)+A(79)*GX(3+IG)+A(80)*GX(5+IG)
      X2(18+I2)=A(81)*GX(2+IG)+A(82)*GX(4+IG)+A(83)*GX(6+IG)
      X2(19+I2)=A(84)*GX(1+IG)+A(85)*GX(3+IG)+A(86)*GX(5+IG)+
     $          A(87)*GX(7+IG)
      X2(20+I2)=A(88)*GX(2+IG)+A(89)*GX(4+IG)+A(90)*GX(6+IG)+
     $          A(91)*GX(8+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2( 5+I2)=A(16)*GY(2+IG)
      Y2( 6+I2)=A(17)*GY(3+IG)
      Y2( 7+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 8+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Y2( 9+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(10+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(11+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(12+I2)=A(39)*GY(2+IG)+A(40)*GY(4+IG)+A(41)*GY(6+IG)
      Y2(13+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(14+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Y2(15+I2)=A(59)*GY(2+IG)+A(60)*GY(4+IG)+A(61)*GY(6+IG)
      Y2(16+I2)=A(62)*GY(3+IG)+A(63)*GY(5+IG)+A(64)*GY(7+IG)
      Y2(17+I2)=A(78)*GY(1+IG)+A(79)*GY(3+IG)+A(80)*GY(5+IG)
      Y2(18+I2)=A(81)*GY(2+IG)+A(82)*GY(4+IG)+A(83)*GY(6+IG)
      Y2(19+I2)=A(84)*GY(1+IG)+A(85)*GY(3+IG)+A(86)*GY(5+IG)+
     $          A(87)*GY(7+IG)
      Y2(20+I2)=A(88)*GY(2+IG)+A(89)*GY(4+IG)+A(90)*GY(6+IG)+
     $          A(91)*GY(8+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2( 5+I2)=A(16)*GZ(2+IG)
      Z2( 6+I2)=A(17)*GZ(3+IG)
      Z2( 7+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 8+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
      Z2( 9+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(10+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2(11+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      Z2(12+I2)=A(39)*GZ(2+IG)+A(40)*GZ(4+IG)+A(41)*GZ(6+IG)
      Z2(13+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      Z2(14+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
      Z2(15+I2)=A(59)*GZ(2+IG)+A(60)*GZ(4+IG)+A(61)*GZ(6+IG)
      Z2(16+I2)=A(62)*GZ(3+IG)+A(63)*GZ(5+IG)+A(64)*GZ(7+IG)
      Z2(17+I2)=A(78)*GZ(1+IG)+A(79)*GZ(3+IG)+A(80)*GZ(5+IG)
      Z2(18+I2)=A(81)*GZ(2+IG)+A(82)*GZ(4+IG)+A(83)*GZ(6+IG)
      Z2(19+I2)=A(84)*GZ(1+IG)+A(85)*GZ(3+IG)+A(86)*GZ(5+IG)+
     $          A(87)*GZ(7+IG)
  440 Z2(20+I2)=A(88)*GZ(2+IG)+A(89)*GZ(4+IG)+A(90)*GZ(6+IG)+
     $          A(91)*GZ(8+IG)
      GO TO 500
C
C.....  CASE (4,4)
CDIR$ IVDEP
144   I2=-25
      IG=-9
      DO 450 I=1,NZERO
      I2=I2+25
      IG=IG+9
      X2( 1+I2)=      GX(1+IG)
      X2( 2+I2)=A( 1)*GX(2+IG)
      X2( 3+I2)=A( 2)*GX(1+IG)+A( 3)*GX(3+IG)
      X2( 4+I2)=A( 4)*GX(2+IG)+A( 5)*GX(4+IG)
      X2( 5+I2)=A(6)*GX(1+IG)+A( 7)*GX(3+IG)+A( 8)*GX(5+IG)
      X2( 6+I2)=A(16)*GX(2+IG)
      X2( 7+I2)=A(17)*GX(3+IG)
      X2( 8+I2)=A(18)*GX(2+IG)+A(19)*GX(4+IG)
      X2( 9+I2)=A(20)*GX(3+IG)+A(21)*GX(5+IG)
      X2(10+I2)=A(22)*GX(2+IG)+A(23)*GX(4+IG)+A(24)*GX(6+IG)
      X2(11+I2)=A(32)*GX(1+IG)+A(33)*GX(3+IG)
      X2(12+I2)=A(34)*GX(2+IG)+A(35)*GX(4+IG)
      X2(13+I2)=A(36)*GX(1+IG)+A(37)*GX(3+IG)+A(38)*GX(5+IG)
      X2(14+I2)=A(39)*GX(2+IG)+A(40)*GX(4+IG)+A(41)*GX(6+IG)
      X2(15+I2)=A(42)*GX(1+IG)+A(43)*GX(3+IG)+A(44)*GX(5+IG)+
     $          A(45)*GX(7+IG)
      X2(16+I2)=A(55)*GX(2+IG)+A(56)*GX(4+IG)
      X2(17+I2)=A(57)*GX(3+IG)+A(58)*GX(5+IG)
      X2(18+I2)=A(59)*GX(2+IG)+A(60)*GX(4+IG)+A(61)*GX(6+IG)
      X2(19+I2)=A(62)*GX(3+IG)+A(63)*GX(5+IG)+A(64)*GX(7+IG)
      X2(20+I2)=A(65)*GX(2+IG)+A(66)*GX(4+IG)+A(67)*GX(6+IG)+
     $          A(68)*GX(8+IG)
      X2(21+I2)=A(78)*GX(1+IG)+A(79)*GX(3+IG)+A(80)*GX(5+IG)
      X2(22+I2)=A(81)*GX(2+IG)+A(82)*GX(4+IG)+A(83)*GX(6+IG)
      X2(23+I2)=A(84)*GX(1+IG)+A(85)*GX(3+IG)+A(86)*GX(5+IG)+
     $          A(87)*GX(7+IG)
      X2(24+I2)=A(88)*GX(2+IG)+A(89)*GX(4+IG)+A(90)*GX(6+IG)+
     $          A(91)*GX(8+IG)
      X2(25+I2)=A(92)*GX(1+IG)+A(93)*GX(3+IG)+A(94)*GX(5+IG)+
     $          A(95)*GX(7+IG)+A(96)*GX(9+IG)
      Y2( 1+I2)=      GY(1+IG)
      Y2( 2+I2)=A( 1)*GY(2+IG)
      Y2( 3+I2)=A( 2)*GY(1+IG)+A( 3)*GY(3+IG)
      Y2( 4+I2)=A( 4)*GY(2+IG)+A( 5)*GY(4+IG)
      Y2( 5+I2)=A(6)*GY(1+IG)+A( 7)*GY(3+IG)+A( 8)*GY(5+IG)
      Y2( 6+I2)=A(16)*GY(2+IG)
      Y2( 7+I2)=A(17)*GY(3+IG)
      Y2( 8+I2)=A(18)*GY(2+IG)+A(19)*GY(4+IG)
      Y2( 9+I2)=A(20)*GY(3+IG)+A(21)*GY(5+IG)
      Y2(10+I2)=A(22)*GY(2+IG)+A(23)*GY(4+IG)+A(24)*GY(6+IG)
      Y2(11+I2)=A(32)*GY(1+IG)+A(33)*GY(3+IG)
      Y2(12+I2)=A(34)*GY(2+IG)+A(35)*GY(4+IG)
      Y2(13+I2)=A(36)*GY(1+IG)+A(37)*GY(3+IG)+A(38)*GY(5+IG)
      Y2(14+I2)=A(39)*GY(2+IG)+A(40)*GY(4+IG)+A(41)*GY(6+IG)
      Y2(15+I2)=A(42)*GY(1+IG)+A(43)*GY(3+IG)+A(44)*GY(5+IG)+
     $          A(45)*GY(7+IG)
      Y2(16+I2)=A(55)*GY(2+IG)+A(56)*GY(4+IG)
      Y2(17+I2)=A(57)*GY(3+IG)+A(58)*GY(5+IG)
      Y2(18+I2)=A(59)*GY(2+IG)+A(60)*GY(4+IG)+A(61)*GY(6+IG)
      Y2(19+I2)=A(62)*GY(3+IG)+A(63)*GY(5+IG)+A(64)*GY(7+IG)
      Y2(20+I2)=A(65)*GY(2+IG)+A(66)*GY(4+IG)+A(67)*GY(6+IG)+
     $          A(68)*GY(8+IG)
      Y2(21+I2)=A(78)*GY(1+IG)+A(79)*GY(3+IG)+A(80)*GY(5+IG)
      Y2(22+I2)=A(81)*GY(2+IG)+A(82)*GY(4+IG)+A(83)*GY(6+IG)
      Y2(23+I2)=A(84)*GY(1+IG)+A(85)*GY(3+IG)+A(86)*GY(5+IG)+
     $          A(87)*GY(7+IG)
      Y2(24+I2)=A(88)*GY(2+IG)+A(89)*GY(4+IG)+A(90)*GY(6+IG)+
     $          A(91)*GY(8+IG)
      Y2(25+I2)=A(92)*GY(1+IG)+A(93)*GY(3+IG)+A(94)*GY(5+IG)+
     $          A(95)*GY(7+IG)+A(96)*GY(9+IG)
      Z2( 1+I2)=      GZ(1+IG)
      Z2( 2+I2)=A( 1)*GZ(2+IG)
      Z2( 3+I2)=A( 2)*GZ(1+IG)+A( 3)*GZ(3+IG)
      Z2( 4+I2)=A( 4)*GZ(2+IG)+A( 5)*GZ(4+IG)
      Z2( 5+I2)=A(6)*GZ(1+IG)+A( 7)*GZ(3+IG)+A( 8)*GZ(5+IG)
      Z2( 6+I2)=A(16)*GZ(2+IG)
      Z2( 7+I2)=A(17)*GZ(3+IG)
      Z2( 8+I2)=A(18)*GZ(2+IG)+A(19)*GZ(4+IG)
      Z2( 9+I2)=A(20)*GZ(3+IG)+A(21)*GZ(5+IG)
      Z2(10+I2)=A(22)*GZ(2+IG)+A(23)*GZ(4+IG)+A(24)*GZ(6+IG)
      Z2(11+I2)=A(32)*GZ(1+IG)+A(33)*GZ(3+IG)
      Z2(12+I2)=A(34)*GZ(2+IG)+A(35)*GZ(4+IG)
      Z2(13+I2)=A(36)*GZ(1+IG)+A(37)*GZ(3+IG)+A(38)*GZ(5+IG)
      Z2(14+I2)=A(39)*GZ(2+IG)+A(40)*GZ(4+IG)+A(41)*GZ(6+IG)
      Z2(15+I2)=A(42)*GZ(1+IG)+A(43)*GZ(3+IG)+A(44)*GZ(5+IG)+
     $          A(45)*GZ(7+IG)
      Z2(16+I2)=A(55)*GZ(2+IG)+A(56)*GZ(4+IG)
      Z2(17+I2)=A(57)*GZ(3+IG)+A(58)*GZ(5+IG)
      Z2(18+I2)=A(59)*GZ(2+IG)+A(60)*GZ(4+IG)+A(61)*GZ(6+IG)
      Z2(19+I2)=A(62)*GZ(3+IG)+A(63)*GZ(5+IG)+A(64)*GZ(7+IG)
      Z2(20+I2)=A(65)*GZ(2+IG)+A(66)*GZ(4+IG)+A(67)*GZ(6+IG)+
     $          A(68)*GZ(8+IG)
      Z2(21+I2)=A(78)*GZ(1+IG)+A(79)*GZ(3+IG)+A(80)*GZ(5+IG)
      Z2(22+I2)=A(81)*GZ(2+IG)+A(82)*GZ(4+IG)+A(83)*GZ(6+IG)
      Z2(23+I2)=A(84)*GZ(1+IG)+A(85)*GZ(3+IG)+A(86)*GZ(5+IG)+
     $          A(87)*GZ(7+IG)
      Z2(24+I2)=A(88)*GZ(2+IG)+A(89)*GZ(4+IG)+A(90)*GZ(6+IG)+
     $          A(91)*GZ(8+IG)
  450 Z2(25+I2)=A(92)*GZ(1+IG)+A(93)*GZ(3+IG)+A(94)*GZ(5+IG)+
     $          A(95)*GZ(7+IG)+A(96)*GZ(9+IG)
      GO TO 500
C
  500 RETURN
      END
      SUBROUTINE F03CTR(NZERO)
C***********************************************************************
C     ROUTINE TO FORM TWO-DIMENSION 3-CENTER INTEGRALS
C     GIVEN 2-CENTER INTEGRALS AND CCQX, CCQY AND CCQZ.
C     THIS ROUTINE FILLS LDMAX*LCMAX*LPMAX*NZERO LOCATIONS
C     IN EACH OF X3, Y3 AND Z3.
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/CCPQ/CCPX(48),CCPY(48),CCPZ(48),CCQX(48),CCQY(48),CCQZ(48)
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/FST2D0/GX(45),GY(45),GZ(45),
     $              X2(125),Y2(125),Z2(125),
     $              X3(225),Y3(225),Z3(225),
     $              X4(405),Y4(405),Z4(405)
C
C     (LCMAX,LDMAX) INDEPENDENT INITIALIZATION.
      N2=LPMAX*LQMAX
      N3=LCMAX*LDMAX
      I2=-N2
      I3=-N3
      NZP=NZERO*LPMAX
C
C     BRANCH TO APPROPRIATE CODE.
      IGO=(LCMAX-1)*3+LDMAX
      GO TO(100,101,102,110,111,112,120,121,122),IGO
C
C.....  CASE (0,0) ... COPY INTEGRALS.
100   DO 210 I=1,NZP
      X3(I)=X2(I)
      Y3(I)=Y2(I)
  210 Z3(I)=Z2(I)
      GO TO 300
C
C.....  CASE (0,1).
CDIR$ IVDEP
101   I2=-2
      I3=-2
      DO 220 I=1,NZP
      I2=I2+2
      I3=I3+2
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 1)*X2( 1+I2)+         X2( 2+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 1)*Y2( 1+I2)+         Y2( 2+I2)
      Z3( 1+I3)=         Z2( 1+I2)
  220 Z3( 2+I3)=CCQZ( 1)*Z2( 1+I2)+         Z2( 2+I2)
      GO TO 300
C
C.....  CASE (0,2).
CDIR$ IVDEP
102   I2=-3
      I3=-3
      DO 230 I=1,NZP
      I2=I2+3
      I3=I3+3
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 1)*X2( 1+I2)+         X2( 2+I2)
      X3( 3+I3)=CCQX( 2)*X2( 1+I2)+CCQX( 3)*X2( 2+I2)+
     $                   X2( 3+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 1)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 3+I3)=CCQY( 2)*Y2( 1+I2)+CCQY( 3)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Z3( 1+I3)=         Z2( 1+I2)
      Z3( 2+I3)=CCQZ( 1)*Z2( 1+I2)+         Z2( 2+I2)
  230 Z3( 3+I3)=CCQZ( 2)*Z2( 1+I2)+CCQZ( 3)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      GO TO 300
C
C.....  CASE (1,0).
CDIR$ IVDEP
110   I2=-2
      I3=-2
      DO 240 I=1,NZP
      I2=I2+2
      I3=I3+2
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 7)*X2( 1+I2)+         X2( 2+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 7)*Y2( 1+I2)+         Y2( 2+I2)
      Z3( 1+I3)=         Z2( 1+I2)
  240 Z3( 2+I3)=CCQZ( 7)*Z2( 1+I2)+         Z2( 2+I2)
      GO TO 300
C
C.....  CASE (1,1).
CDIR$ IVDEP
111   I2=-3
      I3=-4
      DO 250 I=1,NZP
      I2=I2+3
      I3=I3+4
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 1)*X2( 1+I2)+         X2( 2+I2)
      X3( 3+I3)=CCQX( 7)*X2( 1+I2)+         X2( 2+I2)
      X3( 4+I3)=CCQX( 8)*X2( 1+I2)+CCQX( 9)*X2( 2+I2)+
     $                   X2( 3+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 1)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 3+I3)=CCQY( 7)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 4+I3)=CCQY( 8)*Y2( 1+I2)+CCQY( 9)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Z3( 1+I3)=         Z2( 1+I2)
      Z3( 2+I3)=CCQZ( 1)*Z2( 1+I2)+         Z2( 2+I2)
      Z3( 3+I3)=CCQZ( 7)*Z2( 1+I2)+         Z2( 2+I2)
  250 Z3( 4+I3)=CCQZ( 8)*Z2( 1+I2)+CCQZ( 9)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      GO TO 300
C
C.....  CASE (1,2).
CDIR$ IVDEP
112   I2=-4
      I3=-6
      DO 260 I=1,NZP
      I2=I2+4
      I3=I3+6
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 1)*X2( 1+I2)+         X2( 2+I2)
      X3( 3+I3)=CCQX( 2)*X2( 1+I2)+CCQX( 3)*X2( 2+I2)+
     $                   X2( 3+I2)
      X3( 4+I3)=CCQX( 7)*X2( 1+I2)+         X2( 2+I2)
      X3( 5+I3)=CCQX( 8)*X2( 1+I2)+CCQX( 9)*X2( 2+I2)+
     $                   X2( 3+I2)
      X3( 6+I3)=CCQX(10)*X2( 1+I2)+CCQX(11)*X2( 2+I2)+
     $          CCQX(12)*X2( 3+I2)+         X2( 4+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 1)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 3+I3)=CCQY( 2)*Y2( 1+I2)+CCQY( 3)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Y3( 4+I3)=CCQY( 7)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 5+I3)=CCQY( 8)*Y2( 1+I2)+CCQY( 9)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Y3( 6+I3)=CCQY(10)*Y2( 1+I2)+CCQY(11)*Y2( 2+I2)+
     $          CCQY(12)*Y2( 3+I2)+         Y2( 4+I2)
      Z3( 1+I3)=         Z2( 1+I2)
      Z3( 2+I3)=CCQZ( 1)*Z2( 1+I2)+         Z2( 2+I2)
      Z3( 3+I3)=CCQZ( 2)*Z2( 1+I2)+CCQZ( 3)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      Z3( 4+I3)=CCQZ( 7)*Z2( 1+I2)+         Z2( 2+I2)
      Z3( 5+I3)=CCQZ( 8)*Z2( 1+I2)+CCQZ( 9)*Z2( 2+I2)+
     $                   Z2( 3+I2)
  260 Z3( 6+I3)=CCQZ(10)*Z2( 1+I2)+CCQZ(11)*Z2( 2+I2)+
     $          CCQZ(12)*Z2( 3+I2)+         Z2( 4+I2)
      GO TO 300
C
C.....  CASE (2,0).
CDIR$ IVDEP
120   I2=-3
      I3=-3
      DO 270 I=1,NZP
      I2=I2+3
      I3=I3+3
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 7)*X2( 1+I2)+         X2( 2+I2)
      X3( 3+I3)=CCQX(17)*X2( 1+I2)+CCQX(18)*X2( 2+I2)+
     $                   X2( 3+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 7)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 3+I3)=CCQY(17)*Y2( 1+I2)+CCQY(18)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Z3( 1+I3)=         Z2( 1+I2)
      Z3( 2+I3)=CCQZ( 7)*Z2( 1+I2)+         Z2( 2+I2)
  270 Z3( 3+I3)=CCQZ(17)*Z2( 1+I2)+CCQZ(18)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      GO TO 300
C
C.....  CASE (2,1).
CDIR$ IVDEP
121   I2=-4
      I3=-6
      DO 280 I=1,NZP
      I2=I2+4
      I3=I3+6
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 1)*X2( 1+I2)+         X2( 2+I2)
      X3( 3+I3)=CCQX( 7)*X2( 1+I2)+         X2( 2+I2)
      X3( 4+I3)=CCQX( 8)*X2( 1+I2)+CCQX( 9)*X2( 2+I2)+
     $                   X2( 3+I2)
      X3( 5+I3)=CCQX(17)*X2( 1+I2)+CCQX(18)*X2( 2+I2)+
     $                   X2( 3+I2)
      X3( 6+I3)=CCQX(19)*X2( 1+I2)+CCQX(20)*X2( 2+I2)+
     $          CCQX(21)*X2( 3+I2)+         X2( 4+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 1)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 3+I3)=CCQY( 7)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 4+I3)=CCQY( 8)*Y2( 1+I2)+CCQY( 9)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Y3( 5+I3)=CCQY(17)*Y2( 1+I2)+CCQY(18)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Y3( 6+I3)=CCQY(19)*Y2( 1+I2)+CCQY(20)*Y2( 2+I2)+
     $          CCQY(21)*Y2( 3+I2)+         Y2( 4+I2)
      Z3( 1+I3)=         Z2( 1+I2)
      Z3( 2+I3)=CCQZ( 1)*Z2( 1+I2)+         Z2( 2+I2)
      Z3( 3+I3)=CCQZ( 7)*Z2( 1+I2)+         Z2( 2+I2)
      Z3( 4+I3)=CCQZ( 8)*Z2( 1+I2)+CCQZ( 9)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      Z3( 5+I3)=CCQZ(17)*Z2( 1+I2)+CCQZ(18)*Z2( 2+I2)+
     $                   Z2( 3+I2)
  280 Z3( 6+I3)=CCQZ(19)*Z2( 1+I2)+CCQZ(20)*Z2( 2+I2)+
     $          CCQZ(21)*Z2( 3+I2)+         Z2( 4+I2)
      GO TO 300
C
C.....  CASE (2,2).
CDIR$ IVDEP
122   I2=-5
      I3=-9
      DO 290 I=1,NZP
      I2=I2+5
      I3=I3+9
      X3( 1+I3)=         X2( 1+I2)
      X3( 2+I3)=CCQX( 1)*X2( 1+I2)+         X2( 2+I2)
      X3( 3+I3)=CCQX( 2)*X2( 1+I2)+CCQX( 3)*X2( 2+I2)+
     $                   X2( 3+I2)
      X3( 4+I3)=CCQX( 7)*X2( 1+I2)+         X2( 2+I2)
      X3( 5+I3)=CCQX( 8)*X2( 1+I2)+CCQX( 9)*X2( 2+I2)+
     $                   X2( 3+I2)
      X3( 6+I3)=CCQX(10)*X2( 1+I2)+CCQX(11)*X2( 2+I2)+
     $          CCQX(12)*X2( 3+I2)+         X2( 4+I2)
      X3( 7+I3)=CCQX(17)*X2( 1+I2)+CCQX(18)*X2( 2+I2)+
     $                   X2( 3+I2)
      X3( 8+I3)=CCQX(19)*X2( 1+I2)+CCQX(20)*X2( 2+I2)+
     $          CCQX(21)*X2( 3+I2)+         X2( 4+I2)
      X3( 9+I3)=CCQX(22)*X2( 1+I2)+CCQX(23)*X2( 2+I2)+
     $          CCQX(24)*X2( 3+I2)+CCQX(25)*X2( 4+I2)+
     $                   X2( 5+I2)
      Y3( 1+I3)=         Y2( 1+I2)
      Y3( 2+I3)=CCQY( 1)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 3+I3)=CCQY( 2)*Y2( 1+I2)+CCQY( 3)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Y3( 4+I3)=CCQY( 7)*Y2( 1+I2)+         Y2( 2+I2)
      Y3( 5+I3)=CCQY( 8)*Y2( 1+I2)+CCQY( 9)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Y3( 6+I3)=CCQY(10)*Y2( 1+I2)+CCQY(11)*Y2( 2+I2)+
     $          CCQY(12)*Y2( 3+I2)+         Y2( 4+I2)
      Y3( 7+I3)=CCQY(17)*Y2( 1+I2)+CCQY(18)*Y2( 2+I2)+
     $                   Y2( 3+I2)
      Y3( 8+I3)=CCQY(19)*Y2( 1+I2)+CCQY(20)*Y2( 2+I2)+
     $          CCQY(21)*Y2( 3+I2)+         Y2( 4+I2)
      Y3( 9+I3)=CCQY(22)*Y2( 1+I2)+CCQY(23)*Y2( 2+I2)+
     $          CCQY(24)*Y2( 3+I2)+CCQY(25)*Y2( 4+I2)+
     $                   Y2( 5+I2)
      Z3( 1+I3)=         Z2( 1+I2)
      Z3( 2+I3)=CCQZ( 1)*Z2( 1+I2)+         Z2( 2+I2)
      Z3( 3+I3)=CCQZ( 2)*Z2( 1+I2)+CCQZ( 3)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      Z3( 4+I3)=CCQZ( 7)*Z2( 1+I2)+         Z2( 2+I2)
      Z3( 5+I3)=CCQZ( 8)*Z2( 1+I2)+CCQZ( 9)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      Z3( 6+I3)=CCQZ(10)*Z2( 1+I2)+CCQZ(11)*Z2( 2+I2)+
     $          CCQZ(12)*Z2( 3+I2)+         Z2( 4+I2)
      Z3( 7+I3)=CCQZ(17)*Z2( 1+I2)+CCQZ(18)*Z2( 2+I2)+
     $                   Z2( 3+I2)
      Z3( 8+I3)=CCQZ(19)*Z2( 1+I2)+CCQZ(20)*Z2( 2+I2)+
     $          CCQZ(21)*Z2( 3+I2)+         Z2( 4+I2)
  290 Z3( 9+I3)=CCQZ(22)*Z2( 1+I2)+CCQZ(23)*Z2( 2+I2)+
     $          CCQZ(24)*Z2( 3+I2)+CCQZ(25)*Z2( 4+I2)+
     $                   Z2( 5+I2)
      GO TO 300
C
  300 RETURN
      END
      SUBROUTINE F04CTR(NZERO)
C***********************************************************************
C     ROUTINE TO SPLIT THREE-CENTER INTEGRALS INTO FOUR-CENTER
C     INTEGRALS USING CCP.
C     THIS ROUTINE WILL FILL UP TO NZERO*LEN2D LOCATIONS IN
C     THE ARRAYS X4, Y4 AND Z4.
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/CCPQ/CCPX(48),CCPY(48),CCPZ(48),CCQX(48),CCQY(48),CCQZ(48)
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/FST2D0/GX(45),GY(45),GZ(45),
     $              X2(125),Y2(125),Z2(125),
     $              X3(225),Y3(225),Z3(225),
     $              X4(405),Y4(405),Z4(405)
C
C     INITIALIZATION.
      NCD=LCMAX*LDMAX
      N3=LPMAX*LCMAX*LDMAX
      N4=LAMAX*LBMAX*LCMAX*LDMAX
      I3Z=-N3
      I4Z=-N4
C
C     BRANCH TO APPROPRIATE PIECE OF CODE.
      IGO=(LAMAX-1)*3+LBMAX
      GO TO(100,101,102,110,111,112,120,121,122),IGO
C
C.....  CASE (0,0) ... COPY INTEGRALS.
100   NZCD=NZERO*NCD
      DO 210 I=1,NZCD
      X4(I)=X3(I)
      Y4(I)=Y3(I)
  210 Z4(I)=Z3(I)
      GO TO 300
C
C.....  CASE (0,1).
101   NCD1=NCD
      DO 220 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 220 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 1)*X3(I3     )+         X3(I3+NCD1)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 1)*Y3(I3     )+         Y3(I3+NCD1)
      Z4(I4     )=         Z3(I3     )
  220 Z4(I4+NCD1)=CCPZ( 1)*Z3(I3     )+         Z3(I3+NCD1)
      GO TO 300
C
C.....  CASE (0,2).
102   NCD1=NCD
      NCD2=NCD+NCD1
      DO 230 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 230 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 1)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD2)=CCPX( 2)*X3(I3     )+CCPX( 3)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 1)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD2)=CCPY( 2)*Y3(I3     )+CCPY( 3)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Z4(I4     )=         Z3(I3     )
      Z4(I4+NCD1)=CCPZ( 1)*Z3(I3     )+         Z3(I3+NCD1)
  230 Z4(I4+NCD2)=CCPZ( 2)*Z3(I3     )+CCPZ( 3)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      GO TO 300
C
C.....  CASE (1,0).
110   NCD1=NCD
      DO 240 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 240 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 7)*X3(I3     )+         X3(I3+NCD1)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 7)*Y3(I3     )+         Y3(I3+NCD1)
      Z4(I4     )=         Z3(I3     )
  240 Z4(I4+NCD1)=CCPZ( 7)*Z3(I3     )+         Z3(I3+NCD1)
      GO TO 300
C
C.....  CASE (1,1).
111   NCD1=NCD
      NCD2=NCD+NCD1
      NCD3=NCD+NCD2
      DO 250 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 250 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 1)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD2)=CCPX( 7)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD3)=CCPX( 8)*X3(I3     )+CCPX( 9)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 1)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD2)=CCPY( 7)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD3)=CCPY( 8)*Y3(I3     )+CCPY( 9)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Z4(I4     )=         Z3(I3     )
      Z4(I4+NCD1)=CCPZ( 1)*Z3(I3     )+         Z3(I3+NCD1)
      Z4(I4+NCD2)=CCPZ( 7)*Z3(I3     )+         Z3(I3+NCD1)
  250 Z4(I4+NCD3)=CCPZ( 8)*Z3(I3     )+CCPZ( 9)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      GO TO 300
C
C.....  CASE (1,2).
112   NCD1=NCD
      NCD2=NCD+NCD1
      NCD3=NCD+NCD2
      NCD4=NCD+NCD3
      NCD5=NCD+NCD4
      DO 260 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 260 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 1)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD2)=CCPX( 2)*X3(I3     )+CCPX( 3)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      X4(I4+NCD3)=CCPX( 7)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD4)=CCPX( 8)*X3(I3     )+CCPX( 9)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      X4(I4+NCD5)=CCPX(10)*X3(I3     )+CCPX(11)*X3(I3+NCD1)+
     $            CCPX(12)*X3(I3+NCD2)+         X3(I3+NCD3)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 1)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD2)=CCPY( 2)*Y3(I3     )+CCPY( 3)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Y4(I4+NCD3)=CCPY( 7)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD4)=CCPY( 8)*Y3(I3     )+CCPY( 9)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Y4(I4+NCD5)=CCPY(10)*Y3(I3     )+CCPY(11)*Y3(I3+NCD1)+
     $            CCPY(12)*Y3(I3+NCD2)+         Y3(I3+NCD3)
      Z4(I4     )=         Z3(I3     )
      Z4(I4+NCD1)=CCPZ( 1)*Z3(I3     )+         Z3(I3+NCD1)
      Z4(I4+NCD2)=CCPZ( 2)*Z3(I3     )+CCPZ( 3)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      Z4(I4+NCD3)=CCPZ( 7)*Z3(I3     )+         Z3(I3+NCD1)
      Z4(I4+NCD4)=CCPZ( 8)*Z3(I3     )+CCPZ( 9)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
  260 Z4(I4+NCD5)=CCPZ(10)*Z3(I3     )+CCPZ(11)*Z3(I3+NCD1)+
     $            CCPZ(12)*Z3(I3+NCD2)+         Z3(I3+NCD3)
      GO TO 300
C
C.....  CASE (2,0).
120   NCD1=NCD
      NCD2=NCD+NCD1
      DO 270 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 270 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 7)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD2)=CCPX(17)*X3(I3     )+CCPX(18)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 7)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD2)=CCPY(17)*Y3(I3     )+CCPY(18)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Z4(I4     )=         Z3(I3     )
      Z4(I4+NCD1)=CCPZ( 7)*Z3(I3     )+         Z3(I3+NCD1)
  270 Z4(I4+NCD2)=CCPZ(17)*Z3(I3     )+CCPZ(18)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      GO TO 300
C
C.....  CASE (2,1).
121   NCD1=NCD
      NCD2=NCD+NCD1
      NCD3=NCD+NCD2
      NCD4=NCD+NCD3
      NCD5=NCD+NCD4
      DO 280 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 280 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 1)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD2)=CCPX( 7)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD3)=CCPX( 8)*X3(I3     )+CCPX( 9)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      X4(I4+NCD4)=CCPX(17)*X3(I3     )+CCPX(18)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      X4(I4+NCD5)=CCPX(19)*X3(I3     )+CCPX(20)*X3(I3+NCD1)+
     $            CCPX(21)*X3(I3+NCD2)+         X3(I3+NCD3)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 1)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD2)=CCPY( 7)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD3)=CCPY( 8)*Y3(I3     )+CCPY( 9)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Y4(I4+NCD4)=CCPY(17)*Y3(I3     )+CCPY(18)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Y4(I4+NCD5)=CCPY(19)*Y3(I3     )+CCPY(20)*Y3(I3+NCD1)+
     $            CCPY(21)*Y3(I3+NCD2)+         Y3(I3+NCD3)
      Z4(I4     )=         Z3(I3     )
      Z4(I4+NCD1)=CCPZ( 1)*Z3(I3     )+         Z3(I3+NCD1)
      Z4(I4+NCD2)=CCPZ( 7)*Z3(I3     )+         Z3(I3+NCD1)
      Z4(I4+NCD3)=CCPZ( 8)*Z3(I3     )+CCPZ( 9)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      Z4(I4+NCD4)=CCPZ(17)*Z3(I3     )+CCPZ(18)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
  280 Z4(I4+NCD5)=CCPZ(19)*Z3(I3     )+CCPZ(20)*Z3(I3+NCD1)+
     $            CCPZ(21)*Z3(I3+NCD2)+         Z3(I3+NCD3)
      GO TO 300
C
C.....  CASE (2,2).
122   NCD1=NCD
      NCD2=NCD+NCD1
      NCD3=NCD+NCD2
      NCD4=NCD+NCD3
      NCD5=NCD+NCD4
      NCD6=NCD+NCD5
      NCD7=NCD+NCD6
      NCD8=NCD+NCD7
      DO 290 IZ=1,NZERO
      I3Z=I3Z+N3
      I4Z=I4Z+N4
CDIR$ IVDEP
      DO 290 I=1,NCD
      I3=I+I3Z
      I4=I+I4Z
      X4(I4     )=         X3(I3     )
      X4(I4+NCD1)=CCPX( 1)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD2)=CCPX( 2)*X3(I3     )+CCPX( 3)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      X4(I4+NCD3)=CCPX( 7)*X3(I3     )+         X3(I3+NCD1)
      X4(I4+NCD4)=CCPX( 8)*X3(I3     )+CCPX( 9)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      X4(I4+NCD5)=CCPX(10)*X3(I3     )+CCPX(11)*X3(I3+NCD1)+
     $            CCPX(12)*X3(I3+NCD2)+         X3(I3+NCD3)
      X4(I4+NCD6)=CCPX(17)*X3(I3     )+CCPX(18)*X3(I3+NCD1)+
     $                     X3(I3+NCD2)
      X4(I4+NCD7)=CCPX(19)*X3(I3     )+CCPX(20)*X3(I3+NCD1)+
     $            CCPX(21)*X3(I3+NCD2)+         X3(I3+NCD3)
      X4(I4+NCD8)=CCPX(22)*X3(I3     )+CCPX(23)*X3(I3+NCD1)+
     $            CCPX(24)*X3(I3+NCD2)+CCPX(25)*X3(I3+NCD3)+
     $                     X3(I3+NCD4)
      Y4(I4     )=         Y3(I3     )
      Y4(I4+NCD1)=CCPY( 1)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD2)=CCPY( 2)*Y3(I3     )+CCPY( 3)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Y4(I4+NCD3)=CCPY( 7)*Y3(I3     )+         Y3(I3+NCD1)
      Y4(I4+NCD4)=CCPY( 8)*Y3(I3     )+CCPY( 9)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Y4(I4+NCD5)=CCPY(10)*Y3(I3     )+CCPY(11)*Y3(I3+NCD1)+
     $            CCPY(12)*Y3(I3+NCD2)+         Y3(I3+NCD3)
      Y4(I4+NCD6)=CCPY(17)*Y3(I3     )+CCPY(18)*Y3(I3+NCD1)+
     $                     Y3(I3+NCD2)
      Y4(I4+NCD7)=CCPY(19)*Y3(I3     )+CCPY(20)*Y3(I3+NCD1)+
     $            CCPY(21)*Y3(I3+NCD2)+         Y3(I3+NCD3)
      Y4(I4+NCD8)=CCPY(22)*Y3(I3     )+CCPY(23)*Y3(I3+NCD1)+
     $            CCPY(24)*Y3(I3+NCD2)+CCPY(25)*Y3(I3+NCD3)+
     $                     Y3(I3+NCD4)
      Z4(I4     )=         Z3(I3     )
      Z4(I4+NCD1)=CCPZ( 1)*Z3(I3     )+         Z3(I3+NCD1)
      Z4(I4+NCD2)=CCPZ( 2)*Z3(I3     )+CCPZ( 3)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      Z4(I4+NCD3)=CCPZ( 7)*Z3(I3     )+         Z3(I3+NCD1)
      Z4(I4+NCD4)=CCPZ( 8)*Z3(I3     )+CCPZ( 9)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      Z4(I4+NCD5)=CCPZ(10)*Z3(I3     )+CCPZ(11)*Z3(I3+NCD1)+
     $            CCPZ(12)*Z3(I3+NCD2)+         Z3(I3+NCD3)
      Z4(I4+NCD6)=CCPZ(17)*Z3(I3     )+CCPZ(18)*Z3(I3+NCD1)+
     $                     Z3(I3+NCD2)
      Z4(I4+NCD7)=CCPZ(19)*Z3(I3     )+CCPZ(20)*Z3(I3+NCD1)+
     $            CCPZ(21)*Z3(I3+NCD2)+         Z3(I3+NCD3)
  290 Z4(I4+NCD8)=CCPZ(22)*Z3(I3     )+CCPZ(23)*Z3(I3+NCD1)+
     $            CCPZ(24)*Z3(I3+NCD2)+CCPZ(25)*Z3(I3+NCD3)+
     $                     Z3(I3+NCD4)
      GO TO 300
C
  300 RETURN
      END
      SUBROUTINE F0CFIL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NS=#NS, NP=#NP)
C##
      PARAMETER (NS=120, NP=300)
C###
      INTEGER SHELLA,SHELLN,SHELLT,SHELLC,AOS
      INTEGER UBOUND
C***********************************************************************
C     ROUTINE TO PRE-TABULATE CONTRACTION COEFFICIENTS.
C***********************************************************************
      COMMON/B/IXYZ(NS),SHELLA(NS),SHELLN(NS),SHELLT(NS),
     1 SHELLC(NS),AOS(NS),NSHELL,MAXTYP,EXX(NP),C1(NP),C2(NP)
      COMMON/GCLOOP/IGAUSS,IGBEG,IGEND,IGDF,JGAUSS,JGBEG,JGEND,JGDF,
     $              KGAUSS,KGBEG,KGEND,KGDF,LGAUSS,LGBEG,LGEND,LGDF
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/LIMIT2/ITYPE,JTYPE,KTYPE,LTYPE,
     1              IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,LSTART,
     2              IEND,JEND,KEND,LEND,LENTQ,
     3              R3OV2,ROOT3,ROOT5,ROOT15,R1,R2,R4,Z1,Z2,Z3
      COMMON/FSTCTR/CKLDAT(400),CIJ(100),CKL(100),KLCOND(100)
C
      DIMENSION CC(20),CD(20),UBOUND(4)
C
      DATA UBOUND/1,4,10,20/
C
C     FILL KLCOND (USING CLEVER CODE FROM HONDO).
      INTC=0
      KL=0
      KLIM=UBOUND(LCMAX)
      DO 5 K=KSTART,KLIM
      GO TO(1,1,2,2,1,2,2,1,2,2),K
    1 KLSAVE=KL
    2 KL=KLSAVE
      LEND=UBOUND(LDMAX)
      IF(KML.EQ.0)LEND=K
      DO 5 L=LSTART,LEND
      GO TO(3,3,4,4,3,4,4,3,4,4),L
    3 KL=KL+1
    4 INTC=INTC+1
    5 KLCOND(INTC)=KL
C
C     PACK OUT DESIRED CONTRACTION COEFFICIENTS.
      KLIND=0
      DO 60 KGAUSS=KGBEG,KGEND
      CALL FILLC(KTYPE,KGAUSS,C1,C2,CC)
      DO 60 LGAUSS=LGBEG,LGEND
      KL=KLIND*4
      KLIND=KLIND+1
      CALL FILLC(LTYPE,LGAUSS,C1,C2,CD)
      KLIM=UBOUND(LCMAX)
      DO 55 K=KSTART,KLIM
      GO TO(30,30,55,55,30,55,55,30,55,55),K
   30 LEND=UBOUND(LDMAX)
      IF(KML.EQ.0)LEND=K
      DO 50 L=LSTART,LEND
      GO TO(40,40,50,50,40,50,50,40,50,50),L
   40 KL=KL+1
      CKLDAT(KL)=CC(K)*CD(L)
   50 CONTINUE
   55 CONTINUE
   60 CONTINUE
C
      RETURN
      END
      SUBROUTINE F0CLD1(CA,CB)
C***********************************************************************
C     ROUTINE TO LOAD (IGAUSS,JGUASS) CONTRACTION COEFFICIENTS
C     INTO CIJ.
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/LIMIT2/ITYPE,JTYPE,KTYPE,LTYPE,
     1              IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,LSTART,
     2              IEND,JEND,KEND,LEND,LENTQ,
     3              R3OV2,ROOT3,ROOT5,ROOT15,R1,R2,R4,Z1,Z2,Z3
      COMMON/FSTCTR/CKLDAT(400),CIJ(100),CKL(100),KLCOND(100)
C
      DIMENSION CA(20),CB(20)
C
C     LOAD OUT (IGAUSS,JGAUSS) COEFICIENTS.
      IJ=0
      DO 10 I=ISTART,IEND
      IF(IMJ.EQ.0)JEND=I
      DO 10 J=JSTART,JEND
      IJ=IJ+1
      CIJ(IJ)=CA(I)*CB(J)
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE F0CLD2(KLIND)
C***********************************************************************
C     ROUTINE TO LOAD (KGAUSS,LGAUSS) CONTRACTION COEFICIENTS
C     OUT OF CKLDAT INTO CKL.
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/FSTIND/IJXDAT(100),IJYDAT(100),IJZDAT(100),
     $              KLXDAT(100),KLYDAT(100),KLZDAT(100),
     $              NKLDAT(100),NIJ,MAXNKL
      COMMON/FSTCTR/CKLDAT(400),CIJ(100),CKL(100),KLCOND(100)
C
C     LOAD OUT (KGAUSS,LGAUSS) CONTRACTION COEFFICIENTS.
      KL4=(KLIND-1)*4
      DO 10 KL=1,MAXNKL
      II=KLCOND(KL)+KL4
   10 CKL(KL)=CKLDAT(II)
C
      RETURN
      END
      LOGICAL FUNCTION F0INIT(KNTFST)
C
C     --------------
C     GAUSSIAN 82
C     U OF T VERSION
C     FEBRUARY 1987
C     --------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER UBOUND
C***********************************************************************
C     ROUTINE TO
C     1.  DETERMINE IF 'FAST' CODE CAN BE USED, AND IF SO
C     2.  PERFORM VARIOUS INITIALIZATION CHORES
C         A.  INITIALIZE POINTER ARRAYS.
C***********************************************************************
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/LIMIT2/ITYPE,JTYPE,KTYPE,LTYPE,
     1              IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,LSTART,
     2              IEND,JEND,KEND,LEND,LENTQ,
     3              R3OV2,ROOT3,ROOT5,ROOT15,R1,R2,R4,Z1,Z2,Z3
      COMMON /SHLCOM/ ISHELL,JSHELL,KSHELL,LSHELL,IRANGE,JRANGE,KRANGE,
     1 LRANGE
      COMMON/FSTIND/IJXDAT(100),IJYDAT(100),IJZDAT(100),
     $              KLXDAT(100),KLYDAT(100),KLZDAT(100),
     $              NKLDAT(100),NIJ,MAXNKL
C     COMMON/FSTAO/INTX(1296),INTY(1296),INTZ(1296),NIJKL
C
      DIMENSION INDLX(20),INDLY(20),INDLZ(20),UBOUND(4)
C
C     LOGICAL FAST
C
      DATA INDLX/1,2,1,1,3,1,1,2,2,1,4,1,1,2,3,3,2,1,1,2/
      DATA INDLY/1,1,2,1,1,3,1,2,1,2,1,4,1,3,2,1,1,2,3,2/
      DATA INDLZ/1,1,1,2,1,1,3,1,2,2,1,1,4,1,1,2,3,3,2,2/
      DATA UBOUND/1,4,10,20/
C     DATA FAST/.FALSE./
C
C     DETERMINE ACTION.
      F0INIT=.FALSE.
      IF(LPMAX.GT.5)RETURN
      IF(LQMAX.GT.5)RETURN
      IF(LAMAX.GT.3)RETURN
      IF(LBMAX.GT.3)RETURN
      IF(LCMAX.GT.3)RETURN
      IF(LDMAX.GT.3)RETURN
      IF(IRANGE.GE.10)RETURN
      IF(JRANGE.GE.10)RETURN
      IF(KRANGE.GE.10)RETURN
      IF(LRANGE.GE.10)RETURN
      F0INIT=.TRUE.
      KNTFST=KNTFST+1
C
C     FILL INDEXING ARRAYS TO BE USED IN SUMMING UP INTEGRALS.
C     FIRST, KLDAT ...
      NKL=0
      KLIM=UBOUND(LCMAX)
      DO 10 K=KSTART,KLIM
      LEND=UBOUND(LDMAX)
      IF(KML.EQ.0)LEND=K
      DO 10 L=LSTART,LEND
      NKL=NKL+1
      KLXDAT(NKL)=(INDLX(K)-1)*LDMAX+INDLX(L)
      KLYDAT(NKL)=(INDLY(K)-1)*LDMAX+INDLY(L)
   10 KLZDAT(NKL)=(INDLZ(K)-1)*LDMAX+INDLZ(L)
      MAXNKL=NKL
C     THEN, IJDAT ...
      KLR=LCMAX*LDMAX
      NIJ=0
      DO 20 I=ISTART,IEND
      IF(IMJ.EQ.0)JEND=I
      DO 20 J=JSTART,JEND
      NIJ=NIJ+1
      IJXDAT(NIJ)=(((INDLX(I)-1)*LBMAX+INDLX(J))-1)*KLR
      IJYDAT(NIJ)=(((INDLY(I)-1)*LBMAX+INDLY(J))-1)*KLR
   20 IJZDAT(NIJ)=(((INDLZ(I)-1)*LBMAX+INDLZ(J))-1)*KLR
C     DETERMINE LIMIT FOR KL PASSES.
      IF(IMKJML.EQ.0)GO TO 40
C        (ISH,JSH) DIFFERENT FROM (KSH,LSH).
         DO 30 IJ=1,NIJ
   30    NKLDAT(IJ)=NKL
         GO TO 60
C        (ISH,JSH) = (KSH,LSH).
   40    DO 50 IJ=1,NIJ
   50    NKLDAT(IJ)=IJ
   60 CONTINUE
C
C     PREPARE MACHINE-DEPENDENT GATHER ARRAYS.
C     CALL MACHINE-DEPENDENT ROUTINE TO PACK THE INFORMATION
C     AS APPROPRIATE FOR THE LOCAL ENVIRONMENT.
C
C     IF(.NOT.FAST) RETURN
C     NIJKL=0
C     DO 70 IJ=1,NIJ
C     IJX=IJXDAT(IJ)
C     IJY=IJYDAT(IJ)
C     IJZ=IJZDAT(IJ)
C     NKL=NKLDAT(IJ)
C     DO 70 KL=1,NKL
C     NIJKL=NIJKL+1
C     INTX(NIJKL)=KLXDAT(KL)+IJX
C     INTY(NIJKL)=KLYDAT(KL)+IJY
C  70 INTZ(NIJKL)=KLZDAT(KL)+IJZ
C     CALL IPTXYZ(NIJKL,INTX,INTY,INTZ)
      RETURN
      END
