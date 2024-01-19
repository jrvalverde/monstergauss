C     GL0307A      07 SEP 89                                         MRP
C?IBM/GLD/GBR/VAX/UNX
      SUBROUTINE UNCON
C??
C?CDC
C     PROGRAM UNCON
C??
C1INTRODUCTION
C     LINK 0307
C*
C     ----------------
C     GAUSSIAN 82
C     U OF T VERSION
C     SEPTEMBER 1989
C     ----------------
C*
C     TWO-ELECTRON INTEGRAL PACKAGE FOR SPD FUNCTIONS.
C     NORMALLY NO INTEGRALS ARE EVALUATED HERE AS THE PHOENIX ROUTINE
C     (LINK 0308), USING RYS POLYNOMIALS, IS FASTER.
C
C     THE MODE OF OPERATION IS DETERMINED BY IOP(66) - IF PROGRAM SHELL
C     WAS CALLED DURING THIS OVERLAY, UNCON ASSUMES ALL NON-D INTEGRALS
C     HAVE ALREADY BEEN EVALUATED. HOWEVER, IF IOP(66) IS 0, THEN ALL
C     THE INTEGRALS ARE EVALUATED HERE - THUS IT IS VERY IMPORTANT TO
C     CALL LINKS 0306 AND 0307 FROM THE SAME OVERLAY CARD IN THE ROUTE.
C     FOR BASIS SETS CONTAINING F ORBITALS, NO INTEGRALS MAY BE
C     EVALUATED HERE - ALL D/F INTEGRALS MUST BE DONE IN LINK 0308.
C*
C1OPTIONS
C     ******************************************************************
C     OPTIONS ... IOP() ... SEE PROGRAM GINPUT (LINK 0301).
C     ******************************************************************
C*
C/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA)
C     PARAMETER (NS=#NS, NP=#NP)
C##
      PARAMETER (NA= 36)
      PARAMETER (NS=120, NP=300)
C###
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
      INTEGER SHELLA,SHELLN,SHELLT,SHELLC,AOS,UBOUND
C
      COMMON /A/ IOP(99),ID1(NA),ID2(7),C(NA,3),CD(4),ID3(401)
C
      COMMON/B/IXYZ(NS),SHELLA(NS),SHELLN(NS),SHELLT(NS),
     1 SHELLC(NS),AOS(NS),NSHELL,MAXTYP,EXX(NP),C1(NP),C2(NP)
C
      COMMON/C0307A/F(9)
      COMMON/C0307A/XYZP(3),EP,RAB2,OFAB,
     1              XYZQ(3),EQ,RCD2,OFCD,
     2              PQX,PQY,PQZ,RPQ
      COMMON/C0307A/SAB(9),SCD(9),R(3,3),D(6,6),ET(35),SA(73),SB(73),
     1 SC(73),SD(73),TW(197),TQ(1296),FT(2215),CFILL(1296)
      COMMON/C0307B/R3OV2,IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,
     1 LSTART,IEND,JEND,KEND,LEND,MEND,NUMD,IT,JT,KT,LT,INTC,
     2 IS1(30),SAS(3),SBS(3),SCS(3),SDS(3),XYZA(3),XYZB(3),
     3 XYZC(3),XYZD(3),RA(10)
      COMMON /SHLCOM/ ISHELL,JSHELL,KSHELL,LSHELL,IRANGE,JRANGE,KRANGE,
     1 LRANGE
      COMMON/IO/INN,IOUT,IODUM(215)
C*
      DIMENSION IJKL(10,10),IJKLAO(10),LBOUND(3),UBOUND(3)
C
      EQUIVALENCE (XA,XYZA(1)),(YA,XYZA(2)),(ZA,XYZA(3))
      EQUIVALENCE (XB,XYZB(1)),(YB,XYZB(2)),(ZB,XYZB(3))
      EQUIVALENCE (XC,XYZC(1)),(YC,XYZC(2)),(ZC,XYZC(3))
      EQUIVALENCE (XD,XYZD(1)),(YD,XYZD(2)),(ZD,XYZD(3))
      EQUIVALENCE (PX,XYZP(1)),(PY,XYZP(2)),(PZ,XYZP(3))
      EQUIVALENCE (QX,XYZQ(1)),(QY,XYZQ(2)),(QZ,XYZQ(3))
C*
      DATA IJKL/   1,   2,   6,  10,  14,  24,  34,  44,  54,  64,  74,
     1  78,  88,  98, 108, 128, 148, 168, 188, 208, 228, 232, 242, 252,
     2 262, 282, 302, 322, 342, 362, 382, 386, 396, 406, 416, 436, 456,
     3 476, 496, 516, 536, 546, 566, 586, 606, 641, 676, 711, 746, 781,
     4 816, 826, 846, 866, 886, 921, 956, 991,1026,1061,1096,1106,1126,
     51146,1166,1201,1236,1271,1306,1341,1376,1386,1406,1426,1446,1481,
     61516,1551,1586,1621,1656,1666,1686,1706,1726,1761,1796,1831,1866,
     71901,1936,1946,1966,1986,2006,2041,2076,2111,2146,2181/
      DATA IJKLAO/1,3*2,6*3/
      DATA LBOUND/1,1,5/, UBOUND/1,4,10/
      DATA PI/3.14159265358979D0/
      DATA HALF/0.5D0/
      DATA TWOPT5/2.5D0/, THREE/3.0D0/, CUT1/499.0D0/, TENM15/1.0D-15/
C*
 1000 FORMAT('0PROGRAM UNCON CANNOT BE USED FOR F ORBITAL BASIS SETS')
 1010 FORMAT ('0FAILURE IN FMTGEN (IN UNCON) FOR T =',1PD14.6)
C*
      IF(MAXTYP.LE.2)GO TO 1
      WRITE(IOUT,1000)
      IOP(1)=-2
      RETURN
C     INITIALIZE CONSTANTS.
    1 IF(IOP(10).EQ.0)RETURN
      CONST=TWO*PI**TWOPT5
      R3=DSQRT(THREE)
      DO 6 I=1,7
    6 RA(I)=ONE
      RA( 8)=R3
      RA( 9)=R3
      RA(10)=R3
      R3OV2=HALF*R3
      DO 7 I=1,9
    7 F(I) = ZERO
      CALL FMTSET (3, 1, 0)
      IOP66=IOP(66)
      IF(IOP66.EQ.0)CALL SHLOUT(-1,TQ,1296)
      IF (IOP(1) .NE. 0) RETURN
C
      DO 405 ISHELL=1,NSHELL
      I=IXYZ(ISHELL)
      XA=C(I,1)
      YA=C(I,2)
      ZA=C(I,3)
      IA=SHELLA(ISHELL)
      IN=SHELLN(ISHELL)+IA-1
      IT=SHELLT(ISHELL)
      ISTART=LBOUND(IT+1)
C
      DO 401 JSHELL=1,ISHELL
      J=IXYZ(JSHELL)
      XB=C(J,1)
      YB=C(J,2)
      ZB=C(J,3)
      JA=SHELLA(JSHELL)
      JN=SHELLN(JSHELL)+JA-1
      JT=SHELLT(JSHELL)
      JSTART=LBOUND(JT+1)
      RABSQ=(XB-XA)*(XB-XA)+(YB-YA)*(YB-YA)+(ZB-ZA)*(ZB-ZA)
C
      DO 402 KSHELL=1,ISHELL
      K=IXYZ(KSHELL)
      XC=C(K,1)
      YC=C(K,2)
      ZC=C(K,3)
      KA=SHELLA(KSHELL)
      KN=SHELLN(KSHELL)+KA-1
      KT=SHELLT(KSHELL)
      KSTART=LBOUND(KT+1)
      MAXL=KSHELL
      IF (ISHELL.EQ.KSHELL) MAXL=JSHELL
C
      DO 403 LSHELL=1,MAXL
      LT=SHELLT(LSHELL)
C     CALCULATE THE NUMBER OF D SHELLS.
      NUMD=IT/2+JT/2+KT/2+LT/2
      IF(IOP66.NE.0.AND.NUMD.EQ.0)GO TO 403
      IMJ=ISHELL-JSHELL
      IMK=ISHELL-KSHELL
      JML=IABS(JSHELL-LSHELL)
      KML=KSHELL-LSHELL
      IMKJML=IMK+JML
      L=IXYZ(LSHELL)
      XD=C(L,1)
      YD=C(L,2)
      ZD=C(L,3)
      LA=SHELLA(LSHELL)
      LN=SHELLN(LSHELL)+LA-1
      LSTART=LBOUND(LT+1)
      RCDSQ=(XD-XC)*(XD-XC)+(YD-YC)*(YD-YC)+(ZD-ZC)*(ZD-ZC)
      IEND=UBOUND(IT+1)
      JEND=UBOUND(JT+1)
      KEND=UBOUND(KT+1)
      LEND=UBOUND(LT+1)
      IRANGE=IEND-ISTART+1
      JRANGE=JEND-JSTART+1
      KRANGE=KEND-KSTART+1
      LRANGE=LEND-LSTART+1
      MEND=IRANGE*JRANGE*KRANGE*LRANGE
      DO 730 I=1,MEND
  730 TQ(I)=ZERO
C
C     COMMENCE LOOP OVER GAUSSIAN EXPANSION
C
      INTCNT=0
      DO 300 IGAUSS=IA,IN
      AS=EXX(IGAUSS)
      CALL UFILLC(IT,IGAUSS,C1,C2,SAS)
C
      DO 301 JGAUSS=JA,JN
      BS=EXX(JGAUSS)
      CALL UFILLC(JT,JGAUSS,C1,C2,SBS)
      EP=AS+BS
      EPI=ONE/EP
      PX=(AS*XA+BS*XB)*EPI
      PY=(AS*YA+BS*YB)*EPI
      PZ=(AS*ZA+BS*ZB)*EPI
      OFAB=AS*BS*RABSQ*EPI
      IF(OFAB.GE.CUT1)GO TO 301
      OFAB=DEXP(-OFAB)
      OFAB=OFAB*CONST*EPI
C
      DO 302 KGAUSS=KA,KN
      CS=EXX(KGAUSS)
      CALL UFILLC(KT,KGAUSS,C1,C2,SCS)
C
      DO 303 LGAUSS=LA,LN
      DS=EXX(LGAUSS)
      CALL UFILLC(LT,LGAUSS,C1,C2,SDS)
      EQ=CS+DS
      EQI=ONE/EQ
      QX=(CS*XC+DS*XD)*EQI
      QY=(CS*YC+DS*YD)*EQI
      QZ=(CS*ZC+DS*ZD)*EQI
      OFCD=CS*DS*RCDSQ*EQI
      IF(OFCD.GE.CUT1)GO TO 303
      OFCD=OFAB*DEXP(-OFCD)
      OFCD=OFCD/(EQ*DSQRT(EP+EQ))
      K=0
      DO 330 I=1,3
      SASA=SCS(I)*OFCD
      DO 330 J=1,3
      K=K+1
      SAB(K)=SAS(I)*SBS(J)
  330 SCD(K)=SASA*SDS(J)
      PQX=QX-PX
      PQY=QY-PY
      PQZ=QZ-PZ
      RPQSQ=PQX*PQX+PQY*PQY+PQZ*PQZ
C     WATCH FOR SMALL VALUES OF RPQSQ.
      IF(RPQSQ.LE.TENM15)RPQSQ=ZERO
      RPQ=DSQRT(RPQSQ)
      FG=EP*EQ/(EP+EQ)
      T=RPQSQ*FG
      FG=FG+FG
C     COMPUTE FM(T).
C     AFTER GETTING FM(T), CONVERT TO FM'(T).
      M=IT+JT+KT+LT+1
      CALL FMTGEN (F, T, M, ICK)
      IF (ICK .NE. 0) THEN
         WRITE (IOUT,1010) T
         IOP(1) = -2
         RETURN
      END IF
      IF(M.EQ.1)GO TO 281
      TEMP=ONE
      DO 280 I=2,M
      TEMP=TEMP*FG
  280 F(I)=F(I)*TEMP
  281 CALL FORMS(KT+LT)
      CALL ROTATE
      CALL FABCD(SA,XYZA,XYZP,IT)
      CALL FABCD(SB,XYZB,XYZP,JT)
      CALL FABCD(SC,XYZC,XYZQ,KT)
      CALL FABCD(SD,XYZD,XYZQ,LT)
      DO 100 K=KSTART,KEND
      KFT=IJKL(K,1)
      IF (KSHELL.EQ.LSHELL) LEND=K
      DO 100 L=LSTART,LEND
      LFT=IJKL(L,1)
      KL=IJKL(K,L)
  100 CALL FORMQ(SC(KFT),SD(LFT),FT(KL),IJKLAO(K),IJKLAO(L))
      INTC=0
      DO 200 I=ISTART,IEND
      IET=IJKL(I,1)
      IF (ISHELL.EQ.JSHELL) JEND=I
      IF (ISHELL.EQ.KSHELL.AND.JSHELL.EQ.LSHELL) KEND=I
      DO 211 J=JSTART,JEND
      JET=IJKL(J,1)
      CALL FORMP(SA(IET),SB(JET),IJKLAO(I),IJKLAO(J),KT+LT)
      DO 212 K=KSTART,KEND
      LEND=UBOUND(LT+1)
      IF (KSHELL.EQ.LSHELL) LEND=K
      IF (ISHELL.EQ.KSHELL.AND.JSHELL.EQ.LSHELL.AND.I.EQ.K) LEND=J
      KLIMIT=7-IJKLAO(K)
      DO 213 L=LSTART,LEND
      M=IJKL(K,L)-1
      LLIMIT=KLIMIT-IJKLAO(L)
      Q=ET(1)*FT(M+1)
      IF (LLIMIT.EQ.5) GO TO 160
      Q=Q+ET(2)*FT(M+2)+ET(3)*FT(M+3)+ET(4)*FT(M+4)
      IF (LLIMIT.EQ.4) GO TO 160
      Q=Q+ET(5)*FT(M+5)+ET(6)*FT(M+6)+ET(7)*FT(M+7)+ET(8)*FT(M+8)
     1 +ET(9)*FT(M+9)+ET(10)*FT(M+10)
      IF (LLIMIT.EQ.3) GO TO 160
      Q=Q+ET(11)*FT(M+11)+ET(12)*FT(M+12)+ET(13)*FT(M+13)+
     1 ET(14)*FT(M+14)+ET(15)*FT(M+15)+ET(16)*FT(M+16)+ET(18)*FT(M+18)+
     2 ET(19)*FT(M+19)+ET(20)*FT(M+20)+ET(17)*FT(M+17)
      IF (LLIMIT.EQ.2) GO TO 160
      Q=Q+ET(21)*FT(M+21)+ET(22)*FT(M+22)+ET(23)*FT(M+23)+
     1 ET(24)*FT(M+24)+ET(25)*FT(M+25)+ET(26)*FT(M+26)+ET(27)*FT(M+27)+
     2 ET(28)*FT(M+28)+ET(29)*FT(M+29)+ET(30)*FT(M+30)+ET(31)*FT(M+31)+
     3 ET(32)*FT(M+32)+ET(33)*FT(M+33)+ET(34)*FT(M+34)+ET(35)*FT(M+35)
  160 INTC=INTC+1
      TQ(INTC)=TQ(INTC)+Q
  213 CONTINUE
  212 CONTINUE
  211 CONTINUE
  200 CONTINUE
      INTCNT=INTC
  303 CONTINUE
  302 CONTINUE
  301 CONTINUE
  300 CONTINUE
C
C     END OF LOOP OVER GAUSSIANS
C
      IF(INTCNT.EQ.0)GO TO 403
      INTC=0
C     RENORMALIZE CONTRACTED GAUSSIAN INTEGRALS.
      DO 870 I=ISTART,IEND
      RI=RA(I)
      IF(IMJ.EQ.0)JEND=I
      IF(IMK+JML.EQ.0)KEND=I
      DO 870 J=JSTART,JEND
      RJ=RI*RA(J)
      DO 870 K=KSTART,KEND
      RK=RJ*RA(K)
      LEND=UBOUND(LT+1)
      IF(KML.EQ.0)LEND=K
      IF(IMK+JML+IABS(I-K).EQ.0)LEND=J
      DO 870 L=LSTART,LEND
      INTC=INTC+1
  870 TQ(INTC)=TQ(INTC)*RK*RA(L)
C     RESTORE JEND, KEND AND LEND.
      JEND=UBOUND(JT+1)
      KEND=UBOUND(KT+1)
      LEND=UBOUND(LT+1)
C     RESTORE SHELL DUPLICATES; TRANSFORM 6D TO 5D.
      CALL URD65
      IF (IOP(1) .NE. 0) RETURN
  403 CONTINUE
  402 CONTINUE
  401 CONTINUE
  405 CONTINUE
C     EMPTY LAST BUFFER.
      CALL SHLOUT(0,TQ,1296)
      RETURN
      END
      SUBROUTINE UFILLC(IT,IG,C1,C2,SAS)
C*
C     --------------
C     U OF T VERSION
C     FEBRUARY 1987
C     --------------
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NP=#NP)
C##
      PARAMETER (NP=300)
C###
      PARAMETER (ZERO=0.0D0)
C
      DIMENSION C1(NP),C2(NP),SAS(3)
C
      IF(IT-1)10,20,30
C     S SHELL.
   10 SAS(1)=C1(IG)
      SAS(2)=ZERO
      SAS(3)=ZERO
      RETURN
C     SP SHELL.
   20 SAS(1)=C1(IG)
      SAS(2)=C2(IG)
      SAS(3)=ZERO
      RETURN
C     D SHELL.
   30 SAS(1)=ZERO
      SAS(2)=ZERO
      SAS(3)=C1(IG)
      RETURN
      END
      SUBROUTINE URD65
C
C     ----------------
C     QCPE GAUSSIAN 80
C     U OF T VERSION
C     FEBRUARY 1987
C     ----------------
C
C     ROUTINE TO TRANSFORM (IN A STEP-WISE FASHION) THE 6D INTEGRALS
C     TO THE FIVE PURE D INTEGRALS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA)
C##
      PARAMETER (NA= 36)
C###
      PARAMETER (ZERO=0.0D0)
C
      COMMON /A/ IOP(99),IC1(NA),IC2(7),CD1(NA,3),CD2(4),IC3(401)
C
      COMMON/C0307A/CFILL1(612),TQ(1296),CFILL2(2215),TQOUT(1296)
      COMMON/C0307B/R3OV2,IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,
     1 LSTART,IEND,JEND,KEND,LEND,MEND,NUMD,IT,JT,KT,LT,INTC,
     2 IIND(10),JIND(10),KIND(10),DUM(34)
      COMMON /SHLCOM/ ISHELL,JSHELL,KSHELL,LSHELL,IRANGE,JRANGE,KRANGE,
     1 LRANGE
C
      DIMENSION TQNEW(1296)
C
      EQUIVALENCE (TQNEW(1),TQ(1))
C
      DATA PT5/0.5D0/
C
C     INDEXING RANGES ARE PICKED UP FROM /SHLCOM/.
C     COMPUTE RANGE PRODUCTS AND INDEXING BIASES.
C
      KSL=KRANGE*LRANGE
      JSKSL=JRANGE*KSL
      ISJ=IRANGE*JRANGE
      ISTP1=ISTART-1
      JSTP1=JSTART-1
      KSTP1=KSTART-1
      LSTP1=LSTART-1
C
C     COMPUTE INDEXING ARRAYS.
C     THESE ARRAYS GIVE THE PROPER INDEX INTO THE ARRAY AS DIMENSIONED
C     BY THE MAXIMUM RANGES IN A PARTICULAR SHELL COMBINATION.
C
      ITEMP=0
      DO 100 K=KSTART,KEND
      KIND(K)=ITEMP*LRANGE
  100 ITEMP=ITEMP+1
      ITEMP=0
      DO 110 J=JSTART,JEND
      JIND(J)=ITEMP*KSL
  110 ITEMP=ITEMP+1
      ITEMP=0
      DO 120 I=ISTART,IEND
      IIND(I)=ITEMP*JSKSL
  120 ITEMP=ITEMP+1
C
C     FILL OUT TQ TO RESTORE SHELL DUPLICATES.
C
C     IN WHAT FOLLOWS, THERE ARE 6 CASES TO BE CONSIDERED.
C
C        CASE  CONDITION          RELATION OF ISHELL, ETC.
C        ----  ---------          ------------------------
C
C          1   NONE               ALL DISTINCT, IE. NO TWO SHELLS ARE
C                                 THE SAME.
C
C          2   A                  ISHELL=JSHELL ONLY (EG. (33,21)).
C
C          3   B                  KSHELL=LSHELL ONLY (EG. (32,11)).
C
C          4   C                  ISHELL=KSHELL AND JSHELL=LSHELL ONLY
C                                 (EG. (21,21)).
C
C          5   AB                 ISHELL=JSHELL AND KSHELL=LSHELL,
C                                 (EG. (22,11)).
C
C          6   ABC, AC, BC        ALL EQUAL (EG. (22,22)).
C
C
C     CASE 1 IS THE MOST FREQUENTLY OCCURRING CASE, AND IS TESTED FOR
C     FIRST.  THE OTHER CASES ARE TESTED FOR ESSENTIALLY IN THAT
C     MANNER THAT REDUCES THE NUMBER OF IF-STATEMENTS.
C
C     IF THE PRODUCT IS NON-ZERO, WE HAVE CASE 1.
C
      IF(IMJ*KML*IMKJML.NE.0)GO TO 1000
C
C     ENTRY AT THIS POINT IMPLIES THAT AT LEAST ONE OF THE 3
C     POSSIBLE SHELL COINCIDENCE FLAGS IS ZERO.
C     THE ALL EQUAL CASE IS EASIEST TO DETECT.
C
      IF(IMJ+IMKJML.EQ.0)GO TO 6
C
C     ENTRY AT THIS POINT IMPLIES THAT ONE OR BOTH OF IMJ, IMKJML
C     IS NON-ZERO.
C     THEREFORE, IMKJML=0 IMPLIES CASE 4.
C
      IF(IMKJML.EQ.0)GO TO 4
C
C     AT THIS POINT, EITHER IMJ OR KML (OR BOTH) IS (ARE) ZERO.
C     TEST FOR BOTH SIMULTANEOUSLY ZERO.
C
      IF(IMJ+KML.EQ.0)GO TO 5
C
C     AT THIS POINT, BOTH ARE NOT ZERO TOGETHER.
C     (IMJ*IMKJML*KML)=0 AND THE PREVIOUS TWO TESTS  MEAN THAT AT
C     LEAST ONE IS ZERO.
C     THEREFORE, ONLY TEST FOR ONE OF THEM (IMJ).
C
      IF(IMJ.NE.0)GO TO 3
C
C     ALL POSSIBLE CASES HAVE BEEN ELIMINATED.
C
C     ******************************************************************
C     CASE 1, NO OPERATION REQUIRED, INTEGRALS ARE READY FOR TRANS-
C     FORMATION.
C     ******************************************************************
C
C     ******************************************************************
C     CASE 2, ISHELL=JSHELL.
C     ******************************************************************
C
C     TQ CONTAINS (IRANGE*(IRANGE+1))/2 RECTANGULAR MATRICES
C     WITH EACH CONTAINING (LRANGE*KRANGE) INTEGRALS.
C     WE DESIRE TO EXPAND (IJ) THROUGHOUT INTO TQNEW.
C     TO DO THIS, WE FALL BACKWARDS THROUGH TQ.
C     THE PLOY IS AS FOLLOWS,
C        SELECT A PAIR (IJ) IN DECREMENTING MODE.
C        THERE ARE TWO POSSIBILITIES,
C             1. I.GE.J, HERE, THE DESIRED RECTANGULAR ARRAY
C             IS LOCATED IN TQ.
C             IT IS COPIED TO TQNEW IN DESCENDING MODE, USING RUNNING
C             INDICES.
C             2. I.LT.J, HERE, THE DESIRED ARRAY IS ALREADY IN TQNEW.
C             IT IS COPIED, AGAIN IN DESCENDING MODE, TO ITS NEW
C             RESTING PLACE IN TQNEW.
C
C     NOTE THAT IF KSHELL=LSHELL, THIS CASE DEGENERATES TO A SERIES
C     OF SINGLE ELEMENT COPIES.
C
C     BYPASS EXPANSION IF IRANGE=1.
C     SINCE ISHELL=JSHELL, IT IS ONLY NECESSARY TO TEST IRANGE.
C
      IF(IRANGE.EQ.1)GO TO 1000
C
C     PERFORM INITIALIZATION.
C     LWATQ IS PICKED UP FROM INTC   (INDX2).
C     LWATQN IS PICKED UP FROM MEND   (INDX1).
C     THESE ARE MAINTAINED AS RUNNING INDICES.
C     THUS, IN ANY GIVEN PASS, THE ONLY INDEX THAT MUST BE
C     COMPUTED  IS INDX3 (BASED ON IND AND JND).
C
      INDX1=MEND+1
      INDX2=INTC+1
      IND=IEND+1
      DO 350 I=ISTART,IEND
      IND=IND-1
      JND=JEND+1
      DO 350 J=JSTART,JEND
      JND=JND-1
C     NOW HAVE A PAIR (IND,JND).  MAKE TESTS AND BRANCH TO
C     PARTICULAR COPY CODE.
      IF(IND.LT.JND)GO TO 330
C     PREFERRED CASE, IND.GE.JND.
C     HERE, WE COPY FROM TQ TO TQNEW SEQUENTIALLY BACKWARDS.
      DO 320 KL=1,KSL
      INDX1=INDX1-1
      INDX2=INDX2-1
  320 TQNEW(INDX1)=TQ(INDX2)
      GO TO 350
C     JND.GT.IND.
C     HERE, COPY FROM TQNEW TO TQNEW, USING COMPUTED INDEX INDX3.
  330 INDX3=IIND(JND)+JIND(IND)+KSL
      DO 340 KL=1,KSL
      INDX1=INDX1-1
      TQNEW(INDX1)=TQNEW(INDX3)
  340 INDX3=INDX3-1
  350 CONTINUE
      GO TO 1000
C
C     ******************************************************************
C     CASE 3, KSHELL=LSHELL ONLY.
C     ******************************************************************
C
C     CURRENTLY, TQ CONTAINS IRANGE*JRANGE SYMMETRIC MATRICES
C     ((KRANGE*(KRANGE+1))/2 ELEMENTS IN EACH).
C     WE DESIRE TO EXPAND K.GE.L FOR EACH PAIR (IJ) THROUGHOUT.
C     THUS, TQNEW WILL CONTAIN IRANGE*JRANGE SQUARE MATRICES.
C
C     THE TACTIC HERE IS TO FALL BACKWARDS THROUGH ALL PAIRS (IJ)
C     AND TO PERFORM A SIMPLE LINEAR TO SQUARE CONVERSION
C     FOR EACH.
C     THIS PROCESS IS ACCOMPLISHED IN 2 STAGES FOR EACH PAIR (IJ).
C     FIRST, THE EXISTING ELEMENTS OF THE SYMMETRIC MATRIX ARE
C     TRANSFERRED TO THE APPROPRIATE PLACES IN TQNEW.
C     SECOND, THE MATRIX IS EXPANDED TO SQUARE FORM IN PLACE IN
C     TQNEW.
C
C     IN ANY PASS THROUGH THE IJ-LOOP, THE LWATQN FOR THE NTT WORD
C     TRANSFER IS COMPUTED FROM INDX1, A RUNNING INDEX STARTING AT MEND.
C     SIMILARLY, LWATQ IS MAINTAINED IN INDX2, A RUNNING INDEX THAT
C     STARTS AT INTC.
C
C     THIS UNPACK IS BYPASSED IF KRANGE AND LRANGE=1.
C     SINCE KSHELL=LSHELL, IT FOLLOWS THAT KRANGE=LRANGE, AND
C     IT IS SUFFICIENT TO TEST JUST KRANGE.
C
    3 IF(KRANGE.EQ.1)GO TO 1000
C     PERFORM NECESSARY INITIALIZATION.
      NTT=(KRANGE*(KRANGE+1))/2
      KRP1=KRANGE+1
      KRM1=KRANGE-1
      LWATQ=INTC
      LWATQN=MEND
      DO 450 IJ=1,ISJ
      INDX1=LWATQN
      INDX2=LWATQ
      LWATQ=LWATQ-NTT
      LWATQN=LWATQN-KSL
C     PERFORM NTT-WORD TRANSFER.
      DO 420 K=1,KRANGE
      LLIM=KRP1-K
      DO 410 L=1,LLIM
      TQNEW(INDX1)=TQ(INDX2)
      INDX2=INDX2-1
  410 INDX1=INDX1-1
  420 INDX1=INDX1-K
C     PERFORM EXPANSION OF (KL) IN PLACE IN TQNEW.
      INDX1=LWATQN+2
      INDX2=LWATQN+KRP1
      DO 440 K=1,KRM1
C     USE INDS1 AND INDS2 TO PRESERVE INDX1 AND INDX2.
      INDS1=INDX1
      INDS2=INDX2
      KRMK=KRANGE-K
      DO 430 L=1,KRMK
      TQNEW(INDX1)=TQNEW(INDX2)
      INDX1=INDX1+1
  430 INDX2=INDX2+KRANGE
      INDX1=INDS1+KRP1
  440 INDX2=INDS2+KRP1
  450 CONTINUE
      GO TO 1000
C
C     ******************************************************************
C     CASE 4, ISHELL=KSHELL AND JSHELL=LSHELL.
C     ******************************************************************
C
C     THIS CASE (AND ALSO CASE 6) IS SOMEWHAT COMPLICATED AND A
C     SUBSTANTIAL AMOUNT OF OVERHEAD IS INCURRED.
C
C     THE PROCEDURE IS AS FOLLOWS.
C
C     WE STEP BACKWARDS THROUGH (I,J,K,L), USING THE FULL POSSIBLE
C     RANGE.
C     INSIDE THE L-LOOP, THE FOUR INDICES ARE EXAMINED.  THE FOLLOWING
C     LIMITS APPLY TO THE INTEGRALS IN TQ,
C
C        I.GE.K
C        WHEN I=K, J.GE.L
C
C     THE FOUR INDICES ARE TESTED AGAINST THESE CONDITIONS AND
C     EITHER
C
C        1.  THE CONDITIONS ARE MET, IN WHICH CASE, THE NEXT
C            SEQUENTIALLY DECREMENTING INTEGRAL IS COPIED FROM TQ TO
C            TQNEW AND THE APPROPRIATE COUNTERS (INDICES) ARE
C            DECREMENTED.
C
C        2.  THE CONDITIONS ARE NOT SATISFIED.  IN THIS CASE, DUE TO THE
C            NATURE OF THE COPY, THE DESIRED INTEGRAL ALREADY RESIDES
C            IN TQNEW.
C            IT IS TRANSFERRED, AND AGAIN
C            INDEXING IS DONE.
C
C     INDX1 INDEXES IN TQNEW, AND STARTS AT MEND.
C     INDX2 INDEXES IN TQ, AND STARTS AT INTC.
C     INDX3 INDEXES IN TQNEW WHEN COPYING FORM TQNEW TO TQNEW.  THIS
C     INDEX MUST BE COMPUTED.
C
    4 INDX1=MEND
      INDX2=INTC+1
      IND=IEND+1
      DO 560 I=ISTART,IEND
      IND=IND-1
      ITEMP=KIND(IND)-JSTP1
      JND=JEND+1
      DO 560 J=JSTART,JEND
      JND=JND-1
      JTEMP=ITEMP+JND
      KND=KEND+1
      DO 560 K=KSTART,KEND
      KND=KND-1
      KTEMP=JTEMP+IIND(KND)
      LND=LEND+1
      DO 560 L=LSTART,LEND
      LND=LND-1
C     NOW IN L-LOOP, PERFORM TESTS.
      IF(IND-KND)530,520,550
  520 IF(JND.GE.LND)GO TO 550
C     COPY FROM TQNEW TO TQNEW AFTER COMPUTING INDX3.
  530 INDX3=KTEMP+JIND(LND)
      TQNEW(INDX1)=TQNEW(INDX3)
      GO TO 560
C     COPY NEXT SEQUENTIAL INTEGRAL FROM TQ TO TQNEW.
  550 INDX2=INDX2-1
      TQNEW(INDX1)=TQ(INDX2)
  560 INDX1=INDX1-1
      GO TO 1000
C
C     ******************************************************************
C     CASE 5, ISHELL=JSHELL AND KSHELL=LSHELL.
C     ******************************************************************
C
C     THIS CASE IS ESSENTIALLY THE SAME AS CASES 2 AND 3 COMBINED.
C     LWA IN TQ IS MOST CONVENIENTLY SPECIFIED BY INTC (MAINTAINED IN
C     INDX2).
C     LWA IN TQNEW IS IIND(IND)+JIND(JND)+KSL, OR STARTING AT
C     MEND, IT IS DECREMENTED BY KSL ON EACH PASS.
C
C     FOR EACH PAIR (IND,JND) WE EITHER
C        (IND.GE.JND) COPY FROM TQ TO TQNEW AND THEN EXPAND (KL) IN
C                     PLACE IN TQNEW.
C        (IND.LT.JND) COPY (STILL DECREMENTING) DIRECTLY FROM TQNEW
C                     TO TQNEW.
C
C     NOTE THAT IF  KRANGE AND LRANGE ARE 1, AN ALTERNATIVE PROCEDURE
C     MUST BE USED .
C
    5 IF(KRANGE.EQ.1)GO TO 710
      LWATQ=INTC
      LWATQN=MEND
      NTT=(KRANGE*(KRANGE+1))/2
      KRP1=KRANGE+1
      KRM1=KRANGE-1
      IND=IEND+1
      DO 690 I=ISTART,IEND
      IND=IND-1
      JND=JEND+1
      DO 690 J=JSTART,JEND
      JND=JND-1
C     TEST (IND,JND) TO DETERMINE THE ACTION TO BE TAKEN.
      IF(IND.LT.JND)GO TO 670
C     IND.GE.JND, DO AS IN CASE 3.
C     COPY NTT INTEGRALS FROM TQ TO TQNEW WITH CORRECT PLACEMENT.
      INDX2=LWATQ
      INDX1=LWATQN
      LWATQ=LWATQ-NTT
      LWATQN=LWATQN-KSL
      DO 640 K=1,KRANGE
      LLIM=KRP1-K
      DO 630 L=1,LLIM
      TQNEW(INDX1)=TQ(INDX2)
      INDX2=INDX2-1
  630 INDX1=INDX1-1
  640 INDX1=INDX1-K
C     EXECUTE EXPANSION OF (KL).
      INDX1=LWATQN+2
      INDX2=LWATQN+KRP1
      DO 660 K=1,KRM1
      INDS1=INDX1
      INDS2=INDX2
      KRMK=KRANGE-K
      DO 650 L=1,KRMK
      TQNEW(INDX1)=TQNEW(INDX2)
      INDX1=INDX1+1
  650 INDX2=INDX2+KRANGE
      INDX1=INDS1+KRP1
  660 INDX2=INDS2+KRP1
      GO TO 690
C
C     (IND.LT.JND), DESIRED ARRAY ALREADY RESIDES IN TQNEW.
C                   FIND IT AND COPY IT OUT.
C     THIS IS DONE BY COPY SEQUENTIALLY IN REVERSE FROM TQNEW TO
C     TQNEW.
C     INDX1 STEPS INPUT INTO TQNEW.
C     INDX2 STEPS OUTPUT FROM TQNEW.
C     DETERMINE INDX1 FROM LWATQN, AND DECREMENT LWATQN.
C     LWATQN (OUTPUT) IS COMPUTED FROM THE STANDARD INDEXING ARRAYS.
C
  670 INDX1=LWATQN
      LWATQN=LWATQN-KSL
      INDX2=IIND(JND)+JIND(IND)+KSL
C     COPY OVER KL.
      DO 680 KL=1,KSL
      TQNEW(INDX1)=TQNEW(INDX2)
      INDX1=INDX1-1
  680 INDX2=INDX2-1
  690 CONTINUE
      GO TO 1000
C
C     IN CASE 5 (ISHELL=JSHELL AND KSHELL=LSHELL), IT IS NECESSARY
C     TO PROCEED DIFFERENTLY IF KRANGE AND LRANGE EQUAL ONE.
C     IN THIS CASE, WE HAVE WHAT AMOUNTS TO A SYMMETRIC
C     MATRIX IN (IJ).  THE DIMENSION IS EITHER 1, 4 OR 6.
C     IN ANY EVENT, WE MERELY PERFORM A LINEAR TO SQUARE CONVERSION.
C
C     THE LWA IN TQ IS CLEARLY SPECIFIED BY INTC.   (INDX2)
C     THE LWA IN TQNEW IS MEND.   (INDX1)
C
C     IF IRANGE (AND THEREFORE JRANGE) IS ALSO 1, SKIP THE EXPANSION.
C
  710 IF(IRANGE.EQ.1)GO TO 1000
      INDX1=MEND
      INDX2=INTC
      IRP1=IRANGE+1
      IRM1=IRANGE-1
C     ARRANGE THE NTT EXISTING ELEMENTS.
      DO 730 I=1,IRANGE
      JLIM=IRP1-I
      DO 720 J=1,JLIM
      TQNEW(INDX1)=TQ(INDX2)
      INDX2=INDX2-1
  720 INDX1=INDX1-1
  730 INDX1=INDX1-I
C     EXPAND OVER (IJ).
C     THE FWA IS ASSUMED TO BE 1.
      INDX1=2
      INDX2=IRP1
      DO 750 I=1,IRM1
      INDS1=INDX1
      INDS2=INDX2
      IRMI=IRANGE-I
      DO 740 J=1,IRMI
      TQNEW(INDX1)=TQNEW(INDX2)
      INDX1=INDX1+1
  740 INDX2=INDX2+JRANGE
      INDX1=INDS1+IRP1
  750 INDX2=INDS2+IRP1
      GO TO 1000
C
C     ******************************************************************
C     CASE 6, ALL SHELL INDICES ARE EQUAL.
C     ******************************************************************
C
C     THE TECHNIQUES OF CASE 4 ARE EMPLOYED.
C     THIS CASE INVOLVES THE MOST OVERHEAD, BUT IS THE LEAST FRE-
C     QUENTLY EXECUTED PART (LESS THAN OR EQUAL TO NSHELL TIMES
C     PER INTEGRAL EVALUATION).
C
    6 INDX1=MEND
      INDX2=INTC+1
      IND=IEND+1
      DO 890 I=ISTART,IEND
      IND=IND-1
      JND=JEND+1
      DO 890 J=JSTART,JEND
      JND=JND-1
      KND=KEND+1
      DO 890 K=KSTART,KEND
      KND=KND-1
      LND=LEND+1
      DO 890 L=LSTART,LEND
      LND=LND-1
C     NOW HAVE ALL FOUR INDICES (IND,JND,KND,LND).
C     IFLAG IS FOR DETERMINING WHERE TO PULL THE NEXT INTEGRAL FROM.
C     IFLAG=0, INTEGRAL COMES FROM TQNEW.
C     IFLAG=1, INTEGRAL COMES FROM TQ.
C     INDX2 STEPS BACKWARDS THROUGH TQ.  IT IS DECREMENTED ONLY WHEN
C     AN INTEGRAL IS ACTUALLY COPIED.
C     INDX1 STEPS INPUT INTO TQNEW.  IT IS ALWAYS DECREMENTED.
C     INDX3 STEPS OUTPUT FROM TQNEW INTO TQNEW.
C     THIS INDEX IS COMPUTED ON DEMAND.
C     NOTE THAT IN THE FOLLOWING, IND, ETC. ARE SORTED.
      IFLAG=1
      INEW=IND
      JNEW=JND
      KNEW=KND
      LNEW=LND
      IF(INEW.GE.JNEW)GO TO 810
C     SWITCH (IJ), AND SET IFLAG.
      ITEMP=INEW
      INEW=JNEW
      JNEW=ITEMP
      IFLAG=0
C     SWITCH (KL) AND SET IFLAG.
  810 IF(KNEW.GE.LNEW)GO TO 830
      ITEMP=KNEW
      KNEW=LNEW
      LNEW=ITEMP
      IFLAG=0
  830 IF(INEW-KNEW)850,840,860
  840 IF(JNEW.GE.LNEW)GO TO 860
C     SWITCH (IJ) AND (KL) AND SET IFLAG.
  850 ITEMP=INEW
      INEW=KNEW
      KNEW=ITEMP
      ITEMP=JNEW
      JNEW=LNEW
      LNEW=ITEMP
      GO TO 880
  860 IF(IFLAG.EQ.0)GO TO 880
C     COPY FROM TQ TO TQNEW.
      INDX2=INDX2-1
      TQNEW(INDX1)=TQ(INDX2)
      GO TO 890
C     COPY FROM TQNEW TO TQNEW.
  880 INDX3=IIND(INEW)+JIND(JNEW)+KIND(KNEW)+LNEW-LSTP1
      TQNEW(INDX1)=TQ(INDX3)
C     ALWAYS DECREMENT INDX1.
  890 INDX1=INDX1-1
C
C***********************************************************************
C     TRANSFORMATION SECTION.
C***********************************************************************
C
C     WE ARE FINALLY READY TO DO, IN A STEPWISE FASHION, THE 6D TO 5D
C     CONVERSION.  AT THIS POINT, IT IS APPROPRIATE TO CLEAN UP THE
C     NOTATION.
C
C     THE LOOP ORDER IS (OUTERMOST) I, J, K, L (INNERMOST).
C
C     THE I-LOOP GOES OVER THE FUNCTIONS AT CENTER A.
C     THE J-LOOP GOES OVER THE FUNCTIONS AT CENTER B.
C     THE K-LOOP GOES OVER THE FUNCTIONS AT CENTER C.
C     THE L-LOOP GOES OVER THE FUNCTIONS AT CENTER D.
C
C     THE TRANSFORMATION IS DONE STEP-WISE, USING THE FOLLOWING
C     TRANSFORMATION MATRIX,
C
C        ( 1  0  0  0  0  0  0  0  0  0 ) ( S    )   ( S           )
C        (                              ) (      )   (             )
C        ( 0  1  0  0  0  0  0  0  0  0 ) ( X    )   ( X           )
C        (                              ) (      )   (             )
C        ( 0  0  1  0  0  0  0  0  0  0 ) ( Y    )   ( Y           )
C        (                              ) (      )   (             )
C        ( 0  0  0  1  0  0  0  0  0  0 ) ( Z    )   ( Z           )
C        (                              ) (      )   (             )
C        ( 0  0  0  0 -H -H  1  0  0  0 ) ( X**2 ) = ( 3*Z**2-R**2 )
C        (                              ) (      )   (             )
C        ( 0  0  0  0  R -R  0  0  0  0 ) ( Y**2 )   ( X**2-Y**2   )
C        (                              ) (      )   (             )
C        ( 0  0  0  0  0  0  0  1  0  0 ) ( Z**2 )   ( XY          )
C        (                              ) (      )   (             )
C        ( 0  0  0  0  0  0  0  0  1  0 ) ( XY   )   ( XZ          )
C        (                              ) (      )   (             )
C        ( 0  0  0  0  0  0  0  0  0  1 ) ( XZ   )   ( YZ          )
C                                         (      )
C                                         ( YZ   )
C
C     WHERE H=0.5, AND R=SQRT(3.0)/2.0.
C     SINCE THIS TRANSFORMATION IS CLOSE TO AN IDENTITY TRANSFORMATION,
C     ONLY THE REQUIRED STEPS ARE ACTUALLY CARRIED OUT.
C
C     THE NUMBER OF WORDS IN TQNEW IS EQUAL TO (IRANGE*JRANGE*KRANGE*
C     LRANGE).
C
C     ALL FURTHER OPERATIONS TAKE PLACE ENTIRELY IN TQNEW.
C     THE CURRENT ORDER OF FUNCTIONS IS
C
C         1,2,3,4,5,   6,   7,   8, 9, 10
C        (S,X,Y,Z,X**2,Y**2,Z**2,XY,XZ,YZ)
C
C     INITIALIZATION SECTION.
 1000 IEND5D=IEND
      JEND5D=JEND
      KEND5D=KEND
      LEND5D=LEND
      IF(IOP(8).EQ.0)GO TO 1010
C     TRANSFORM 6D TO 5D.
      IF(IT.EQ.2)IEND5D=9
      IF(JT.EQ.2)JEND5D=9
      IF(KT.EQ.2)KEND5D=9
      IF(LT.EQ.2)LEND5D=9
 1010 IRNG5D=IEND5D-ISTP1
      JRNG5D=JEND5D-JSTP1
      KRNG5D=KEND5D-KSTP1
      LRNG5D=LEND5D-LSTP1
      K5SL5=KRNG5D*LRNG5D
      J6K5L5=JRANGE*K5SL5
      J5K5L5=JRNG5D*K5SL5
C     BY-PASS 6D TO 5D TRANSFORMATION IF POSSIBLE.
      IF(IOP(8).EQ.0.OR.NUMD.EQ.0)GO TO 1600
C
C     BYPASS TRANSFORMATION AT CENTERS C AND D IF THERE ARE NO SECOND
C     ORDER FUNCTIONS AT EITHER CENTER.
      IF((KT-2)*(LT-2).NE.0)GO TO 1400
C
C     EXECUTE TRANSFORMATION AT CENTERS C AND D.
C
C     DESCRIPTION OF INDEXING, TRANSFORMATION AT CENTER D (L-INDEX).
C
C        IIND(I)+JIND(J)+KIND(K)
C
C     GIVES THE (FWA-1) OF A VECTOR (OVER L) THAT IS LRANGE WORDS LONG.
C     X**2 CORRESPONDS TO L=5.  THUS, LNEW=5-LSTP1.  HENCE,
C
C        IIND(I)+JIND(J)+KIND(K)+5-LSTP1
C
C     POINTS TO AN X**2 FUNCTION ON CENTER D.  SINCE WE START AT THE
C     BEGINNING, THE INITIAL INDEX IS 5-LSTP1 AND IS SUBSEQUENTLY IN-
C     CREMENTED BY LRANGE.
C
C     DESCRIPTION OF INDEXING, TRANSFORMATION AT CENTER C (K-INDEX).
C
C     THE ARGUMENT IS SIMILAR TO THE ABOVE,
C     INDEXING STARTS AT
C
C        IIND(I)+JIND(J)+KIND(5)+1
C
C     AND IS INCREMENTED BY 1 ON EACH PASS THROUGH THE L-LOOP.
C
C     INDX1 IS USED FOR THE L-TRANSFORMATION.
C     INDX2 IS USED FOR THE K-TRANSFORMATION.
C
C     INITIALIZE INDX1.
      INDX1=5-LSTP1
C     COLLECT FACTORS OF LRANGE.
      LRANG2=LRANGE+LRANGE
      LRANG3=LRANG2+LRANGE
      LRANG4=LRANG3+LRANGE
      LRANG5=LRANG4+LRANGE
      K5P1=KIND(5)+1
C
C     COMMENCE LOOPS OVER THE FUNCTIONS AT CENTERS A AND B.
C
      DO 1140 I=ISTART,IEND
      ITEMP=IIND(I)+K5P1
      DO 1140 J=JSTART,JEND
C     HAVE ONE PAIR (IJ).
C     BYPASS TRANSFORMATION AT D IF POSSIBLE.
      IF(LT.NE.2)GO TO 1120
C     LOOP OVER ALL POSSIBLE VALUES OF K (FUNCTIONS AT CENTER C).
      DO 1110 K=KSTART,KEND
C     HAVE ONE TRIPLE (IJK).
C     DO TRANSFORMATION AT CENTER D FOR THIS TRIPLE.
C     SAVE X**2 AT D.
      TEMP=TQNEW(INDX1)
C     COMPUTE 3*Z**2-R**2 AT D.
      TQNEW(INDX1)=TQNEW(INDX1+2)-PT5*(TEMP+TQNEW(INDX1+1))
C     COMPUTE X**2-Y**2 AT D.
      TEMP=R3OV2*(TEMP-TQNEW(INDX1+1))
C     SHIFT REMAINING FUNCTIONS AT D FOR THIS TRIPLE (IJK).
      TQNEW(INDX1+1)=TQNEW(INDX1+4)
      TQNEW(INDX1+2)=TQNEW(INDX1+5)
      TQNEW(INDX1+4)=TQNEW(INDX1+3)
      TQNEW(INDX1+3)=TEMP
 1110 INDX1=INDX1+LRANGE
C     TRANSFORMATION COMPLETE AT CENTER D FOR THE PAIR (IJ) AND ONE
C     RANGE OF K.
C     THE FUNCTIONS AT D HAVE THE ORDER
C        (S,X,Y,Z,3*Z**2-R**2,XZ,YZ,X**2-Y**2,XY)
C     THE TENTH SLOT IS TO BE CONSIDERED BLANK. (IT IS DONE AWAY
C     WITH LATER).
C
C     BYPASS TRANSFORMATION AT C IF POSSIBLE.
      IF(KT.NE.2)GO TO 1140
C     PERFORM TRANSFORMATION AT C.
C     OBTAIN STARTING INDEX.
 1120 INDX2=JIND(J)+ITEMP
C     LOOP OVER FUNCTIONS AT CENTER D.  HERE, WE NEED ONLY GO UP TO
C     LEND5D.
      DO 1130 L=LSTART,LEND5D
C     HAVE ONE TRIPLE (IJL).
C     DO TRANSFORMATION AT CENTER C.
      TEMP=TQ(INDX2)
      TQNEW(INDX2)=TQNEW(INDX2+LRANG2)-PT5*(TEMP+TQNEW(INDX2+LRANGE))
      TEMP=R3OV2*(TEMP-TQNEW(INDX2+LRANGE))
      TQNEW(INDX2+LRANGE)=TQNEW(INDX2+LRANG4)
      TQNEW(INDX2+LRANG2)=TQNEW(INDX2+LRANG5)
      TQNEW(INDX2+LRANG4)=TQNEW(INDX2+LRANG3)
      TQNEW(INDX2+LRANG3)=TEMP
 1130 INDX2=INDX2+1
 1140 CONTINUE
C     THE ORDER  OF THE FUNCTIONS AT A AND B IS UNCHANGED.
C     THE ORDER OF THE FUNCTIONS AT C AND D IS NOW
C        (S,X,Y,Z,3*Z**2-R**2,XZ,YZ,X**2-Y**2,XY).
C     DUE TO THESE ALTERATIONS, THERE ARE A NUMBER OF GAPS IN TQNEW.
C     THESE GAPS ARE COMPRESSED BEFORE PROCEEDING.
C
C     COMPUTE AUXILLIARY INCREMENTS FOR COMPRESSION.
      LINCR=1
      IF(LT.NE.2)LINCR=0
      KINCR=LRANGE
      IF(KT.NE.2)KINCR=0
C
C     INITIALIZE INDICES.  INDX1 IS FOR COMPRESSED ARRAY.
C     INDX2 COUNTS IN FULL ARRAY.
      INDX1=0
      INDX2=0
C     LOOP OVER COMBINED FUNCTIONS AT A AND B.
      DO 1210 IJ=1,ISJ
C     LOOP OVER 5D FUNCTIONS AT CENTER C.
      DO 1200 K=KSTART,KEND5D
C     LOOP OVER 5D FUNCTIONS AT CENTER D.
      DO 1190 L=LSTART,LEND5D
C     INCREMENT COUNTERS AND COPY.
      INDX1=INDX1+1
      INDX2=INDX2+1
 1190 TQNEW(INDX1)=TQNEW(INDX2)
C     TERMINATE LOOP AT CENTER C BY INCREMENTING INDX2 BY 1 (LT=2)
C     OR 0 (LT.NE.2).
C     THIS STEPS PAST THE EXTRANEOUS FUNCTION AT D.
 1200 INDX2=INDX2+LINCR
C     TERMINATE THE IJ LOOP BY INCREMENTING INDX2 BY LRANGE (KT=2) OR
C     0 (KT.NE.2).
C     THIS STEPS PAST THE FUNCTIONS CORRESPONDING TO K=10 THAT ARE
C     DUPLICATE.
 1210 INDX2=INDX2+KINCR
C     RE-COMPUTE IIND, JIND, AND KIND.
C     BYPASS THE RE-COMPUTATION OF KIND IF POSSIBLE.
      IF(LT.NE.2)GO TO 1250
      ITEMP=0
      DO 1240 K=KSTART,KEND5D
      KIND(K)=ITEMP*LRNG5D
 1240 ITEMP=ITEMP+1
 1250 ITEMP=0
      DO 1260 J=JSTART,JEND
      JIND(J)=ITEMP*K5SL5
 1260 ITEMP=ITEMP+1
      ITEMP=0
      DO 1270 I=ISTART,IEND
      IIND(I)=ITEMP*J6K5L5
 1270 ITEMP=ITEMP+1
C
C     TRANSFORMATION COMPLETE AT CENTERS C AND D.
C     THE ARRAY TQNEW CONTAINS IRANGE*JRANGE*KRNG5D*LRNG5D INTEGRALS.
C     THE ONLY DUPLICATE INTEGRALS CORRESPOND TO SHELL DUPLICATES.
C
C     TRANSFORM THE FUNCTIONS AT CENTER B.
C     HERE WE LOOP OVER ALL I, AND THE PRODUCT OF THE RANGES AT K AND L.
C
C     FORM FACTORS OF KSL.
 1400 K5SL52=K5SL5+K5SL5
      K5SL53=K5SL52+K5SL5
      K5SL54=K5SL53+K5SL5
      K5SL55=K5SL54+K5SL5
C
C     BYAPSS TRANSFORMATION AT CENTER B IF POSSIBLE.
      IF(JT.NE.2)GO TO 1500
C     PERFORM TRANSFORMATION AT CENTER B.
      DO 1430 I=ISTART,IEND
      INDX1=IIND(I)+JIND(5)
      DO 1430 KL=1,K5SL5
C     HAVE ONE PSEUDO-TRIPLE (IKL).  PERFORM TRANSFORMATION AT
C     CENTER B.
      INDX1=INDX1+1
      TEMP=TQNEW(INDX1)
      TQNEW(INDX1)=TQNEW(INDX1+K5SL52)-PT5*(TEMP+TQNEW(INDX1+K5SL5))
      TEMP=R3OV2*(TEMP-TQNEW(INDX1+K5SL5))
      TQNEW(INDX1+K5SL5)=TQNEW(INDX1+K5SL54)
      TQNEW(INDX1+K5SL52)=TQNEW(INDX1+K5SL55)
      TQNEW(INDX1+K5SL54)=TQNEW(INDX1+K5SL53)
 1430 TQNEW(INDX1+K5SL53)=TEMP
C     TRANSFORMATION AT CENTER B IS COMPLETE.
C
C     PASS THROUGH TQNEW AS BEFORE, AND COMPRESS OUT UNNECESSARY
C     INTEGRALS.
C     FOR EACH I, THERE IS A BLOCK (K5SL5 WORDS LONG) AFTER J=9
C     THAT MUST BE COMPRESSED OUT.
C     INDX1 COUNTS IN THE COMPRESSED ARRAY, INDX2 COUNTS IN THE FULL
C     ARRAY.
      INDX1=0
      INDX2=0
C     LOOP OVER FUNCTIONS AT CENTER A.
      DO 1460 I=ISTART,IEND
C     LOOP OVER COMBINED FUNCTIONS AT KL.
      DO 1450 KL=1,J5K5L5
      INDX1=INDX1+1
      INDX2=INDX2+1
 1450 TQNEW(INDX1)=TQNEW(INDX2)
C     INCREMENT COUNTER AN ADDITIONAL K5SL5 AT END OF EACH I-LOOP.
 1460 INDX2=INDX2+K5SL5
C     RE-COMPUTE IIND.
      ITEMP=0
      DO 1470 I=ISTART,IEND
      IIND(I)=ITEMP*J5K5L5
 1470 ITEMP=ITEMP+1
C
C     TRANSFORMATION AND COMPRESSION COMPLETE AT CENTERS B, C, AND D.
C     THE ORDER OF THE FUNCTIONS AT A IS UNCHANGED.
C     THE ORDER OF THE FUNCTIONS AT CENTERS B, C, AND D IS
C        (S,X,Y,Z,3*Z**2-R**2,XZ,YZ,X**2-Y**2,XY)
C     THERE ARE CURRENTLY IRANGE*JRNG5D*KRNG5D*LRNG5D INTEGRALS
C     IN THE ARRAY TQNEW.
C     THE ONLY REDUNDANCIES CORRESPOND TO SHELL DUPLICATES.
C
C     TRANSFORM THE FUNCTIONS AT CENTER A.
C     SKIP THE TRANFORMATION AT CENTER A IF POSSIBLE.
 1500 IF(IT.NE.2)GO TO 1600
C     INITIALIZE INDX1.
      INDX1=IIND(5)
C     FORM FACTORS OF J5K5L5.
      JKL52=J5K5L5+J5K5L5
      JKL53=JKL52+J5K5L5
      JKL54=JKL53+J5K5L5
      JKL55=JKL54+J5K5L5
C     PASS OVER COMBINED FUNCTIONS AT B, C, AND D.
      DO 1520 JKL=1,J5K5L5
      INDX1=INDX1+1
      TEMP=TQNEW(INDX1)
      TQNEW(INDX1)=TQNEW(INDX1+JKL52)-PT5*(TEMP+TQNEW(INDX1+J5K5L5))
      TEMP=R3OV2*(TEMP-TQNEW(INDX1+J5K5L5))
      TQNEW(INDX1+J5K5L5)=TQNEW(INDX1+JKL54)
      TQNEW(INDX1+JKL52)=TQNEW(INDX1+JKL55)
      TQNEW(INDX1+JKL54)=TQNEW(INDX1+JKL53)
 1520 TQNEW(INDX1+JKL53)=TEMP
C
C     TRANSFORMATION AT CENTER A IS COMPLETE.
C     THE ORDERING AT ALL CENTERS IS
C
C        (S,X,Y,Z,3*Z**2-R**2,XZ,YZ,X**2-Y**2,XY).
C
C     NO COMPRESSION IS REQUIRED AFTER THE TRANSFORMATION AT
C     CENTER A.
C     TQNEW CURRENTLY CONTAINS IRNG5D*JRNG5D*KRNG5D*LRNG5D INTEGRALS.
C*
C     UPDATE IRANGE TO LRANGE IN /SHLCOM/.
 1600 IRANGE=IRNG5D
      JRANGE=JRNG5D
      KRANGE=KRNG5D
      LRANGE=LRNG5D
      MEND=IRANGE*JRANGE*KRANGE*LRANGE
C     IT IS FINALLY TIME TO GET RID OF THE SHELL DUPLICATES.
C     BYPASS THIS SECTION IF THERE ARE NO SHELL DUPLICATES.
      IF(IMJ*KML*IMKJML.EQ.0)GO TO 1610
      CALL SHLOUT(MEND,TQNEW,MEND)
      RETURN
C     SHELL DUPLICATES EXIST - COPY UNIQUE INTEGRALS TO TQOUT.
C     AT THIS STAGE, IT IS NECESSARY TO RE-COMPUTE IIND, JIND, AND KIND.
 1610 ITEMP=0
      DO 1620 K=KSTART,KEND5D
      KIND(K)=ITEMP*LRNG5D
 1620 ITEMP=ITEMP+1
      ITEMP=0
      DO 1650 J=JSTART,JEND5D
      JIND(J)=ITEMP*K5SL5
 1650 ITEMP=ITEMP+1
      ITEMP=0
      DO 1680 I=ISTART,IEND5D
      IIND(I)=ITEMP*J5K5L5
 1680 ITEMP=ITEMP+1
C     CLEAR TQOUT.
      DO 1690 I=1,MEND
 1690 TQOUT(I)=ZERO
C
C     LOOP OVER INDICES AND SET SHELL DUPLICATES TO ZERO.
      JLIM=JEND5D
      KLIM=KEND5D
      DO 1770 I=ISTART,IEND5D
      ITEMP=IIND(I)-LSTP1
      IF(IMJ.EQ.0)JLIM=I
      IF(IMKJML.EQ.0)KLIM=I
      DO 1770 J=JSTART,JLIM
      JTEMP=ITEMP+JIND(J)
      DO 1770 K=KSTART,KLIM
      KTEMP=JTEMP+KIND(K)
      LLIM=LEND5D
      IF(KML.EQ.0)LLIM=K
      IF(IMKJML+IABS(I-K).EQ.0)LLIM=J
      DO 1770 L=LSTART,LLIM
      INDX1=KTEMP+L
 1770 TQOUT(INDX1)=TQNEW(INDX1)
      CALL SHLOUT(MEND,TQOUT,MEND)
      RETURN
      END
