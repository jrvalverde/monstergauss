C     GL0308       09 APR 87                                         MRP
C?IBM/GLD/GBR/VAX/UNX
      SUBROUTINE PHOENX
C??
C?CDC
C     PROGRAM PHOENX
C??
C1INTRODUCTION
C     LINK 0308
C*
C     ----------------
C     GAUSSIAN 82
C     U OF T VERSION
C     APRIL 1987
C     ----------------
C*
C     TWO-ELECTRON INTEGRAL PACKAGE FOR SPDF FUNCTIONS, USING THE
C     METHOD OF RYS POLYNOMIALS.
C     NORMALLY ONLY INTEGRALS INVOLVING D/F SHELLS ARE EVALUATED HERE
C     BUT THIS ROUTINE IS CAPABLE OF EVALUATING ALL THE INTEGRALS.
C
C     THE MODE OF OPERATION IS DETERMINED BY IOP(66) - IF PROGRAM SHELL
C     WAS CALLED DURING THIS OVERLAY, PHOENX ASSUMES ALL SP INTEGRALS
C     HAVE ALREADY BEEN EVALUATED. HOWEVER, IF IOP(66) IS 0, THEN ALL
C     THE INTEGRALS ARE EVALUATED HERE - THUS IT IS VERY IMPORTANT TO
C     CALL LINKS 0306 AND 0308 FROM THE SAME OVERLAY CARD IN THE ROUTE.
C     UNCON AND PHOENX (LINKS 0307 AND 0308 RESPECTIVELY) CANNOT BE
C     CALLED TOGETHER - USE ONE OR THE OTHER FOR THE D INTEGRALS.
C*
C     THE ORDER OF THE FUNCTIONS PRODUCED IS:
C      1      S
C      2      X
C      3      Y
C      4      Z
C      5      X**2      D(0)   3*Z**2-R**2
C      6      Z**2      D(1,+) X*Z
C      7      Z**2      D(1,-) Y*Z
C      8      X*Y       D(2,+) X**2-Y**2
C      9      X*Z       D(2,-) X*Y
C     10      Y*Z
C     11      X**3                             F(0)   Z*(5*Z**2-3*R**2)
C     12      Y**3                             F(1,+) X*(5*Z**2-R**2)
C     13      Z**3                             F(1,-) Y*(5*Z**2-R**2)
C     14      X*Y**2                           F(2,+) Z*(X**2-Y**2)
C     15      X**2*Y                           F(2,-) Z*(X*Y)
C     16      X**2*Z                           F(3,+) X*(X**2-3*Y**2)
C     17      X*Z**2                           F(3,-) Y*(3*X**2-Y**2)
C     18      Y*Z**2
C     19      Y**2*Z
C     20      X*Y*Z
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
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      INTEGER SHELLA,SHELLN,SHELLT,SHELLC,AOS,UBOUND
C
      LOGICAL FCNTRL, F0INIT, IFAST
C
      COMMON /A/ IOP(99),ID1(NA),ID2(7),C(NA,3),CDUM(4),ID3(401)
      COMMON/B/IXYZ(NS),SHELLA(NS),SHELLN(NS),SHELLT(NS),
     1 SHELLC(NS),AOS(NS),NSHELL,MAXTYP,EXX(NP),C1(NP),C2(NP)
C
      COMMON /C0308A/ CA(20),CB(20),CC(20),CD(20),TP(7),WP(7),
     1 XIP(256),YIP(256),ZIP(256)
C?GBR/IBM/CDC/VAX/UNX
      COMMON/C0308B/TQ(10000),TQPRM(10000)
C??
C?GLD
C     COMMON/C0308B/TQ(10000)
C     EXTENDED BLOCK/C0308C/TQPRM(10000)
C     COMMON /GLDDATA/ MYID(12)
C??
      COMMON/CCPQ/CCPX(48),CCPY(48),CCPZ(48),CCQX(48),CCQY(48),CCQZ(48)
      COMMON/GCLOOP/IGAUSS,IGBEG,IGEND,IGDF,JGAUSS,JGBEG,JGEND,JGDF,
     1              KGAUSS,KGBEG,KGEND,KGDF,LGAUSS,LGBEG,LGEND,LGDF
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1 EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/LIMIT2/ITYPE,JTYPE,KTYPE,LTYPE,
     1              IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,LSTART,
     2              IEND,JEND,KEND,LEND,LENTQ,
     3              R3OV2,ROOT3,ROOT5,ROOT15,R1,R2,R4,Z1,Z2,Z3
      COMMON/QINFO/XA,YA,ZA,XB,YB,ZB,RABSQ,XC,YC,ZC,XD,YD,ZD,RCDSQ,
     1             PQCUT1,PQCUT2,PQCUT3,
     2             EQSAV(100),QXSAV(100),QYSAV(100),QZSAV(100),
     3             QXPSAV(100),EXPARG,PTEST,PEXP,EP,
     4             KLIND,KLCUTQ(100)
      COMMON /SHLCOM/ ISHELL,JSHELL,KSHELL,LSHELL,IRANGE,JRANGE,KRANGE,
     1 LRANGE
      COMMON/IO/INN,IOUT,IODUM(215)
C
      DIMENSION LBOUND(4),UBOUND(4)
C
      DATA LBOUND/1,1,5,11/, UBOUND/1,4,10,20/
      DATA HALF/0.5D0/, ONEPT5/1.5D0/
      DATA TWOPT5/2.5D0/, THREE/3.0D0/, FOUR/4.0D0/, FIVE/5.0D0/
      DATA PI/3.14159265358979D0/
      DATA TENM20/1.0D-20/, SIXTY/60.0D0/
C
 1000 FORMAT('0PHOENX MAY NOT BE USED WHEN UNCON HAS BEEN CALLED')
 1010 FORMAT(' PHOENIX SHELL PASSES:',I9,', NUMBER PROCESSED:',I9,
     1 ', FAST PASSES:' ,I9)
C*
      IF(IOP(67).EQ.0)GO TO 1
      WRITE(IOUT,1000)
      IOP(1)=-2
      RETURN
C     INITIALIZE CONSTANTS.
    1 IF(IOP(10).EQ.0)RETURN
      CALL FILLST
      CALL SETCP
      IPRINT = 0
      IF (IOP(23) .EQ. 1) IPRINT = 1
      CALL RYSSET (IOP(30), IPRINT)
      PQCUT1 = TENM20
      PQCUT2 = TENM20
      PQCUT3 = SIXTY
C
C     SEE IF THE 'FAST' CODE IS TO BE USED WHERE POSSIBLE.
C     GLD VERSION: FAST CODE IS TOO LARGE, SO USE ONLY REGULAR CODE.
C
C?GBR/IBM/VAX/CDC/UNX
      FCNTRL = .TRUE.
      IF (IOP(31) .EQ. 1) FCNTRL = .FALSE.
C??
      ROOT3=DSQRT(THREE)
      R3OV2=HALF*ROOT3
      ROOT5=DSQRT(FIVE)
      ROOT15 = ROOT3 * ROOT5
      R1=HALF*DSQRT(TWOPT5)
      R2=ONEPT5/ROOT5
      R4=HALF*DSQRT(ONEPT5)
      Z1=FOUR/ROOT5
      Z2=ONE/ROOT5
      Z3=THREE/ROOT5
      IOP66=IOP(66)
      IF(IOP66.EQ.0)CALL SHLOUT(-1,TQ,10000)
      IF (IOP(1) .NE. 0) RETURN
      PCONST=(PI+PI)*PI*DSQRT(PI)
      DO 11 I=1,174
   11 A(I)=ZERO
      DO 12 I=1,48
      CCPX(I)=ZERO
      CCPY(I)=ZERO
      CCPZ(I)=ZERO
      CCQX(I)=ZERO
      CCQY(I)=ZERO
   12 CCQZ(I)=ZERO
C
C     FOR GOULD F77+, ARRAY TQPRM IS LOCATED IN EXTENDED MEMORY.
C?GLD
C     NASK = (10000+MYID(12)-1) / MYID(12)
C     CALL X:GDSPCE (NASK, NGET, ISTAT, *9000)
C??
      KNTSHL=0
      KNTEXE=0
      KNTFST=0
C     ******************************************************************
C     COMMENCE LOOPS OVER SHELLS.
C     ******************************************************************
C
C     LOOP OVER REDUCED SET, ISH,... ETC.
      DO 405 ISH=1,NSHELL
      DO 405 JSH=1,ISH
      DO 405 KSH=1,JSH
      DO 404 LSH=1,KSH
C
C     GIVEN THE REDUCED SET ISH,... ETC, DETERMINE THE PRIMED SET
C     ISHP,... ETC THAT REFLECTS THE PRELIMINARY SWITCHES (IF ANY).
C     THIS PIECE OF INFORMATION, ALONG WITH MTYPE (THE MASTER TYPE)
C     AND NSET (THE NUMBER OF PERMUTED SETS) IS RETURNED BY MTGET.
      MTYPE=MTGET(ISH,JSH,KSH,LSH,ISHP,JSHP,KSHP,LSHP,NSET)
C
C     DETERMINE THE NUMBER OF D- OR F-TYPE FUNCTIONS.  UNDER
C     NORMAL CIRCUMSTANCES, THIS LINK COMPUTES ONLY THOSE INTEGRALS
C     THAT HAVE AT LEAST ONE D- OR F-TYPE FUNCTION.  THIS IS
C     READILY DETERMINED FROM THE SHELL TYPES OF THE FOUR SHELLS
C     AND IS INDEPENDENT OF ANY PERMUTING DONE BY THE ISET LOOP.
      NUMDF=SHELLT(ISH)/2+SHELLT(JSH)/2+SHELLT(KSH)/2+SHELLT(LSH)/2
C     IF IOP(66) IS ZERO, DO ALL INTEGRALS HERE.
      IF(IOP66.NE.0.AND.NUMDF.EQ.0)GO TO 404
      KNTEXE=KNTEXE+1
C
C     WITH THE PRIMED SET AT HAND, ALLOCATE THE REQUIRED AMOUNT OF
C     SPACE IN TQ.
      I=SHELLT(ISHP)+1
      NFA=UBOUND(I)-LBOUND(I)+1
      J=SHELLT(JSHP)+1
      NFB=UBOUND(J)-LBOUND(J)+1
      K=SHELLT(KSHP)+1
      NFC=UBOUND(K)-LBOUND(K)+1
      L=SHELLT(LSHP)+1
      NFD=UBOUND(L)-LBOUND(L)+1
      LENTQ=NFA*NFB*NFC*NFD
C     SEE IF ALL SHELLS ARE CONTRACTED FOR THIS SHELL SET.
      NDC=SHELLN(ISH)*SHELLN(JSH)*SHELLN(KSH)*SHELLN(LSH)
C
C     COMMENCE THE LOOP OVER THE EXPANDED (IE. PERMUTED) SETS.
      DO 403 ISET=1,NSET
      KNTSHL=KNTSHL+1
C
C     DETERMINE ISHELL,... ETC FROM ISHP,... ETC.
      GO TO(10,20,30),ISET
C
C     ISET=1 ... NO PERMUTATION, USE ISHP,... ETC DIRECTLY.
   10 ISHELL=ISHP
      JSHELL=JSHP
      KSHELL=KSHP
      LSHELL=LSHP
      GO TO 40
C
C     ISET=2 ... (AD,BC).
   20 ISHELL=ISHP
      JSHELL=LSHP
      KSHELL=JSHP
      LSHELL=KSHP
      GO TO 40
C
C     ISET=3 ... (AC,BD).
   30 ISHELL=ISHP
      JSHELL=KSHP
      KSHELL=JSHP
      LSHELL=LSHP
C
C     OBTAIN INFORMATION FOR ISHELL.
   40 I=IXYZ(ISHELL)
      XA=C(I,1)
      YA=C(I,2)
      ZA=C(I,3)
      IGBEG=SHELLA(ISHELL)
      IGEND=IGBEG+SHELLN(ISHELL)-1
      ITYPE=SHELLT(ISHELL)
      LAMAX=ITYPE+1
      ISTART=LBOUND(LAMAX)
      IEND=UBOUND(LAMAX)
      IRANGE=IEND-ISTART+1
C
C     OBTAIN INFORMATION FOR JSHELL.
      J=IXYZ(JSHELL)
      XB=C(J,1)
      YB=C(J,2)
      ZB=C(J,3)
      JGBEG=SHELLA(JSHELL)
      JGEND=JGBEG+SHELLN(JSHELL)-1
      JTYPE=SHELLT(JSHELL)
      LBMAX=JTYPE+1
      JSTART=LBOUND(LBMAX)
      JEND=UBOUND(LBMAX)
      JRANGE=JEND-JSTART+1
      LPMAX=LAMAX+JTYPE
      IMJ=IABS(ISHELL-JSHELL)
      RABSQ=(XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2
C
C     OBTAIN INFORMATION FOR KSHELL.
      K=IXYZ(KSHELL)
      XC=C(K,1)
      YC=C(K,2)
      ZC=C(K,3)
      KGBEG=SHELLA(KSHELL)
      KGEND=KGBEG+SHELLN(KSHELL)-1
      KTYPE=SHELLT(KSHELL)
      LCMAX=KTYPE+1
      KSTART=LBOUND(LCMAX)
      KEND=UBOUND(LCMAX)
      KRANGE=KEND-KSTART+1
      IMK=IABS(ISHELL-KSHELL)
C
C     OBTAIN INFORMATION FOR LSHELL.
      L=IXYZ(LSHELL)
      XD=C(L,1)
      YD=C(L,2)
      ZD=C(L,3)
      LGBEG=SHELLA(LSHELL)
      LGEND=LGBEG+SHELLN(LSHELL)-1
      LTYPE=SHELLT(LSHELL)
      LDMAX=LTYPE+1
      LSTART=LBOUND(LDMAX)
      LEND=UBOUND(LDMAX)
      LRANGE=LEND-LSTART+1
      LQMAX=LCMAX+LTYPE
      LPQMAX=LPMAX+LQMAX-1
      JML=IABS(JSHELL-LSHELL)
      KML=IABS(KSHELL-LSHELL)
      IMKJML=IMK+JML
      NZERO=((ITYPE+JTYPE+KTYPE+LTYPE)/2)+1
      RCDSQ=(XD-XC)**2+(YD-YC)**2+(ZD-ZC)**2
C
C     CLEAR THE CONTRACTED FUNCTION INTEGRAL ARRAY TQ.
      DO 50 I=1,LENTQ
   50 TQ(I)=ZERO
C?GBR/IBM/VAX/CDC/UNX
      IFAST = FCNTRL .AND. F0INIT(KNTFST)
C??
C     ******************************************************************
C     COMMENCE LOOPS OVER GAUSSIAN EXPANSION.
C     ******************************************************************
C     PRELIMINARY Q-LOOP.
      CALL QINF
C     IF REQUIRED, GET CONTRACTION STUFF FOR ALTERNATE LOGIC.
C?GBR/IBM/VAX/CDC/UNX
      IF(IFAST)CALL F0CFIL
C??
      INTC=0
      DO 300 IGAUSS=IGBEG,IGEND
      AS=EXX(IGAUSS)
      CALL FILLCP(ITYPE,IGAUSS,C1,C2,CA,CMAXI)
      ASXA=AS*XA
      ASYA=AS*YA
      ASZA=AS*ZA
C
      DO 301 JGAUSS=JGBEG,JGEND
      BS=EXX(JGAUSS)
      CALL FILLCP(JTYPE,JGAUSS,C1,C2,CB,CMAXJ)
C
C     IF REQUIRED, LOAD UP CONTRACTION COEFFICIENTS FOR (IGAUSS,JGAUSS).
C?GBR/IBM/VAX/CDC/UNX
      IF(IFAST)CALL F0CLD1(CA,CB)
C??
      EP=AS+BS
      EPI=ONE/EP
      EP2I=ONE/(EP+EP)
      PX=(ASXA+BS*XB)*EPI
      PY=(ASYA+BS*YB)*EPI
      PZ=(ASZA+BS*ZB)*EPI
      EXPARG=AS*BS*RABSQ*EPI
      IF(EXPARG.GE.PQCUT3)GO TO 301
      PEXP=DEXP(-EXPARG)
      PTEMP=PCONST*PEXP
      PTEST=CMAXI*CMAXJ*PEXP
      IF(PTEST.LT.PQCUT1)GO TO 716
      IJCUTP=0
      GO TO 720
  716 IF(PTEST.GE.PQCUT2)GO TO 718
      IJCUTP=2
      GO TO 720
  718 IJCUTP=1
  720 XAP=PX-XA
      XBP=PX-XB
      YAP=PY-YA
      YBP=PY-YB
      ZAP=PZ-ZA
      ZBP=PZ-ZB
      CALL GETCC2(CCPX,XAP,XBP,LAMAX,LBMAX)
      CALL GETCC2(CCPY,YAP,YBP,LAMAX,LBMAX)
      CALL GETCC2(CCPZ,ZAP,ZBP,LAMAX,LBMAX)
C
      KLIND=0
      DO 302 KGAUSS=KGBEG,KGEND
      CALL FILLC (KTYPE,KGAUSS,C1,C2,CC)
C
      DO 303 LGAUSS=LGBEG,LGEND
      KLIND=KLIND+1
      CALL FILLC (LTYPE,LGAUSS,C1,C2,CD)
C
C     IF REQUIRED, LOAD OUT PRE-COMPUTED CONTRACTION COEFFICIENTS
C     AT (KGAUSS,LGAUSS).
C?GBR/IBM/VAX/CDC/UNX
      IF(IFAST)CALL F0CLD2(KLIND)
C??
C     TEST CUTOFFS.
      IPQCUT=IJCUTP+KLCUTQ(KLIND)
      IF(IPQCUT.GE.2)GO TO 303
      EQ=EQSAV(KLIND)
      QX=QXSAV(KLIND)
      QY=QYSAV(KLIND)
      QZ=QZSAV(KLIND)
      QEXP=QXPSAV(KLIND)
      EPEQ=EP*EQ
      EPPEQI=ONE/(EP+EQ)
      RHO=EPEQ*EPPEQI
      TWORHO=RHO+RHO
      ZTEMP=PTEMP*DSQRT(EPPEQI)*QEXP/EPEQ
      XCQ=QX-XC
      XDQ=QX-XD
      YCQ=QY-YC
      YDQ=QY-YD
      ZCQ=QZ-ZC
      ZDQ=QZ-ZD
      CALL GETCC2(CCQX,XCQ,XDQ,LCMAX,LDMAX)
      CALL GETCC2(CCQY,YCQ,YDQ,LCMAX,LDMAX)
      CALL GETCC2(CCQZ,ZCQ,ZDQ,LCMAX,LDMAX)
      PQX=QX-PX
      PQY=QY-PY
      PQZ=QZ-PZ
      RPQSQ=PQX*PQX+PQY*PQY+PQZ*PQZ
      DXYZ=RHO*RPQSQ
      CALL RPOLX(NZERO,DXYZ,TP,WP)
      CALL GETA2
C
C     IF POSSIBLE, GET CONTRIBUTIONS USING "FAST" LOGIC.
C?GBR/IBM/VAX/CDC/UNX
      IF (IFAST) THEN
         CALL F0ROOT(INTC,TQ,PQX,PQY,PQZ,RHO,ZTEMP,NZERO,TP,WP)
         GO TO 303
      END IF
C??
C     CLEAR THE PRIMITIVE ACCUMULATION ARRAY EACH TIME THROUGH
C     THE CONTRACTION LOOPS.  THIS WAY, THE ZEROES LOOP ACCUMULATES
C     INTO A FRESH AREA EACH TIME.
      DO 80 I=1,LENTQ
   80 TQPRM(I)=ZERO
C
C     ******************************************************************
C     COMMENCE LOOP OVER ZEROES OF RYS POLYNOMIAL.
C     ******************************************************************
C
      DO 118 IZERO=1,NZERO
C     OBTAIN ROOT OF RYS POLYNOMIAL.
      RHOT2=TWORHO*TP(IZERO)
      ZCONST=ZTEMP*WP(IZERO)
C     TEST FOR POSSIBLE BYPASS WITHIN ZEROES LOOP.
      IF(ZCONST.LE.PQCUT2)GO TO 118
      CALL GETIP2(XIP,PQX,ONE,CCPX,CCQX)
      CALL GETIP2(YIP,PQY,ONE,CCPY,CCQY)
      CALL GETIP2(ZIP,PQZ,ZCONST,CCPZ,CCQZ)
C
C     ******************************************************************
C     *  COMMENCE LOOP OVER ATOMIC ORBITALS.                           *
C     ******************************************************************
C
      CALL AOSUMF(INTC,TQPRM,XIP,YIP,ZIP)
C
  118 CONTINUE
C
C     USE ROUTINE CNTPRM TO APPLY CONTRACTION COEFFICIENTS.
      CALL CNTPRM(NDC,INTC,TQ,TQPRM)
C
  303 CONTINUE
  302 CONTINUE
  301 CONTINUE
  300 CONTINUE
C
C     ******************************************************************
C     END OF LOOP OVER GAUSSIANS
C     ******************************************************************
C
C     TEST INTC TO SEE IF ANY INTEGRALS ARE IN TQ.  IF THERE ARE ANY,
C     CALL ROUTINE SHLOUT TO HAVE THEM PROCESSED (OUTPUT).
      IF(INTC.EQ.0)GO TO 403
      CALL PURDF2(INTC,IOP,TQ,TQ,TQPRM,INTCP)
C     INTEGRALS FOR OUTPUT ARE LOCATED IN TQ.
      CALL SHLOUT(INTCP,TQ,INTCP)
      IF (IOP(1) .NE. 0) RETURN
C
C     END OF LOOP OVER SETS.
  403 CONTINUE
C
  404 CONTINUE
  405 CONTINUE
C
C     FOR GOULD F77+, DEALLOCATE EXTENDED MEMORY.
C?GLD
C     CALL X:FDSPCE (NASK, NGET, ISTAT, *9010)
C??
      CALL RYSSET (-1, IPRINT)
      IF (IPRINT .EQ. 0) WRITE (IOUT,1010) KNTSHL, KNTEXE, KNTFST
C
C     EMPTY LAST BUFFER.
      CALL SHLOUT(0,TQ,10000)
      RETURN
C
C     ERROR EXITS.
C?GLD
C9000 WRITE (IOUT,9001) NASK, NGET, ISTAT
C9001 FORMAT ('0*** EXTENDED MEMORY ALLOCATION ERROR IN PHOENX: ',
C    1 'NASK =',I4,', NGET =',I4,', ISTAT =',I4,'. ***')
C     IOP(1) = -2
C     RETURN
C*
C9010 WRITE (IOUT,9011) NASK, NGET, ISTAT
C9011 FORMAT ('0*** EXTENDED MEMORY DEALLOCATION ERROR IN PHOENX: ',
C    1 'NASK =',I4,', NGET =',I4,', ISTAT =',I4,'. ***')
C     IOP(1) = -2
C     RETURN
C??
      END
      SUBROUTINE AOSUMF(INTC,TQ,XIP,YIP,ZIP)
C
C     ******************************************************************
C     FORTRAN ROUTINE TO SUM UP THE A. O. INTEGRALS OVER PRIMITIVE
C     GAUSSIANS FROM THE TWO-DIMENSIONAL INTEGRALS PRODUCED BY GETIP2.
C
C     ASSEMBLER VERSION, AOLOOP, IS NOT AVAILABLE.
C     ******************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER UBOUND
C
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/LIMIT2/ITYPE,JTYPE,KTYPE,LTYPE,
     1              IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,LSTART,
     2              IEND,JEND,KEND,LEND,LENTQ,
     3              R3OV2,ROOT3,ROOT5,ROOT15,R1,R2,R4,Z1,Z2,Z3
C
      DIMENSION TQ(10000),XIP(256),YIP(256),ZIP(256),UBOUND(4)
      DIMENSION INDIX(20),INDIY(20),INDIZ(20),INDJX(20),INDJY(20),
     1 INDJZ(20),INDKX(20),INDKY(20),INDKZ(20),INDLX(20),INDLY(20),
     2 INDLZ(20)
C?GLD
C     EXTENDED DUMMY TQ
C??
      DATA INDLX/1,2,1,1,3,1,1,2,2,1,4,1,1,2,3,3,2,1,1,2/
      DATA INDLY/1,1,2,1,1,3,1,2,1,2,1,4,1,3,2,1,1,2,3,2/
      DATA INDLZ/1,1,1,2,1,1,3,1,2,2,1,1,4,1,1,2,3,3,2,2/
      DATA INDKX/0,4,0,0,8,0,0,4,4,0,12,0,0,4,8,8,4,0,0,4/
      DATA INDKY/0,0,4,0,0,8,0,4,0,4,0,12,0,8,4,0,0,4,8,4/
      DATA INDKZ/0,0,0,4,0,0,8,0,4,4,0,0,12,0,0,4,8,8,4,4/
      DATA INDJX/0,16,0,0,32,0,0,16,16,0,48,0,0,16,32,32,16,0,0,16/
      DATA INDJY/0,0,16,0,0,32,0,16,0,16,0,48,0,32,16,0,0,16,32,16/
      DATA INDJZ/0,0,0,16,0,0,32,0,16,16,0,0,48,0,0,16,32,32,16,16/
      DATA INDIX/0,64,0,0,128,0,0,64,64,0,192,0,0,64,128,128,64,0,0,64/
      DATA INDIY/0,0,64,0,0,128,0,64,0,64,0,192,0,128,64,0,0,64,128,64/
      DATA INDIZ/0,0,0,64,0,0,128,0,64,64,0,0,192,0,0,64,128,128,64,64/
      DATA UBOUND/1,4,10,20/
C*
      INTC=0
      DO 185 I=ISTART,IEND
      IF(IMJ.EQ.0)JEND=I
      IF(IMKJML.EQ.0)KEND=I
      IX=INDIX(I)
      IY=INDIY(I)
      IZ=INDIZ(I)
C
      DO 185 J=JSTART,JEND
      JX=INDJX(J)+IX
      JY=INDJY(J)+IY
      JZ=INDJZ(J)+IZ
C
      DO 185 K=KSTART,KEND
      LEND=UBOUND(LDMAX)
      IF(KML.EQ.0)LEND=K
      IF(IMKJML+IABS(I-K).EQ.0)LEND=J
      KX=INDKX(K)+JX
      KY=INDKY(K)+JY
      KZ=INDKZ(K)+JZ
C
      DO 185 L=LSTART,LEND
      INTC=INTC+1
      LX=INDLX(L)+KX
      LY=INDLY(L)+KY
      LZ=INDLZ(L)+KZ
  185 TQ(INTC)=TQ(INTC)+(XIP(LX)*YIP(LY)*ZIP(LZ))
C
      RETURN
      END
      SUBROUTINE CNTPRM(NDC,INTC,TQ,TQPRIM)
C*
C     ----------------
C     GAUSSIAN 82
C     U OF T VERSION
C     FEBRUARY 1987
C     ----------------
C*
C     ******************************************************************
C     ROUTINE TO CONTRACT THE PRIMITIVE TWO-ELECTRON INTEGRALS
C     INTO THE FINAL CONTRACTED INTEGRALS USING THE CONTRACTION
C     COEFFICIENTS STORED IN /C0308A/.  THIS ROUTINE TAKES AS INPUT
C     THE PRIMITIVES (TQPRIM) AND PRODUCES AS OUTPUT THE FINAL
C     INTEGRALS (TQ).
C     ******************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER UBOUND
C
      COMMON /C0308A/ CA(20),CB(20),CC(20),CD(20),TP(7),WP(7),
     1 XIP(256),YIP(256),ZIP(256)
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/LIMIT2/ITYPE,JTYPE,KTYPE,LTYPE,
     1              IMJ,IMK,JML,KML,IMKJML,ISTART,JSTART,KSTART,LSTART,
     2              IEND,JEND,KEND,LEND,LENTQ,
     3              R3OV2,ROOT3,ROOT5,ROOT15,R1,R2,R4,Z1,Z2,Z3
C
      DIMENSION TQ(10000),TQPRIM(10000),UBOUND(4)
C?GLD
C     EXTENDED DUMMY TQPRIM
C??
      DATA UBOUND/1,4,10,20/
C*
C     APPLY CONTRACTION COEFFICIENTS.
C     HOW THIS IS DONE DEPENDS ON WHETHER OR NOT THERE ARE ANY
C     CONTRACTED FUNCTIONS IN THIS PASS.  IF THERE ARE NOT, USE
C     ALTERNATE CODE.
C
      INTC=0
      IF (NDC .EQ. 1) GO TO 320
C
C     SOME/ALL SHELLS CONTRACTED.
C
      DO 310 I=ISTART,IEND
      CC1=CA(I)
      IF(IMJ.EQ.0)JEND=I
      IF(IMKJML.EQ.0)KEND=I
      DO 310 J=JSTART,JEND
      CC2=CC1*CB(J)
      DO 310 K=KSTART,KEND
      CC3=CC2*CC(K)
      LEND=UBOUND(LDMAX)
      IF(KML.EQ.0)LEND=K
      IF(IMKJML+IABS(I-K).EQ.0)LEND=J
      DO 310 L=LSTART,LEND
      INTC=INTC+1
  310 TQ(INTC)=TQ(INTC)+TQPRIM(INTC)*CC3*CD(L)
      RETURN
C
C     ALL SHELLS UNCONTRACTED.
C
  320 DO 330 I=ISTART,IEND
      CC1=CA(I)
      IF(IMJ.EQ.0)JEND=I
      IF(IMKJML.EQ.0)KEND=I
      DO 330 J=JSTART,JEND
      CC2=CC1*CB(J)
      DO 330 K=KSTART,KEND
      CC3=CC2*CC(K)
      LEND=UBOUND(LDMAX)
      IF(KML.EQ.0)LEND=K
      IF(IMKJML+IABS(I-K).EQ.0)LEND=J
      DO 330 L=LSTART,LEND
      INTC=INTC+1
  330 TQ(INTC)=TQPRIM(INTC)*CC3*CD(L)
      RETURN
      END
      SUBROUTINE F0ROOT(INTC,TQ,PQX,PQY,PQZ,RHO,ZTEMP,NZERO,TP,WP)
C***********************************************************************
C     ROUTINE TO COMPUTE REQUIRED TWO-DIMENSIONAL INTEGRALS
C     AND THEN SUM CONTRIBUTIONS INTO OUTPUT TWO-ELECTRON
C     INTEGRAL ARRAY (TQ).
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/FST2D0/GX(45),GY(45),GZ(45),
     $              X2(125),Y2(125),Z2(125),
     $              X3(225),Y3(225),Z3(225),
     $              X4(405),Y4(405),Z4(405)
      COMMON/FSTCTR/CKLDAT(400),CIJ(100),CKL(100),KLCOND(100)
C
      DIMENSION TQ(10000),TP(7),WP(7),XINT(11)
C
      DATA XINT/1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,6.0D0,7.0D0,8.0D0,9.0D0,
     1 10.0D0,11.0D0/
C
C     INITIALIZATION.
C
C     SPECIAL CODE TO HANDLE (SSSS) CASE.
      IF(LPQMAX.GT.1)GO TO 110
      TQ(1)=TQ(1)+CIJ(1)*CKL(1)*WP(1)*ZTEMP
      INTC=1
      RETURN
C
C     DETERMINE 1-CENTER TERMS.
  110 IG=-LPQMAX
      DO 130 IZERO=1,NZERO
      IG=IG+LPQMAX
      RHOT2=(RHO+RHO)*TP(IZERO)
      GX(1+IG)=XINT(1)
      GY(1+IG)=XINT(1)
      GZ(1+IG)=WP(IZERO)*ZTEMP
      GX(2+IG)=RHOT2*PQX*GX(1+IG)
      GY(2+IG)=RHOT2*PQY*GY(1+IG)
      GZ(2+IG)=RHOT2*PQZ*GZ(1+IG)
      IF(LPQMAX.LE.2)GO TO 130
      DO 120 IV=3,LPQMAX
      GX(IV+IG)=RHOT2*(PQX*GX(IV+IG-1)-XINT(IV-2)*GX(IV+IG-2))
      GY(IV+IG)=RHOT2*(PQY*GY(IV+IG-1)-XINT(IV-2)*GY(IV+IG-2))
  120 GZ(IV+IG)=RHOT2*(PQZ*GZ(IV+IG-1)-XINT(IV-2)*GZ(IV+IG-2))
  130 CONTINUE
C
C     CONVERT TO 2-CENTER TERMS.
      CALL F02CTR(NZERO)
C
C     CONVERT TO 3-CENTER TERMS.
      CALL F03CTR(NZERO)
C
C     CONVERT TO 4-CENTER TERMS.
      CALL F04CTR(NZERO)
C
C     SUM UP INTO TQ.
      CALL F0SUM(NZERO,INTC,TQ)
C
      RETURN
      END
      SUBROUTINE F0SUM(NZERO,INTC,TQ)
C***********************************************************************
C     ROUTINE TO SUM-UP TWO-DIMENSIONAL INTEGRAL CONTRIBUTIONS
C     INTO TWO-ELECTRON INTEGRALS.
C
C     ASSEMBLER VERSION, CALLED FROM F0SUMA, IS NOT AVAILABLE.
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA)
C##
      PARAMETER (NA= 36)
C###
C     LOGICAL FAST
C
      COMMON /A/ IOP(99),ID1(NA),ID2(7),C(NA,3),CDUM(4),ID3(401)
      COMMON/IO/INN,IOUT,IODUM(215)
C
      COMMON/I2ECOM/LAMAX,LBMAX,LCMAX,LDMAX,LPMAX,LQMAX,LPQMAX,L2EFLL,
     1              EQ,EP2I,RHOT2,G(13),VALI2P(49),VALI3P(112),A(174)
      COMMON/FST2D0/GX(45),GY(45),GZ(45),
     $              X2(125),Y2(125),Z2(125),
     $              X3(225),Y3(225),Z3(225),
     $              X4(405),Y4(405),Z4(405)
      COMMON/FSTCTR/CKLDAT(400),CIJ(100),CKL(100),KLCOND(100)
      COMMON/FSTIND/IJXDAT(100),IJYDAT(100),IJZDAT(100),
     $              KLXDAT(100),KLYDAT(100),KLZDAT(100),
     $              NKLDAT(100),NIJ,MAXNKL
C     COMMON/FSTAO/IIIDUM(1296,3),NIJKL
C
      DIMENSION TQ(10000)
C
 1000 FORMAT ('0*** IN F0SUM, LINK 0308, NZERO =',I5,' IS INVALID.')
C
      LEN2D=LAMAX*LBMAX*LCMAX*LDMAX
C
C     POSSIBLY EXECUTE MACHINE DEPENDENT CODE.
C     FAST = NIJKL.GT.0
C     IF(.NOT.(FAST.AND.(NZERO.LE.5))) GOTO 1
C        CALL F0SUMA(NZERO,TQ,INTC,LEN2D)
C        RETURN
C   1 CONTINUE
C
C     SELECT CODE BASED ON NUMBER OF ZEROES.
      GO TO(10,30,50,70,90),NZERO
      WRITE (IOUT,1000) NZERO
      IOP(1) = -2
      RETURN
C
C     SUM UP FOR ONE ZERO.
   10 INTC=0
      DO 20 IJ=1,NIJ
      IJX=IJXDAT(IJ)
      IJY=IJYDAT(IJ)
      IJZ=IJZDAT(IJ)
      NKL=NKLDAT(IJ)
      DO 20 KL=1,NKL
      INTC=INTC+1
      KLX=KLXDAT(KL)+IJX
      KLY=KLYDAT(KL)+IJY
      KLZ=KLZDAT(KL)+IJZ
   20 TQ(INTC)=TQ(INTC)+CIJ(IJ)*CKL(KL)*
     $  ( X4(KLX      )*Y4(KLY      )*Z4(KLZ      ) )
      GO TO 110
C
C     SUM UP FOR TWO ZEROES.
   30 INTC=0
      DO 40 IJ=1,NIJ
      IJX=IJXDAT(IJ)
      IJY=IJYDAT(IJ)
      IJZ=IJZDAT(IJ)
      NKL=NKLDAT(IJ)
      DO 40 KL=1,NKL
      INTC=INTC+1
      KLX=KLXDAT(KL)+IJX
      KLY=KLYDAT(KL)+IJY
      KLZ=KLZDAT(KL)+IJZ
   40 TQ(INTC)=TQ(INTC)+CIJ(IJ)*CKL(KL)*
     $  ( X4(KLX      )*Y4(KLY      )*Z4(KLZ      )
     $   +X4(KLX+LEN2D)*Y4(KLY+LEN2D)*Z4(KLZ+LEN2D) )
      GO TO 110
C
C     SUM UP FOR THREE ZEROES.
   50 LEN22=LEN2D+LEN2D
      INTC=0
      DO 60 IJ=1,NIJ
      IJX=IJXDAT(IJ)
      IJY=IJYDAT(IJ)
      IJZ=IJZDAT(IJ)
      NKL=NKLDAT(IJ)
      DO 60 KL=1,NKL
      INTC=INTC+1
      KLX=KLXDAT(KL)+IJX
      KLY=KLYDAT(KL)+IJY
      KLZ=KLZDAT(KL)+IJZ
   60 TQ(INTC)=TQ(INTC)+CIJ(IJ)*CKL(KL)*
     $  ( X4(KLX      )*Y4(KLY      )*Z4(KLZ      )
     $   +X4(KLX+LEN2D)*Y4(KLY+LEN2D)*Z4(KLZ+LEN2D)
     $   +X4(KLX+LEN22)*Y4(KLY+LEN22)*Z4(KLZ+LEN22) )
      GO TO 110
C
C     SUM UP FOR FOUR ZEROES.
   70 LEN22=LEN2D+LEN2D
      LEN23=LEN2D+LEN22
      INTC=0
      DO 80 IJ=1,NIJ
      IJX=IJXDAT(IJ)
      IJY=IJYDAT(IJ)
      IJZ=IJZDAT(IJ)
      NKL=NKLDAT(IJ)
      DO 80 KL=1,NKL
      INTC=INTC+1
      KLX=KLXDAT(KL)+IJX
      KLY=KLYDAT(KL)+IJY
      KLZ=KLZDAT(KL)+IJZ
   80 TQ(INTC)=TQ(INTC)+CIJ(IJ)*CKL(KL)*
     $  ( X4(KLX      )*Y4(KLY      )*Z4(KLZ      )
     $   +X4(KLX+LEN2D)*Y4(KLY+LEN2D)*Z4(KLZ+LEN2D)
     $   +X4(KLX+LEN22)*Y4(KLY+LEN22)*Z4(KLZ+LEN22)
     $   +X4(KLX+LEN23)*Y4(KLY+LEN23)*Z4(KLZ+LEN23) )
      GO TO 110
C
C     SUM UP FOR FIVE ZEROES.
   90 LEN22=LEN2D+LEN2D
      LEN23=LEN2D+LEN22
      LEN24=LEN2D+LEN23
      INTC=0
      DO 100 IJ=1,NIJ
      IJX=IJXDAT(IJ)
      IJY=IJYDAT(IJ)
      IJZ=IJZDAT(IJ)
      NKL=NKLDAT(IJ)
      DO 100 KL=1,NKL
      INTC=INTC+1
      KLX=KLXDAT(KL)+IJX
      KLY=KLYDAT(KL)+IJY
      KLZ=KLZDAT(KL)+IJZ
  100 TQ(INTC)=TQ(INTC)+CIJ(IJ)*CKL(KL)*
     $  ( X4(KLX      )*Y4(KLY      )*Z4(KLZ      )
     $   +X4(KLX+LEN2D)*Y4(KLY+LEN2D)*Z4(KLZ+LEN2D)
     $   +X4(KLX+LEN22)*Y4(KLY+LEN22)*Z4(KLZ+LEN22)
     $   +X4(KLX+LEN23)*Y4(KLY+LEN23)*Z4(KLZ+LEN23)
     $   +X4(KLX+LEN24)*Y4(KLY+LEN24)*Z4(KLZ+LEN24) )
C
  110 RETURN
      END
      INTEGER FUNCTION MTGET(ISH,JSH,KSH,LSH,ISHP,JSHP,KSHP,LSHP,NSET)
C
C     ******************************************************************
C
C     FORTRAN FUNCTION THAT GIVEN THE RAW SHELL INDICES ISH,... ETC,
C     PRODUCES THE CORRECTLY ORDERED SET ISHP,... ETC, THE NUMBER
C     OF DISTINCT SETS, NSET, AND MTYPE, THE MASTER TYPE INDEX.
C
C     THE FOLLOWING TABLE DEFINES MTYPE:
C
C     ===========================================================
C     MTYPE   DESCRIPTION  # OF SETS  COPY TYPE  SECONDARY  SHELL
C                                                  SWITCH    DUPS
C     ===========================================================
C       1     ALL DISTINCT    3          ---         --       --
C       2     C=D ONLY        2          3=2         NO       2
C       3     B=C ONLY        2          3=1         NO       2
C       4     B=C=D ONLY      1          3=2=1       NO       2
C       5     A=B ONLY        2          3=2         YES      2
C       6     A=B, C=D ONLY   2          3=2         NO       3,4
C       7     A=B=C ONLY      1          3=2=1       YES      2
C       8     A=B=C=D         1          3=2=1       NO       5
C     ===========================================================
C
C     THE SHELL DUPLICATE TYPES ARE CONTAINED IN THE FOLLOWING TABLE.
C
C     ====================
C     TYPE     DESCRIPTION
C     ====================
C      1       A=B ONLY
C      2       C=D ONLY
C      3       A=B AND C=D ONLY
C      4       A=C AND B=D ONLY
C      5       A=B=C=D
C     ====================
C
C     THE FOLLOWING TABLE LISTS THE RESULTANT ORDER OF ALL THREE
C     SETS AFTER ANY PRELIMINARY SETS FOR THE LEGAL VALUES OF
C     MTYPE.
C
C     ===============================
C     MTYPE      FINAL ORDER
C     -------------------------------
C              SET1    SET2    SET3
C     ===============================
C       1     (AB,CD) (AD,BC) (AC,BD)
C       2     (AB,CC) (AC,BC) (AC,BC)
C       3     (AB,BD) (AD,BB) (AB,BD)
C       4     (AB,BB) (AB,BB) (AB,BB)
C       5     (DC,AA) (DA,CA) (DA,CA)
C       6     (AA,CC) (AC,AC) (AC,AC)
C       7     (DA,AA) (DA,AA) (DA,AA)
C       8     (AA,AA) (AA,AA) (AA,AA)
C     ===============================
C
C     IN THE ABOVE TABLES, A,... ETC REFER TO CENTERS.
C
C     ******************************************************************
C
C     THE FOLLOWING STATEMENTS MAP MTYPE INTO THE NUMBER OF SETS.
C
      DIMENSION NSETD(8)
C
      DATA NSETD/3,2,2,1,2,2,1,1/
C
C     DETERMINE MTYPE.
      IAB=0
      IF(ISH.EQ.JSH)IAB=4
      IBC=0
      IF(JSH.EQ.KSH)IBC=2
      ICD=0
      IF(KSH.EQ.LSH)ICD=1
      MTYPE=IAB+IBC+ICD+1
C
C     FROM MTYPE, GET NSET.
      NSET=NSETD(MTYPE)
C
C     PROCESS PRELIMINARY SWITCHES.
C     THE CONDITION UNDER WHICH PRELIMINARY SWITCHES IS MADE IS
C     ISH=JSH, KSH.NE.LSH.
C     HERE, MAKE USE OF MTYPE WITH A COMPUTED GOTO.
      ISHP=ISH
      JSHP=JSH
      KSHP=KSH
      LSHP=LSH
      GO TO(1,1,1,1,5,1,5,1),MTYPE
C     FOR MTYPE = 5 OR 7 PERFORM PRELIMINARY SWITCH.
C     NOTE THAT WE SWITCH NOT ONLY (I,J) WITH (J,L), BUT WE ALSO
C     PERMUTE WITHIN THE PAIRS.
    5 ISHP=LSH
      JSHP=KSH
      KSHP=JSH
      LSHP=ISH
    1 MTGET=MTYPE
      RETURN
      END
