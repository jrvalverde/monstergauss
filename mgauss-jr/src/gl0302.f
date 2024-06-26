C     GL0302       15 JUL 87                                         MRP
C?IBM/GLD/GBR/VAX/UNX
      SUBROUTINE STVINT
C??
C?CDC
C     PROGRAM STVINT
C??
C1INTRODUCTION
C     ----------------
C     GAUSSIAN 82
C     U OF T VERSION
C     JULY 1987
C     ----------------
C*
C     CALCULATION OF OVERLAP AND CORE HAMILTONIAN INTEGRALS.
C*
C1OPTIONS
C     ******************************************************************
C     OPTIONS ... IOP() ... SEE PROGRAM GINPUT (LINK 0301).
C     ******************************************************************
C*
C/
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C     PARAMETER (NS=#NS, NP=#NP)
C##
      PARAMETER (NA= 36, NB=200)
      PARAMETER (NS=120, NP=300)
C###
      PARAMETER (NBB=NB*(NB+1)/2)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
      INTEGER SHELLA,SHELLN,SHELLT,SHELLC,AOS
C
      COMMON /A/ IOP(99)
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3)
      COMMON /A/ CDUM(4),ICDUM(401)
C
      COMMON/B/IXYZ(NS),SHELLA(NS),SHELLN(NS),SHELLT(NS),
     1 SHELLC(NS),AOS(NS),NSHELL,MAXTYP,EXX(NP),C1(NP),C2(NP)
      COMMON/BUF/S(NBB),T(NBB)
      COMMON/LIMIT/IMJ,ISTART,JSTART,IEND,JEND,IRANGE,JRANGE,
     1 ITYPE,JTYPE,IAOS,JAOS,IOP8,JOP8,LIND(NB)
      COMMON/C302A/IXYZA(NS,6),
     1 NSHELA,MAXTPA,EXXA(NP,3)
      COMMON/C302B/TP(4),WP(4),TWOCX(9),TWOCY(9),TWOCZ(9),INDIX(20),
     1 INDIY(20),INDIZ(20),SX(32),SY(32),SZ(32),CA(20),CB(20),
     2 XIP(16),YIP(16),ZIP(16),A(45),SS(100),EEK(100),EEP(100),
     3 CCX(192),CCY(192),CCZ(192),S1C(6)
      COMMON /I1ECOM/ LAMAX, LBMAX, LPMAX
      COMMON/IO/IN,IOUT,IODUM(215)
C
      DIMENSION INDJX(20),INDJY(20),INDJZ(20),IJST(4),IJEND(4)
      DIMENSION SH(NB,NB)
C
      EQUIVALENCE (SH(1,1),S(1))
C
      DATA INDJX/1,2,1,1,3,1,1,2,2,1,4,1,1,2,3,3,2,1,1,2/
      DATA INDJY/1,1,2,1,1,3,1,2,1,2,1,4,1,3,2,1,1,2,3,2/
      DATA INDJZ/1,1,1,2,1,1,3,1,2,2,1,1,4,1,1,2,3,3,2,2/
      DATA IJST/1,1,5,11/, IJEND/1,4,10,20/
      DATA PI/3.14159265358979D0/
      DATA CUT1/-600.0D0/, PT5/0.5D0/
C
 2000 FORMAT('1OVERLAP INTEGRALS')
 2010 FORMAT('1CORE HAMILTONIAN (H=T+V)')
 2020 FORMAT (/I1,'MATRIX DUMP (AO BY AO):')
 2030 FORMAT('0',I3,6(1PD21.12)/(4X,6D21.12))
C*
      IF (IOP(10).EQ.0 .AND. IOP(22).NE.1) RETURN
C     FOR ED/BSSE, IS THIS LINK REQUIRED FOR THIS PASS.
      IF (IOP(4).EQ.1 .OR. IOP(4).EQ.5 .OR. IOP(4).EQ.6 .OR.
     1 IOP(4).EQ.8 .OR. IOP(4).EQ.11) RETURN
      IF ((IOP(4).EQ.9 .OR. IOP(4).EQ.12) .AND. IOP(34).EQ.3) RETURN
      CALL PURSET
      CALL FILLST
      IPRINT = 0
      IF (IOP(23) .EQ. 1) IPRINT = 1
      CALL RYSSET (0, IPRINT)
      IF (IOP(20).NE.0 .AND. IOP(4).LE.3) CALL OVRLAP
      TWOPI=PI+PI
      NTT=(NBASIS*(NBASIS+1))/2
      IOP8 = IOP(8)
      JOP8 = IOP(8)
C     FILL INDIX.
      DO 30 I=1,20
      INDIX(I)=4*(INDJX(I)-1)
      INDIY(I)=4*(INDJY(I)-1)
   30 INDIZ(I)=4*(INDJZ(I)-1)
C
C     FILL INDEXING ARRAY FOR FILMAT.
C
      DO 60 I=1,NBASIS
   60 LIND(I)=(I*(I-1))/2
C
C     LOOP OVER SHELLS.
C
C     LOOP OVER ISHELL.
C
      DO 1000 ISHELL=1,NSHELL
      I=IXYZ(ISHELL)
      XA=C(I,1)
      YA=C(I,2)
      ZA=C(I,3)
      IGBEGN=SHELLA(ISHELL)
      IGEND=IGBEGN+SHELLN(ISHELL)-1
      ITYPE=SHELLT(ISHELL)
      LAMAX=ITYPE+1
      I=IJST(LAMAX)
      ISTART=I+SHELLC(ISHELL)
      IEND=IJEND(LAMAX)
      IRANGE=IEND-ISTART+1
      IAOS=AOS(ISHELL)-I
C
C     LOOP OVER JSHELL.
C
      DO 1000 JSHELL=1,ISHELL
      J=IXYZ(JSHELL)
      XB=C(J,1)
      YB=C(J,2)
      ZB=C(J,3)
      JGBEGN=SHELLA(JSHELL)
      JGEND=JGBEGN+SHELLN(JSHELL)-1
      JTYPE=SHELLT(JSHELL)
      LBMAX=JTYPE+1
      J=IJST(LBMAX)
      JSTART=J+SHELLC(JSHELL)
      JEND=IJEND(LBMAX)
      JRANGE=JEND-JSTART+1
      JAOS=AOS(JSHELL)-J
C
      LPMAX=LAMAX+JTYPE
      LIM1DS=(LPMAX+3)/2
      LENTQ=IRANGE*JRANGE
      IMJ=ISHELL-JSHELL
      NZERO=(ITYPE+JTYPE)/2+1
      ABX=XB-XA
      ABY=YB-YA
      ABZ=ZB-ZA
      RABSQ=ABX*ABX+ABY*ABY+ABZ*ABZ
      DO 90 I=1,LENTQ
      SS(I)=ZERO
   90 EEK(I)=ZERO
C
C     LOOP OVER PRIMITIVE GAUSSIANS.
C
      DO 900 IGAUSS=IGBEGN,IGEND
      AS=EXX(IGAUSS)
      TWOASQ=TWO*AS*AS
      ASXA=AS*XA
      ASYA=AS*YA
      ASZA=AS*ZA
      ARABSQ=AS*RABSQ
      CALL FILLC(ITYPE,IGAUSS,C1,C2,CA)
C
      DO 900 JGAUSS=JGBEGN,JGEND
      BS=EXX(JGAUSS)
      CALL FILLC(JTYPE,JGAUSS,C1,C2,CB)
C
      EP=AS+BS
      EPI=ONE/EP
      EPIO2=PT5*EPI
      TWOP=EP+EP
      ARG=-BS*ARABSQ*EPI
      PEXP=ZERO
      IF(ARG.GT.CUT1)PEXP=DEXP(ARG)
      ZTEMP=TWOPI*EPI*PEXP
      PX=(ASXA+BS*XB)*EPI
      PY=(ASYA+BS*YB)*EPI
      PZ=(ASZA+BS*ZB)*EPI
C
      XAP=PX-XA
      XBP=PX-XB
      YAP=PY-YA
      YBP=PY-YB
      ZAP=PZ-ZA
      ZBP=PZ-ZB
C
      CALL GETCC1(CCX,XAP,XBP,2)
      CALL GETCC1(CCY,YAP,YBP,2)
      CALL GETCC1(CCZ,ZAP,ZBP,2)
C
C     ZERO ACCUMULATION AREA.
C
      DO 110 I=1,LENTQ
  110 EEP(I)=ZERO
C
C     LOOP OVER ATOMS.
C
      DO 800 IATOM=1,NATOMS
C
      IA=IAN(IATOM)
      IF(IA.LE.0)GO TO 800
      XC=C(IATOM,1)
      YC=C(IATOM,2)
      ZC=C(IATOM,3)
      ZT=ZTEMP*DFLOAT(IA)
      PCX=XC-PX
      PCY=YC-PY
      PCZ=ZC-PZ
      RPCSQ=PCX*PCX+PCY*PCY+PCZ*PCZ
      ARG=EP*RPCSQ
      CALL RPOLX(NZERO,ARG,TP,WP)
      CALL GETA1(A,EPIO2,0)
C
C     LOOP OVER ZEROES OF RYS POLYNOMIAL.
C
      DO 700 IZERO=1,NZERO
C
      TWOPT2=TWOP*TP(IZERO)
      ZCONST=ZT*WP(IZERO)
C
      CALL GET2C(TWOCX,PCX,ONE,A,TWOPT2,0)
      CALL GET2C(TWOCY,PCY,ONE,A,TWOPT2,0)
      CALL GET2C(TWOCZ,PCZ,ZCONST,A,TWOPT2,0)
C
      CALL GET3C(XIP,TWOCX,CCX)
      CALL GET3C(YIP,TWOCY,CCY)
      CALL GET3C(ZIP,TWOCZ,CCZ)
C
C     LOOP OVER ATOMIC ORBITALS.
C
      INTC=0
      DO 600 I=ISTART,IEND
      IX=INDIX(I)
      IY=INDIY(I)
      IZ=INDIZ(I)
C
      DO 600 J=JSTART,JEND
      JX=INDJX(J)
      JY=INDJY(J)
      JZ=INDJZ(J)
C
      INTC=INTC+1
      EEP(INTC)=EEP(INTC)+XIP(IX+JX)*YIP(IY+JY)*ZIP(IZ+JZ)
C
  600 CONTINUE
C     END OF AO LOOP.
C
  700 CONTINUE
C     END OF LOOP OVER RYS ZEROES.
C
  800 CONTINUE
C     END OF LOOP OVER ATOMS.
C
C     APPLY THE CONTRACTION COEFFICIENTS.
C     THE POTENTIAL INTEGRALS, ALREADY IN EEP, WILL BE COMBINED WITH
C     THE KINETIC INTEGRALS TO FORM THE CORE HAMILTONIAN (H=T+V).
C
C     CALCULATE THE OVERLAP AND KINETIC ENERGY INTEGRALS.
      STERM=DSQRT(EPI*PI)
      CALL GET1CS(S1C,STERM,EPIO2,2)
      CALL GET2CS(SX,S1C,CCX,2)
      CALL GET2CS(SY,S1C,CCY,2)
      DO 200 I=1,LIM1DS
  200 S1C(I)=S1C(I)*PEXP
      CALL GET2CS(SZ,S1C,CCZ,2)
C
C     BEGIN LOOP OVER ATOMIC ORBITALS FOR OVERLAP, KINETIC ENERGY INTS.
C
      INTC=0
      DO 210 I=ISTART,IEND
      IX=INDIX(I)
      IY=INDIY(I)
      IZ=INDIZ(I)
      IXP=INDJX(I)
      IYP=INDJY(I)
      IZP=INDJZ(I)
      XK=DFLOAT(2*(IXP+IYP+IZP)-3)*AS
      XIIM1=DFLOAT((IXP-1)*(IXP-2))
      YIIM1=DFLOAT((IYP-1)*(IYP-2))
      ZIIM1=DFLOAT((IZP-1)*(IZP-2))
      DO 210 J=JSTART,JEND
      JX=INDJX(J)
      JY=INDJY(J)
      JZ=INDJZ(J)
      IJX=IX+JX
      IJY=IY+JY
      IJZ=IZ+JZ
      INTC=INTC+1
C
      COEF=CA(I)*CB(J)
      SYZ=SY(IJY)*SZ(IJZ)
      TEMPS=SX(IJX)*SYZ
      SS(INTC)=SS(INTC)+TEMPS*COEF
      TEMPK=TEMPS*XK
      TEMPK=TEMPK-TWOASQ*(SX(IJX+8)*SYZ+SX(IJX)*(SY(IJY+8)*SZ(IJZ)+
     1 SY(IJY)*SZ(IJZ+8)))
      TK=ZERO
      IF(IXP.GT.2)TK=XIIM1*SX(IJX-8)*SYZ
      IF(IYP.GT.2)TK=TK+YIIM1*SX(IJX)*SY(IJY-8)*SZ(IJZ)
      IF(IZP.GT.2)TK=TK+ZIIM1*SX(IJX)*SY(IJY)*SZ(IJZ-8)
      TEMPK=TEMPK-PT5*TK
  210 EEK(INTC)=EEK(INTC)+(TEMPK-EEP(INTC))*COEF
C
  900 CONTINUE
C     END OF LOOP OVER GAUSSIANS.
C
C     FILMAT TAKES THE INTEGRALS IN SS AND EEK AND STORES THEM IN
C     THEIR PROPER PLACES IN S AND T.
C
      CALL FILMAT(SS,S,NBB)
      CALL FILMAT(EEK,T,NBB)
C
 1000 CONTINUE
C     END OF LOOP OVER SHELLS.
C
      CALL TWRITE(8,S,NBB,1,NTT,1,0)
      CALL TWRITE(13,T,NBB,1,NTT,1,0)
      CALL RYSSET (-1, IPRINT)
C
C     CHECK THE S/H PRINT OPTION.
C
      IF (IOP(16) .EQ. 0) RETURN
      IF (IOP(16) .EQ. 3) GO TO 1100
      CALL TREAD(8,SH,NB,NB,NBASIS,NBASIS,1)
      WRITE(IOUT,2000)
      CALL GBSOUT(SH,SS,NB,NB,NBASIS,0)
      IF (IOP(16) .EQ. 1) RETURN
C
      IF (IOP(16) .NE. 7) GO TO 1100
      I = 0
      IF (NBASIS .GT. 10) I = 1
      WRITE (IOUT,2020) I
      DO 1010 I=1,NBASIS
 1010 WRITE (IOUT,2030) I,(SH(J,I),J=1,NBASIS)
C
 1100 CALL TREAD(13,SH,NB,NB,NBASIS,NBASIS,1)
      WRITE(IOUT,2010)
      CALL GBSOUT(SH,SS,NB,NB,NBASIS,0)
C
      IF (IOP(16) .NE. 7) RETURN
      I = 0
      IF (NBASIS .GT. 10) I = 1
      WRITE (IOUT,2020) I
      DO 1110 I=1,NBASIS
 1110 WRITE (IOUT,2030) I,(SH(J,I),J=1,NBASIS)
      RETURN
      END
      SUBROUTINE OVRLAP
C*
C     --------------
C     U OF T VERSION
C     FEBRUARY 1987
C     --------------
C*
C     OVERLAP OF BASIS TO PROJECT WITH ITSELF, AND THE CURRENT BASIS.
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C     PARAMETER (NS=#NS, NP=#NP)
C##
      PARAMETER (NA= 36, NB=200)
      PARAMETER (NS=120, NP=300)
C###
      PARAMETER (NB2=NB*NB)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      INTEGER SHELAA,SHELNA,SHELTA,SHELCA,AOSA
      INTEGER SHELAB,SHELNB,SHELTB,SHELCB,AOSB
C
      COMMON /A/ IOP(99)
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3)
      COMMON /A/ CDUM(4),ICDUM(401)
C
      COMMON/C302A/IXYZA(NS),SHELAA(NS),SHELNA(NS),SHELTA(NS),
     1 SHELCA(NS),AOSA(NS),NSHELA,MAXTPA,EXXA(NP),C1A(NP),C2A(NP)
      COMMON/B/IXYZB(NS),SHELAB(NS),SHELNB(NS),SHELTB(NS),
     1 SHELCB(NS),AOSB(NS),NSHELB,MAXTPB,EXXB(NP),C1B(NP),C2B(NP)
      COMMON/BUF/S(NB2),SFILL(NB)
      COMMON/LIMIT/IMJ,ISTART,JSTART,IEND,JEND,IRANGE,JRANGE,
     1 ITYPE,JTYPE,IAOS,JAOS,IOP8,JOP8,LIND(NB)
      COMMON/C302B/TP(4),WP(4),TWOCX(9),TWOCY(9),TWOCZ(9),INDIX(20),
     1 INDIY(20),INDIZ(20),SX(32),SY(32),SZ(32),CA(20),CB(20),
     2 XIP(16),YIP(16),ZIP(16),A(45),SS(100),EEK(100),EEP(100),
     3 CCX(192),CCY(192),CCZ(192),S1C(6)
      COMMON /I1ECOM/ LAMAX, LBMAX, LPMAX
      COMMON/IO/IN,IOUT,IODUM(15),NFILE(100,2)
C
      DIMENSION INDJX(20),INDJY(20),INDJZ(20),IJST(4),IJEND(4)
C
      EQUIVALENCE (IXYZA(1),XYZA), (IXYZB(1),XYZB)
C
      DATA INDJX/1,2,1,1,3,1,1,2,2,1,4,1,1,2,3,3,2,1,1,2/
      DATA INDJY/1,1,2,1,1,3,1,2,1,2,1,4,1,3,2,1,1,2,3,2/
      DATA INDJZ/1,1,1,2,1,1,3,1,2,2,1,1,4,1,1,2,3,3,2,2/
      DATA IJST/1,1,5,11/, IJEND/1,4,10,20/
      DATA PI/3.14159265358979D0/
      DATA CUT1/-600.0D0/, PT5/0.5D0/
C
 1000 FORMAT('1OVERLAP BETWEEN THE BASES (CURRENT X PROJECT)')
C
      MBASIS=IOP(93)
      DO 10 I=1,20
      INDIX(I)=4*(INDJX(I)-1)
      INDIY(I)=4*(INDJY(I)-1)
   10 INDIZ(I)=4*(INDJZ(I)-1)
C
C     GET SELF-OVERLAP FOR BASIS TO BE PROJECTED.
C
      IPASS=0
C     SET INDEXING ARRAY FOR FILMAT, FOR UPPER TRIANGULAR PACKING.
      DO 20 I=1,MBASIS
   20 LIND(I)=(I*(I-1))/2
      CALL TREAD(31,XYZA,NFILE(6,1),1,NFILE(6,1),1,0)
      CALL TREAD(31,XYZB,NFILE(6,1),1,NFILE(6,1),1,0)
      NBB=NB*(NB+1)/2
      NTT=(MBASIS*(MBASIS+1))/2
      IOP8 = 1
      JOP8 = 1
      GO TO 60
C
C     GET OVERLAP BETWEEN THE BASES.
C
   30 IPASS=1
C     SET INDEXING ARRAY FOR FILMAT, FOR RECTANGULAR PACKING.
      DO 40 I=1,MBASIS
   40 LIND(I)=NBASIS*(I-1)
      CALL TREAD(6,XYZB,NFILE(6,1),1,NFILE(6,1),1,0)
      NBB=NB*NB
      NTT=NBASIS*MBASIS
      JOP8 = IOP(8)
C
C     BEGIN LOOP OVER SHELLS.
C
C     LOOP OVER ISHELL.
C
   60 DO 900 ISHELL=1,NSHELA
      I=IXYZA(ISHELL)
      AX=C(I,1)
      AY=C(I,2)
      AZ=C(I,3)
      IGBEGN=SHELAA(ISHELL)
      IGEND=IGBEGN+SHELNA(ISHELL)-1
      ITYPE=SHELTA(ISHELL)
      LAMAX=ITYPE+1
      I=IJST(LAMAX)
      ISTART=I+SHELCA(ISHELL)
      IEND=IJEND(LAMAX)
      IRANGE=IEND-ISTART+1
      IAOS=AOSA(ISHELL)-I
C
C     LOOP OVER JSHELL.
C
      JLIMIT=NSHELB*IPASS+ISHELL*(1-IPASS)
      DO 900 JSHELL=1,JLIMIT
      J=IXYZB(JSHELL)
      BX=C(J,1)
      BY=C(J,2)
      BZ=C(J,3)
      JGBEGN=SHELAB(JSHELL)
      JGEND=JGBEGN+SHELNB(JSHELL)-1
      JTYPE=SHELTB(JSHELL)
      LBMAX=JTYPE+1
      J=IJST(LBMAX)
      JSTART=J+SHELCB(JSHELL)
      JEND=IJEND(LBMAX)
      JRANGE=JEND-JSTART+1
      JAOS=AOSB(JSHELL)-J
C
      LPMAX=LAMAX+JTYPE
      LIM1DS=(LPMAX+1)/2
      LENTQ=IRANGE*JRANGE
      IMJ=IABS(ISHELL-JSHELL)+IPASS
      ABX=BX-AX
      ABY=BY-AY
      ABZ=BZ-AZ
      RABSQ=ABX*ABX+ABY*ABY+ABZ*ABZ
      DO 90 I=1,LENTQ
   90 SS(I)=ZERO
C
C     BEGIN LOOP OVER PRIMITIVE GAUSSIANS.
C
      DO 800 IGAUSS=IGBEGN,IGEND
      AS=EXXA(IGAUSS)
      ASXA=AS*AX
      ASYA=AS*AY
      ASZA=AS*AZ
      ARABSQ=AS*RABSQ
      CALL FILLC(ITYPE,IGAUSS,C1A,C2A,CA)
C
      DO 800 JGAUSS=JGBEGN,JGEND
      BS=EXXB(JGAUSS)
      CALL FILLC(JTYPE,JGAUSS,C1B,C2B,CB)
C
      EP=AS+BS
      EPI=ONE/EP
      EPIO2=PT5*EPI
      ARGP=-BS*ARABSQ*EPI
      PEXP=ZERO
      IF(ARGP.GT.CUT1)PEXP=DEXP(ARGP)
      PX=(ASXA+BS*BX)*EPI
      PY=(ASYA+BS*BY)*EPI
      PZ=(ASZA+BS*BZ)*EPI
      XAP=PX-AX
      XBP=PX-BX
      YAP=PY-AY
      YBP=PY-BY
      ZAP=PZ-AZ
      ZBP=PZ-BZ
      CALL GETCC1(CCX,XAP,XBP,0)
      CALL GETCC1(CCY,YAP,YBP,0)
      CALL GETCC1(CCZ,ZAP,ZBP,0)
      STERM=DSQRT(EPI*PI)
      CALL GET1CS(S1C,STERM,EPIO2,0)
      CALL GET2CS(SX,S1C,CCX,0)
      CALL GET2CS(SY,S1C,CCY,0)
      DO 200 I=1,LIM1DS
  200 S1C(I)=S1C(I)*PEXP
      CALL GET2CS(SZ,S1C,CCZ,0)
C
C     BEGIN LOOP OVER ATOMIC ORBITALS.
C
      INTC=0
      DO 700 I=ISTART,IEND
      IX=INDIX(I)
      IY=INDIY(I)
      IZ=INDIZ(I)
      DO 700 J=JSTART,JEND
      JX=INDJX(J)
      JY=INDJY(J)
      JZ=INDJZ(J)
      INTC=INTC+1
  700 SS(INTC)=SS(INTC)+SX(IX+JX)*SY(IY+JY)*SZ(IZ+JZ)*CA(I)*CB(J)
C     END OF LOOP OVER ATOMIC ORBITALS.
C
  800 CONTINUE
C     END OF LOOP OVER PRIMITIVE GAUSSIANS.
C
      CALL FILMAT(SS,S,NBB)
C
  900 CONTINUE
C     END OF LOOP OVER SHELLS.
C
      IF(IPASS.NE.0)GO TO 910
C
C     STORE SELF-OVERLAP IN FILE 31.
C
      CALL TWRITE(31,S,NBB,1,NTT,1,0)
      GO TO 30
C
C     STORE OVERLAP BETWEEN THE BASES IN FILE 4.
C
  910 CALL TWRITE(4,S,NBB,1,NTT,1,0)
      IF(IOP(16).EQ.0.OR.IOP(16).GE.3)RETURN
      CALL TREAD(4,S,NB,NB,NBASIS,MBASIS,0)
      WRITE(IOUT,1000)
      CALL GBSOUT(S,SS,NB,NB,MBASIS,0)
      RETURN
      END
