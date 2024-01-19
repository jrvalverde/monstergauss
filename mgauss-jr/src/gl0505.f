C     GL0505       14 JUL 87                                         MRP
C?IBM/GLD/GBR/VAX/UNX
      SUBROUTINE ZCLOSE
C??
C?CDC
C     PROGRAM ZCLOSE
C??
C=SCF OPTIONS
C1INTRODUCTION
C     LINK 0505
C*
C     --------------
C     U OF T VERSION
C     JULY 1987
C     --------------
C*
C     SOLUTION OF THE ROOTHAAN EQUATIONS FOR A CLOSED SHELL SYSTEM
C     FOR ENERGY DECOMPOSITION / BSSE (MODIFIED VERSION OF CLOSED).
C
C     ORIGINAL VERSION: K. MOROKUMA, APRIL 1973.
C     MODIFIED: K. KITAURA, FEBRUARY 1975.
C     MODIFIED: H. UMEYAMA, DECEMBER 1975.
C     MODIFIED: R. CAMMI, U. OF PARMA, DECEMBER 1985 / OCTOBER 1986.
C     RE-WRITTEN: M. PETERSON, U. OF TORONTO CHEMISTRY DEPARTMENT,
C     JUNE 1987.
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
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
C
      COMMON /A/ IOP(99)
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3)
      COMMON /A/ DUMC(4),IDUMC(401)
C
      COMMON/SCF1/ACURCY,DD(NB),D(NB,NB)
      COMMON/SCF2/FF(NB),F(NB,NB)
      COMMON/SCF3/VV(NB),V(NB,NB)
C
      COMMON /SCF6/
     1 IBMAP(NB),IAMAP(NA),IANSAV(NA),ICHSAV(3),IMUSAV(3),NAESAV(3),
     2 NBESAV(3),NESAV(3),NBASAV(3),NOCCAS(2),NVACAS(2),NOCCBS(2),
     3 NVACBS(2),NBACPC(2),ZEN(12,2)
C
      COMMON/GEN/E1,E2,E3,EDUM(3),DCONV,SPIN(2),FX,FY,FZ,FE,EMOL,ESOL,
     1 DPOL(4),VCM
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,IODUM1(5),IFORMT,IODUM2(6),
     1 NFILE(100,2)
C
      EQUIVALENCE (AIBMAP,IBMAP(1))
C
      DATA TEN/10.0D0/
C
 1000 FORMAT('0AT TERMINATION TOTAL ENERGY IS',F17.6,'  HARTREES')
 1010 FORMAT('1CLOSED SHELL ED/BSSE SCF',7X,
     1 'NUCLEAR REPULSION ENERGY IS',F17.9,'  HARTREES'/
     2 '0CONVERGENCE ON DENSITY MATRIX REQUIRED TO ',
     3 'EXIT IS',1PD12.4/'0 CYCLE   ELECTRONIC ENERGY',6X,
     4 'TOTAL ENERGY    CONVERGENCE   EXTRAPOLATION'/)
 1020 FORMAT('0THIS CALCULATION WAS PERFORMED WITH AN APPLIED',
     1 ' EXTERNAL ELECTRIC FIELD')
 1030 FORMAT(/'0THIS CALCULATION WAS DONE ACCORDING TO SELF-CONSISTENT',
     1 ' REACTION FIELD THEORY. ENERGY VALUES IN HARTREES ARE'//
     2 10X,'E(Q)',7X,'=',6X,'EMOL',6X,'+',6X,'ESOL'/1X,3F17.6)
 1040 FORMAT ('1INITIAL DENSITY MATRIX (MO BASIS)'//)
 1070 FORMAT('0ASKED FOR CLOSED SHELL BUT MULTIPLICITY > 1')
C*
      IF(MULTIP.NE.1)GO TO 140
      IF(IOP(10).EQ.0)RETURN
C
C     SEE IF THIS LINK IS REQUIRED FOR THIS PASS.
C
      IF (IOP(4).LE.4 .OR. IOP(4).EQ.6) RETURN
C
C     RECOVER MAP BETWEEN BASIS FUNCTIONS AND MOLECULES.
C
      CALL TREAD (44, AIBMAP, NFILE(44,1), 1, NFILE(44,1), 1, 0)
C
C     SET IOP(93) APPROPRIATELY:
C     0 ... NOT USED.
C     1 ... FOR ES(X).
C     2 ... FOR PL(X).
C     3 ... FOR EX.
C     4 ... FOR CT.
C     5 ... REGULAR RHF.
C
C     FOR ED, GET CT INTERACTION.
C
      IF (IOP(4) .EQ. 5) THEN
         IOP(93) = 4
      ELSE
         IOP(93) = 5
      END IF
C
C     SAVE INITIAL COEFFICIENTS IN FILE 27.
C
      CALL TREAD (25, F, NB, NB, NBASIS, NBASIS, 0)
      CALL TWRITE (27, F, NB, NB, NBASIS, NBASIS, 0)
C     CLEAR ENERGY ORIGIN FLAG.
      IOP(94) = 0
C
C     ZFORMV READS IN THE OVERLAP MATRIX FROM DISK ... FILE 8
C     COMPUTES TRANSFORMATION MATRIX AND WRITES IT ON DISK ... FILE 4
C
      CALL ZFORMV (NBASIS, NBACPC, IOP)
C
C     CREATE INITIAL DENSITY MATRIX (MO BASIS) AND WRITE IT ON FILE 21.
C
      DO 20 I=1,NBASIS
      DO 20 J=1,NBASIS
   20 D(J,I)=ZERO
      DO 30 I=1,NAE
   30 D(I,I)=TWO
      CALL TWRITE (21, D, NB, NB, NBASIS, NBASIS, 1)
C     IS PRINTING OF DENSITY MATRIX DESIRED?
      IF (IOP(34) .EQ. 1) THEN
         WRITE(IOUT,1040)
         CALL MATOUT (D, NB, NB, NBASIS, NBASIS)
      END IF
C
      ACURCY=TEN**(-(5+IOP(13)))
      MAXCYC=30+20*(IOP(14)+1)
C     NEED ONLY ONE CYCLE FOR ES(X) AND TWO FOR ESX+EX.
      IF(IOP(93).EQ.1)MAXCYC=1
      IF(IOP(93).EQ.3)MAXCYC=2
C     CALCULATE NUCLEAR REPULSION ENERGY.
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
     1 (C(I,3)-C(J,3))**2)
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
   90 IF(IOP(35).EQ.0)WRITE(IOUT,1010)E1,ACURCY
C     ENTER ITERATION ROUTINE ZCYCCL - RETURN WHEN THE RMS DIFFERENCE
C     BETWEEN 2 SUCCESSIVE DENSITY MATRICES IS LESS THAN ACURCY
      CALL ZCYCCL (MAXCYC, NBACPC)
C     CONVERGED IF IOP(1) IS 0.
      DCONV=ACURCY
      IF(IOP(35).NE.0)RETURN
      IF(IOP(1).EQ.0)WRITE(IOUT,1000)E3
      IF(IOP(3).EQ.0)RETURN
      IF(IOP(3).EQ.2)GO TO 100
      WRITE(IOUT,1020)
      RETURN
  100 WRITE(IOUT,1030)E3,EMOL,ESOL
      RETURN
C*
  140 WRITE(IOUT,1070)
      IOP(1)=-2
      RETURN
      END
      SUBROUTINE ZFORMV (NBASIS, NBACPC, IOP)
C*
C     --------------
C     U OF T VERSION
C     JULY 1987
C     --------------
C*
C     FORM TRANSFORMATION MATRIX (MODIFIED VERSION OF FORMV).
C
C     INPUT:                        OUTPUT:
C     S(AO) MATRIX ... FILE 8       V(MO) MATRIX ... FILE 4
C     H(AO) MATRIX ... FILE 13      H(MO) MATRIX ... FILE 17
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NB=#NB)
C##
      PARAMETER (NB=200)
C###
      PARAMETER (ONE=1.0D0)
C
      COMMON/SCF1/ACURCY,DD(NB),D(NB,NB)
      COMMON/SCF2/FF(NB),F(NB,NB)
      COMMON/SCF3/VV(NB),V(NB,NB)
      COMMON/IO/IN,IOUT,IODUM(215)
C
      DIMENSION IOP(99), NBACPC(2)
C
 1000 FORMAT ('1OVERLAP MATRIX (INITIAL MO BASIS)'//)
 1010 FORMAT ('1H=T+V MATRIX (INITIAL MO BASIS)'//)
C
C     THE OVERLAP MATRIX IS READ FROM DISK FILE 8
C
      CALL TREAD (8, V, NB, NB, NBASIS, NBASIS, 1)
C
C     TRANSFORM TO MO BASIS AND SELECT.
C
      CALL ZSELCT (V, D, F, 0)
      IF (IOP(32) .EQ. 7) THEN
         WRITE (IOUT,1000)
         CALL MATOUT (V, NB, NB, NBASIS, NBASIS)
      ENDIF
C
C     THE OVERLAP MATRIX IS NOW DIAGONALIZED AND THE V MATRIX FORMED
C     EIGENVECTORS IN D
C     EIGENVALUES IN VV
C     DD IS JUST A SCRATCH ARRAY
C
C     SAVE NBASIS; GET NEW VALUE FOR CPC SCF.
C
      NBASIT = NBASIS
      IF (IOP(4) .GE. 7) NBASIS = NBACPC(1)
      CALL TRED12 (NBASIS, NB, V, D, VV, DD, 2)
      DO 10 J=1,NBASIS
      FACT=ONE/DSQRT(VV(J))
      DO 10 I=1,NBASIS
   10 V(J,I)=D(I,J)*FACT
C
C     RESTORE ORIGINAL VALUE OF NBASIS.
C
      NBASIS = NBASIT
C
C     WRITE OUT TRANSFORMATION MATRIX ON DISK ... FILE 4.
C
      CALL TWRITE (4, V, NB, NB, NBASIS, NBASIS, 0)
C
C     THE CORE MATRIX IS READ FROM FILE 13, TRANSFORMED TO MO BASIS,
C     SELECTED, THEN WRITTEN OUT TO DISK FILE 17.
C
      CALL TREAD (13, V, NB, NB, NBASIS, NBASIS, 1)
      IF (IOP(3) .EQ. 1) CALL EFLD (V, F, NBASIS)
      CALL ZSELCT (V, D, F, 0)
      IF (IOP(32) .EQ. 7) THEN
         WRITE (IOUT,1010)
         CALL MATOUT (V, NB, NB, NBASIS, NBASIS)
      END IF
      CALL TWRITE (17, V, NB, NB, NBASIS, NBASIS, 1)
      RETURN
      END
      SUBROUTINE ZCYCCL (MAXCYC, NBACPC)
C*
C     --------------
C     U OF T VERSION
C     JULY 1987
C     --------------
C*
C     ITERATION ROUTINE (MODIFIED FROM CYCCLO).
C
C     MODIFIED FROM CYCCLO BY KITAURA, FEB. 5 (1975)
C     MONSTERGAUSS VERSION, R. CAMMI, U. DI PARMA, DEC. 1985
C     RE-WRITTEN: M. PETERSON, U. OF TORONTO CHEMISTRY DEPARTMENT,
C     JUNE 1987.
C
C     PARAMETERS AS FOLLOWS:
C     ACURCY ... CONVERGENCE REQUIRED ON DENSITY MATRIX
C     RETURNS WITH ACTUAL CONVERGENCE ACHEIVED
C     MAXCYC ... MAXIMUM NUMBER OF ITERATIONS ALLOWED
C
C     MATRIX FILE USAGE (READ/WRITE UNLESS OTHERWISE NOTED):
C      4: TRANSFORMATION V (ORIGINAL MO BASIS).
C      5: EIGENVALUES (WRITE ONLY).
C      8: OVERLAP S (AO BASIS) (READ ONLY).
C     13: CORE HAMILTONIAN H (AO BASIS) (READ ONLY).
C     17: CORE HAMILTONIAN H (ORIGINAL MO BASIS).
C     19: DENSITY D (AO BASIS) (WRITE ONLY).
C     21: DENSITY D (MO BASIS).
C     25: COEFFICIENTS C (AO / MO BASIS).
C     27: COEFFICIENTS C (AO / ORIGINAL MO BASIS).
C     31-33: DENSITY CONVERGENCE TEST / EXTRAPOLATION.
C     34-36: DIPOLE INTEGRALS X,Y,Z (AO) (READ ONLY).
C     44: ED/BSSE CONTROL DATA (READ ONLY).
C*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C##
      PARAMETER (NA= 36, NB=200)
C###
      PARAMETER (NBB=NB*(NB+1)/2)
      PARAMETER (ZERO=0.0D0)
C
      COMMON /A/ IOP(99)
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3)
      COMMON /A/ DUMC(4),IDUMC(401)
C
      COMMON/SCF1/ACURCY,DD(NB),D(NB,NB)
      COMMON/SCF2/FF(NB),F(NB,NB)
      COMMON/SCF3/VV(NB),V(NB,NB)
C
      COMMON/IO/IN,IOUT,IPUNCH,IMAT,ITWOEL,IODUM(212)
      COMMON/GEN/E1,E2,E3,EDUM(3),DCONV,SPIN(2),FIELD(3),FE,EMOL,ESOL,
     1 DPOL(4),VCM
C
      CHARACTER*8 LABEL
C
      DIMENSION FA(NBB), DA(NBB), NBACPC(2)
C
      EQUIVALENCE (F(1,1),FA(1)), (D(1,1),DA(1))
C
      DATA PT5/0.5D0/, TEN/10.0D0/
C
 1000 FORMAT (1X,I4,1X,2F20.9,1PD14.5,6X,A8)
 1010 FORMAT (///)
 1020 FORMAT ('1COEFFICIENTS (ORIGINAL MO BASIS)'//)
 1030 FORMAT ('0*** SCF DID NOT CONVERGE AFTER',I4,' ITERATIONS ***')
 1040 FORMAT ('1DENSITY MATRIX (ORIGINAL MO BASIS)'//)
 1050 FORMAT ('1DENSITY MATRIX (AO BASIS)')
 1060 FORMAT ('1G=2J+K MATRIX (ORIGINAL MO BASIS)'//)
 1070 FORMAT (' E =',F20.9,10X,'DELTA E =',F15.9,10X,'CONVERGENCE =',
     1 1PD13.5,10X,'RUN ABORTED')
 1080 FORMAT (' E =',F20.9,10X,'DELTA E =',F15.9,10X,'CONVERGENCE =',
     1 1PD13.5,10X,'RUN ALLOWED TO CONTINUE')
 1090 FORMAT ('1FOCK MATRIX (MO BASIS)'//)
 1100 FORMAT ('1COEFFICIENT MATRIX (AO BASIS)')
C*
C     PRESET CONVERGENCE ROUTINES.
C     SETCLO IS AN ENTRY POINT IN CONCLO.
C
      CALL SETCLO
      JCYCLE=1
      NTT=(NBASIS*(NBASIS+1))/2
      IEXT=1
      IPRINT=IOP(33)+IOP(34)
      IOP35=IOP(35)
C
C     PRESET 2-ELECTRON INTEGRAL INPUT ROUTINE.
C
      CALL E2SETC
C
C     SAVE NBASIS FOR LATER.
C
      NBASIT = NBASIS
C
C     JUMP INTO LOOP FOR FIRST TIME ... CALCULATE ENERGY.
C
      E2=ZERO
      GO TO 30
C
C     TOP OF ITERATIONS LOOP.
C
   10 CALL CONCLO (KEY, IEXT, LABEL)
C     PRINT ITERATION RESULTS.
   11 IF(IOP35.NE.0)GO TO 15
      IF(IPRINT.NE.0)WRITE(IOUT,1010)
      IF(IEXT.EQ.0)WRITE(IOUT,1000)JCYCLE,E2,E3,ACURCY,LABEL
      IF(IEXT.NE.0)WRITE(IOUT,1000)JCYCLE,E2,E3
   15 JCYCLE=JCYCLE+1
C
C     MATRIX DA (AO) BACK TO SQUARE FORM D (AO).
C
      CALL SQUARE (DA, D, NBASIS)
C
C     HAS CONVERGENCE BEEN MET ... IF SO ... EXIT
C     IF DENSITY MATRIX HAS JUST BEEN EXTRAPOLATED, CONTINUE FOR 2
C     MORE ITERATIONS TO ALLOW POSSIBLE CONVERGENCE
C
      IF(KEY)30,90,20
C     HAVE WE EXCEEDED THE MAXIMUM ALLOWED NUMBER OF CYCLES ... IF
C     SO ... EXIT
   20 IF(MAXCYC-JCYCLE)91,30,30
   30 IF(IEXT.EQ.0)EPREV=E2
C
C     TRASFORM DENSITY MATRIX (MO BASIS) TO AO BASIS.
C     TRASFORMATION MATRIX IS .... ON DISK ... FILE 27.
C     WRITE  OUT DENSITY MATRIX (AO) ... ON DISK ... FILE 19.
C
      CALL TREAD (27, V, NB, NB, NBASIS, NBASIS, 0)
      CALL MATPAC (V, D, F, NBASIS, 1)
      CALL MATPAC (F, V, D, NBASIS, 3)
      CALL TWRITE (19, D, NB, NB, NBASIS, NBASIS, 1)
      IF (IOP(34) .EQ. 1) THEN
         WRITE (IOUT,1050)
         CALL GBSOUT (D, DD, NB, NB, NBASIS, 0)
      END IF
C
C     LINEARIZE D (AO).
C
      CALL LINEAR (D, DA, NBASIS)
C
C     FORM TWO ELECTRON PART OF THE FOCK MATRIX.
C
      CALL FOFCLO (NTT, IOP, NBASIS)
      IF (IOP(1) .NE. 0) RETURN
C
C     MATRIX FA (AO) BACK TO SQUARE FORM F (AO).
C
      CALL SQUARE (FA, F, NBASIS)
C
C     TRANSFORM F(AO) TO MO BASIS AND SELECT.
C
      CALL ZSELCT (F, D, V, 0)
      IF (IOP(32) .GE. 4) THEN
         WRITE (IOUT,1060)
         CALL MATOUT (F, NB, NB, NBASIS, NBASIS)
      END IF
C
C     THE CORE HAMILTONIAN IS READ FROM DISK INTO MATRIX V.
C
      CALL TREAD (17, V, NB, NB, NBASIS, NBASIS, 1)
C
C     ADD V TO F TO FORM FOCK MATRIX (MO).
C
      DO 59 I=1,NBASIS
      DO 59 J=1,NBASIS
   59 F(J,I)=F(J,I)+V(J,I)
      IF (IOP(32) .GE. 4) THEN
         WRITE (IOUT,1090)
         CALL MATOUT (F, NB, NB, NBASIS, NBASIS)
      END IF
C
C     READ MO DENSITY FROM FILE 21.
C
      CALL TREAD (21, D, NB, NB, NBASIS, NBASIS, 1)
C
C     ENERGY CALCULATION ... E=2*TR(D*H)+TR(D*G).
C
      E2 = ZERO
      DO 60 J=1,NBASIS
      DO 60 I=1,NBASIS
   60 E2 = E2 + (F(I,J)+V(I,J))*D(I,J)
      E2 = PT5 * E2
      E3 = E2 + E1
C
C     SAVE THE FIRST CYCLE ENERGY IN /GEN/.
C
      IF (JCYCLE .EQ. 1) EDUM(1) = E3
C
C     QUIT HERE IF POSSIBLE (ES(X) RUNS).
C
      IF (MAXCYC .EQ. 1) THEN
         IOP35 = 0
         IPRINT = 1
         KEY = 0
         GO TO 11
      END IF
C
C     READ IN THE V MATRIX ... DISK FILE 4.
C
      CALL TREAD (4, V, NB, NB, NBASIS, NBASIS, 0)
C
C     THE F MATRIX IS TRANSFORMED BY THE V MATRIX   TWO STEPS.
C
      CALL MATPAC (V, F, D, NBASIS, 1)
      CALL MATPAC (D, V, F, NBASIS, 3)
C
C     THEN DIAGONALIZED.
C
C     ZERO ALL THE EIGENVALUES BEFORE SUB-SPACE DIAGONALIZED.
C
      DO 70 I=1,NBASIS
   70 FF(I) = ZERO
C
C     GET NEW VALUE OF NBASIS FOR CPC SCF.
C
      IF (IOP(4) .GE. 7) NBASIS = NBACPC(1)
C
      CALL TRED12 (NBASIS, NB, F, D, FF, DD, 2)
C
C     RESTORE ORIGINAL VALUE OF NBASIS.
C
      NBASIS = NBASIT
C
C     WRITE OUT EIGENVALUES ... ARRAY FF ... ON DISK ... FILE 5.
C
      CALL TWRITE (5, FF, NB, 1, NBASIS, 1, 0)
C
C     THE EIGENVECTORS ARE TRANSFORMED BY THE V MATRIX,
C     AND THE MOLECULAR ORBITAL COEFFICIENT MATRIX IS FORMED.
C
      CALL MATPAC (V, D, F, NBASIS, 2)
      IF (IOP(33) .NE. 1) GO TO 78
      WRITE (IOUT,1020)
      CALL MATOUT (F, NB, NB, NBASIS, NBASIS)
C
C     THE DENSITY MATRIX (MO BASIS) IS FORMED FROM THE COEFFICIENTS.
C
   78 DO 80 I=1,NBASIS
      DO 80 J=1,I
      SUM=ZERO
      DO 81 K=1,NAE
   81 SUM=SUM+F(I,K)*F(J,K)
      SUM=SUM+SUM
      D(I,J)=SUM
   80 D(J,I)=SUM
C     WRITE OUT DENSITY MATRIX ON THE DISK ... FILE 21.
      CALL TWRITE (21, D, NB, NB, NBASIS, NBASIS, 1)
      IF(IOP(34).NE.1)GO TO 88
      WRITE(IOUT,1040)
      CALL MATOUT (D, NB, NB, NBASIS, NBASIS)
C
C     THEN TRANSFORMED TO AO BASIS AND WRITTEN OUT ... ON FILE 25.
C
   88 CALL TREAD (27, D, NB, NB, NBASIS, NBASIS, 0)
      CALL MATPAC (D, F, V, NBASIS, 1)
      CALL TWRITE (25, V, NB, NB, NBASIS, NBASIS, 0)
      IF (IOP(33) .EQ. 1) THEN
         WRITE (IOUT,1100)
         CALL GBSOUT (V, FF, NB, NB, NBASIS, 1)
      END IF
C
C     READ DENSITY MATRIX (MO BASIS) FROM FILE 21
C     AND  BACK TO CONVERGENCE ROUTINE.
C
      CALL TREAD (21, DA, NBB, 1, NTT, 1, 0)
      GO TO 10
C
C     MAXIMUM ITERATION COUNT EXCEEDED WITHOUT CONVERGENCE.
C
   91 WRITE(IOUT,1030)MAXCYC
      IOP(1)=-1
      IF(IOP(22).EQ.0)RETURN
C     ALLOW RUN TO CONTINUE IF ENERGY CHANGE IS SMALL ENOUGH.
      EPREV=E2-EPREV
      IF(DABS(EPREV).LE.TEN**(-8+IOP(22)))GO TO 92
      WRITE(IOUT,1070)E3,EPREV,ACURCY
      RETURN
C     ALLOW RUN TO CONTINUE.
   92 WRITE(IOUT,1080)E3,EPREV,ACURCY
      IOP(1)=0
   90 RETURN
      END
      SUBROUTINE ZSELCT (A, B, C, ITYPE)
C*
C     --------------
C     U OF T VERSION
C     JULY 1987
C     --------------
C*
C     TRANSFORM MATRIX A IN AO BASIS TO A IN MO BASIS.
C
C     B AND C ARE SCRATCH MATRICES.
C     ITYPE IS 0 FOR CLOSED SHELL SCF, ALPHA SPIN UHF, OR RHF;
C     1 FOR BETA SPIN UHF.
C
C     THIS PROGRAM APPLIES TO SYSTEMS WITH ONLY TWO SECTIONS.
C
C     ORIGINAL: KITAURA, FEBRUARY 1975.
C     MODIFIED: H. UMEYAMA, DECEMBER 1975.
C     MODIFIED: R. CAMMI, U. OF PARMA, DECEMBER 1985 / OCTOBER 1986.
C     RE-WRITTEN: M. PETERSON, U. OF TORONTO CHEMISTRY DEPARTMENT,
C     JUNE 1987.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA, NB=#NB)
C##
      PARAMETER (NA= 36, NB=200)
C###
      PARAMETER (ZERO=0.0D0)
C
      COMMON /A/ IOP(99),
     1 NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,CDUM(NA,3),
     2 DUMC(4),IDUMC(401)
C
      COMMON /SCF6/
     1 IBMAP(NB),IAMAP(NA),IANSAV(NA),ICHSAV(3),IMUSAV(3),NAESAV(3),
     2 NBESAV(3),NESAV(3),NBASAV(3),NOCCAS(2),NVACAS(2),NOCCBS(2),
     3 NVACBS(2),NBACPC(2),ZEN(12,2)
C
      DIMENSION A(NB,NB),B(NB,NB),C(NB,NB)
C
C
C     SET POINTERS:
C     NOCCA - NUMBER OF OCCUPIED MO IN MOLECULE A.
C     NVACA - NUMBER OF VIRTUAL MO IN MOLECULE A.
C     NOCCB - NUMBER OF OCCUPIED MO IN MOLECULE B.
C     NVACB - NUMBER OF VIRTUAL MO IN MOLECULE B.
C
      NOCCA = NOCCAS(ITYPE+1)
      NOCCB = NOCCBS(ITYPE+1)
      NVACA = NVACAS(ITYPE+1)
      NVACB = NVACBS(ITYPE+1)
      NOCC = NOCCA + NOCCB
      NOCVA = NOCC + NVACA
C
C     READ ISOLATED MO COEFFICIENTS TO B FROM FILE 27.
C
      CALL TREAD (27, B, NB, NB, NBASIS, NBASIS, 0)
C
C     TRANSFORM MATRIX A TO MO BASIS.
C
      CALL MATPAC (B, A, C, NBASIS, 2)
      CALL MATPAC (C, B, A, NBASIS, 1)
C
C     SELECTION OF ELEMENTS BASED ON OPTION 93:
C     1....ES(X).
C     2....PL(X).
C     3....EX.
C     4....CT.
C     5....REGULAR RHF.
C
      IF (IOP(4) .GE. 7) GO TO 700
C
      IF (IOP(93) .EQ. 5) GO TO 600
      IF (IOP(93) .EQ. 2) GO TO 30
C
C     EXCLUDE PL-TERM.
C
      DO 10 J=1,NVACA
      DO 10 I=1,NOCCA
   10 A(I,J+NOCC) = ZERO
      DO 20 J=1,NVACB
      DO 20 I=1,NOCCB
   20 A(I+NOCCA,J+NOCVA) = ZERO
C
   30 IF (IOP(93) .EQ. 3) GO TO 60
C
C     EXCLUDE EX-TERM.
C
      DO 40 J=1,NOCCB
      DO 40 I=1,NOCCA
   40 A(I,J+NOCCA) = ZERO
      DO 50 J=1,NVACB
      DO 50 I=1,NVACA
   50 A(I+NOCC,J+NOCVA) = ZERO
C
   60 IF (IOP(93) .EQ. 4) GO TO 90
C
C     EXCLUDE CT-TERM.
C
      DO 70 J=1,NVACB
      DO 70 I=1,NOCCA
   70 A(I,J+NOCVA) = ZERO
      DO 80 J=1,NVACA
      DO 80 I=1,NOCCB
   80 A(I+NOCCA,J+NOCC) = ZERO
      GO TO 90
C
C     SELECTION OF OPTION BASED ON OPTIONS(16-21).
C
C     IF(IOP(16).NE.0)GO TO 140
C     EXCLUDE POLARIZATION (AOCC-AVAC) TERM BY IOP(16)=0
C     DO 141 I=1,NOCCA
C     DO 141 J=1,NVACA
C 141 A(I,J+NOCC)=ZERO
C 140 CONTINUE
C     IF(IOP(17).NE.0)GO TO 150
C     EXCLUDE POLARIZATION (BOCC-BVAC) TERM BY IOP(17)=0
C     DO 151 I=1,NOCCB
C     DO 151 J=1,NVACB
C 151 A(I+NOCCA,J+NOCVA)=ZERO
C 150 CONTINUE
C     IF(IOP(18).NE.0)GO TO 170
C     EXCLUDE EXCHANGE (AOCC-BOCC) TERM BY IOP(18)=0
C     DO 171 I=1,NOCCA
C     DO 171 J=1,NOCCB
C 171 A(I,J+NOCCA)=ZERO
C 170 CONTINUE
C     IF(IOP(19).NE.0)GO TO 160
C     EXCLUDE EXCHANGE (AVAC-BVAC) TERM BY IOP(19)=0
C     DO 161 I=1,NVACA
C     DO 161 J=1,NVACB
C 161 A(I+NOCC,J+NOCVA)=ZERO
C 160 CONTINUE
C     IF(IOP(20).NE.0)GO TO 120
C     EXCLUDE CT (AOCC-BVAC) TERM BY IOP(20)=0
C     DO 121 I=1,NOCCA
C     DO 121 J=1,NVACB
C 121 A(I,J+NOCVA)=ZERO
C 120 CONTINUE
C     IF(IOP(21).NE.0)GO TO 130
C     EXCLUDE CT (BOCC-AVAC) TERM BY IOP(21)=0
C     DO 131 I=1,NOCCB
C     DO 131 J=1,NVACA
C 131 A(I+NOCCA,J+NOCC)=ZERO
C 130 CONTINUE
C     GO TO 90
C     IOP(15) - IOP(21) TERMINATION
C      SELECT DATA IS USED TO SELECT FOCK ELEMENTS
C 500 CONTINUE
C     IF(NBASIS.EQ.1)GO TO 530
C     DO 520 I=2,NBASIS
C     I1=I-1
C     DO 520 J=1,I1
C     IF (IBMAP(I) .NE. IBMAP(J)) A(J,I)=ZERO
C 520 CONTINUE
C 530 CONTINUE
C     GO TO 90
C
C     CPC SELECTION.
C
  700 NFRST = NBACPC(ITYPE+1) + 1
      IF (NFRST .GT. NBASIS) RETURN
      DO 710 J=NFRST,NBASIS
      DO 710 I=1,NBASIS
  710 A(I,J) = ZERO
C
C     SYMMETRIZE
C
   90 DO 100 I=1,NBASIS
      DO 100 J=1,I
  100 A(I,J) = A(J,I)
C
  600 RETURN
      END
      SUBROUTINE EFLD (F, V, NBASIS)
C*
C     --------------
C     U OF T VERSION
C     JUNE 1987
C     --------------
C*
C     THIS SUBROUTINE EVALUATES THE EXTERNAL ELECTRIC FIELDS PERTUR-
C     BATION MATRIX FOR ED/BSSE.
C
C     ARGUMENTS:
C     F - ON INPUT, CONTAINS MATRIX WHERE FIELD PERTURBATION IS TO
C         BE ADDED; ON OUTPUT, CONTAINS NEW MATRIX (DIMENSION NB*NB).
C     V - SCRATCH MATRIX (DIMENSION NB*NB).
C     NBASIS - NUMBER OF BASIS FUNCTIONS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NB=#NB)
C##
      PARAMETER (NB=200)
C###
      COMMON/GEN/E(6),DCONV,SPIN(2),FIELD(3),FE,EMOL,ESOL,DP(4),VCM
C
      DIMENSION F(NB,NB), V(NB,NB)
C
C     CALCULATE THE EXTERNAL ELECTRIC FIELD PERTURBATION TERMS
C     AND ADD IT TO THE 1-ELECTRON MATRIX.
C
C     CALCULATE X Y AND Z COMPONENTS INDEPENDENTLY:
C     X INTEGRALS ON DISK FILE 34,
C     Y ON 35,
C     Z ON 36.
C
      DO 100 K=1,3
      CALL TREAD (33+K, V, NB, NB, NBASIS, NBASIS, 1)
      FX = FIELD(K)
      DO 100 I=1,NBASIS
      DO 100 J=1,NBASIS
  100 F(I,J) = F(I,J) + FX*V(I,J)
C
      RETURN
      END