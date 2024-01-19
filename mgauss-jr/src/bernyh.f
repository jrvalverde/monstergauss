      PROGRAM TEST                                                      0001.000
C                                                                       0002.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               0003.000
C#                                                                      0004.000
C     PARAMETER (NA=#NA)                                                0005.000
C##                                                                     0006.000
      PARAMETER (NA= 30)                                                0007.000
C###                                                                    0008.000
      PARAMETER (NA3=6*NA, NZ=NA3-6)                                    0009.000
C                                                                       0010.000
      COMMON /A/ IOP(99)                                                0011.000
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3) 0012.000
      COMMON /A/ ANTOAU,DCONST,FCONV,FCCONV,ICDUM(401)                  0013.000
C                                                                       0014.000
      COMMON/ZMAT/IZ(NA,4),ZPARAM(NA,3),IPAR(15,3),NIPAR(3),NPAR,       0015.000
     1 NSTEP,DX(3),NUM,NUMB,LABELS(NA,3),IZMASS(NA),ZMASS(NA)           0016.000
C                                                                       0017.000
      COMMON /IO/ IN,IOUT,IODUM(215)                                    0018.000
C                                                                       0019.000
      CHARACTER TITLE*72                                                0020.000
C                                                                       0021.000
C     INITIALIZATION.                                                   0022.000
C                                                                       0023.000
      WRITE (6,10)                                                      0024.000
   10 FORMAT (' @DOCUMENT(SETUP LANDSCAPE2)')                           0025.000
      I = 2 * NA3 * NA3                                                 0026.000
      CALL M:MAPBLK (NMAP)                                              0027.000
      NASK = (I+NMAP-1) / NMAP                                          0028.000
      CALL X:GDSPCE (NASK, NGET, , )                                    0029.000
C                                                                       0030.000
      IOP = 0                                                           0031.000
      IOP(25) = 2                                                       0032.000
      ICHARG = 0                                                        0033.000
      MULTIP = 1                                                        0034.000
      ANTOAU = 1.0D0 / 0.52917706D0                                     0035.000
      FCCONV = 4.3598149D0                                              0036.000
      IN = 5                                                            0037.000
      IOUT = 6                                                          0038.000
      READ (5,200,END=900) TITLE                                        0039.000
      WRITE (6,210) TITLE                                               0040.000
C                                                                       0041.000
C     GET A Z MATRIX FROM THE 'IN' FILE.                                0042.000
C                                                                       0043.000
  100 CALL FREEZ                                                        0044.000
      IF (IOP(1) .NE. 0) STOP                                           0045.000
      CALL FILMAS (IOP, IN, IOUT, IAN, NATOMS, IZMASS, ZMASS)           0046.000
      CALL BUILDZ                                                       0047.000
C                                                                       0048.000
C     CONVERT C TO ANGSTROMS.                                           0049.000
C                                                                       0050.000
      DO 110 J=1,3                                                      0051.000
      DO 110 I=1,NATOMS                                                 0052.000
  110 C(I,J) = C(I,J) / ANTOAU                                          0053.000
C                                                                       0054.000
C     CALL THE RING HESSIAN ROUTINE.                                    0055.000
C                                                                       0056.000
      CALL BERNYH                                                       0057.000
      WRITE (6,*) 'AFTER BERNYH, IOP(25) =', IOP(25)                    0058.000
      IOP(1) = 0                                                        0059.000
      IOP(25) = 2                                                       0060.000
      READ (IN,200,END=900) TITLE                                       0061.000
  200 FORMAT(A)                                                         0062.000
      WRITE (6,210) TITLE                                               0063.000
  210 FORMAT ('1', A/)                                                  0064.000
      GO TO 100                                                         0065.000
C                                                                       0066.000
  900 STOP                                                              0067.000
      END                                                               0068.000
      SUBROUTINE BERNYH                                                 0069.000
C                                                                       0070.000
C     THIS SUBROUTINE ESTIMATES THE HESSIAN OF RING-TYPE STRUCTURES     0071.000
C     FOR GEOMETRY OPTIMIZATIONS.                                       0072.000
C     THE ALGORITHM FOR THIS SUBROUTINE WAS DESCRIBED IN:               0073.000
C     H.B. SCHLEGEL, THEORET. CHIM. ACTA (BERL), 66, 333-340 (1984).    0074.000
C                                                                       0075.000
C     THE FOLLOWING CHANGES WERE MADE:                                  0076.000
C     1) ALL FORCE CONSTANTS WERE CONVERTED FROM HARTREES/BOHR**2       0077.000
C     TO MDYNE/ANGSTROM.                                                0078.000
C     2) FOR BOND STRETCHES, THE CONSTANTS IN BADGER'S RULE WERE        0079.000
C     TAKEN FROM MONSTERGAUSS ROUTINE 'GBSET' AS THEY COVER THE         0080.000
C     FIRST 5 ROWS OF THE PERIODIC TABLE.                               0081.000
C     3) PROVISION WAS MADE TO HANDLE LINEAR BOND ANGLE BENDS,          0082.000
C     USING THE SAME FORCE CONSTANTS AS REGULAR BENDS.                  0083.000
C     4) NO PROVISION WAS MADE TO HANDLE OUT-OF-PLANE BENDS.            0084.000
C     5) TORSIONS ABOUT A BOND J-K ARE COLLECTED INTO A SINGLE          0085.000
C     TORSION RATHER THAN DEFINING INDIVIDUAL TORSIONS I-J-K-L.         0086.000
C     THE FORCE CONSTANT EXPRESSION WAS CHANGED TO                      0087.000
C     F(TORS) = (AA-BB(R-RCOV))    R<RCOV                               0088.000
C     F(TORS) = AA                 R>=RCOV                              0089.000
C     AND THE CONSTANTS AA AND BB WERE CHANGED BASED ON ESTIMATED       0090.000
C     FORCE CONSTANTS FOR ETHANE AND ETHYLENE AT THE STO-3G AND 3-21G   0091.000
C     LEVELS.                                                           0092.000
C                                                                       0093.000
C     AUTHORS: J.F. MARCOCCIA AND MIKE PETERSON, UNIVERSITY OF TORONTO, 0094.000
C     CHEMISTRY DEPT., TORONTO, ONTARIO, CANADA  M5S 1A1.               0095.000
C     VERSION: JULY 1988.                                               0096.000
C                                                                       0097.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               0098.000
C#                                                                      0099.000
C     PARAMETER (NA=#NA)                                                0100.000
C##                                                                     0101.000
      PARAMETER (NA= 30)                                                0102.000
C###                                                                    0103.000
C ******************************************                            0104.000
C *** BEWARE OF SPECIAL DEFINITION BELOW ***                            0105.000
C ******************************************                            0106.000
      PARAMETER (NA3=6*NA)                                              0107.000
C     DEFINE MAXIMUM NUMBER OF ATOMS FOR EACH END OF A TORSION.         0108.000
      PARAMETER (MAXIL=20)                                              0109.000
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)                                 0110.000
C                                                                       0111.000
      COMMON /A/ IOP(99)                                                0112.000
      COMMON /A/ NATOMS,ICHARG,MULTIP,IAN(NA),NAE,NBE,NE,NBASIS,C(NA,3) 0113.000
      COMMON /A/ ANTOAU,DCONST,FCONV,FCCONV,ICDUM(401)                  0114.000
C                                                                       0115.000
      COMMON /IO/ IN,IOUT,IODUM(215)                                    0116.000
C                                                                       0117.000
      COMMON /C711A/ DISMAT(NA,NA),FVAL(NA3)                            0118.000
      EXTENDED BLOCK /C711B/ B(NA3,NA3)                                 0119.000
C                                                                       0120.000
      DIMENSION A(3), IATOM(MAXIL), LATOM(MAXIL)                        0121.000
      DIMENSION COVRAD(54)                                              0122.000
      DIMENSION FCONST(5,5)                                             0123.000
      DIMENSION WB(NA3)                                                 0124.000
C                                                                       0125.000
      LOGICAL LINEAR                                                    0126.000
C                                                                       0127.000
C     COVDUM IS THE COVALENT RADIUS FOR A DUMMY ATOM (DEFAULT IS 0.5 A).0128.000
C                                                                       0129.000
      DATA COVDUM/0.5D0/                                                0130.000
C                                                                       0131.000
C     COVRAD IS THE COVALENT RADII FOR THE FIRST 54 ELEMENTS OF THE     0132.000
C     PERIODIC TABLE.  THE FIRST 36 RADII ARE FROM H.B. SCHLEGEL        0133.000
C     (REF. ABOVE), THE TRANSITION METAL RADII ARE FROM THE             0134.000
C     SARGENT-WELCH TABLE OF THE ELEMENTS, AND THE REMAINING MAIN GROUP 0135.000
C     RADII ARE FROM: M.R. PETERSON, I.G. CSIZMADIA, THEOCHEM,          0136.000
C     123, 399-412 (1985).                                              0137.000
C                                                                       0138.000
C     FIRST ROW OF THE TABLE:                                           0139.000
      DATA COVRAD/0.32D0, 0.60D0,                                       0140.000
C     SECOND ROW OF THE TABLE:                                          0141.000
     1 1.2D0, 1.05D0, 0.81D0, 0.77D0, 2*0.74D0, 2*0.72D0,               0142.000
C     THIRD ROW OF THE TABLE:                                           0143.000
     2 1.5D0, 1.4D0, 1.3D0, 1.17D0, 1.10D0, 1.04D0, 2*0.99D0,           0144.000
C     FOURTH ROW OF THE TABLE:                                          0145.000
     3 1.8D0, 1.6D0, 1.44D0, 1.32D0, 1.22D0, 1.18D0, 2*1.17D0,          0146.000
     4 1.16D0, 1.15D0, 1.17D0, 1.25D0, 1.4D0, 1.3D0, 2*1.2D0, 2*1.1D0,  0147.000
C     FIFTH ROW OF THE TABLE:                                           0148.000
     5 1.91D0, 1.80D0, 1.62D0, 1.45D0, 1.34D0, 1.30D0, 1.27D0, 2*1.25D0,0149.000
     6 1.28D0, 1.34D0, 1.48D0, 1.42D0, 1.43D0, 1.41D0, 1.39D0, 1.33D0,  0150.000
     7 1.31D0/                                                          0151.000
C                                                                       0152.000
      DATA PT32/0.32D0/, ONEP35/1.35D0/, ONEP3/1.3D0/,                  0153.000
     1 TEN/10.0D0/, PT1/0.1D0/                                          0154.000
C                                                                       0155.000
C     FCONST HAS THE B VALUES USED TO DETERMINE THE FORCE CONSTANTS     0156.000
C     FOR THE BOND STRETCHES IN THE REDUNDANT VALENCE COORDINATES,      0157.000
C     IN MDYNE/ANGSTROM.                                                0158.000
C     VALUES IN FCONST WERE OBTAINED FROM MONSTERGAUSS ROUTINE 'GBSET'. 0159.000
C                                                                       0160.000
      DATA AABL/4.124D0/                                                0161.000
      DATA FCONST/-0.128D0, 0.177D0, 0.368D0, 0.414D0, 0.530D0,         0162.000
     1             0.177D0, 0.562D0, 0.763D0, 0.853D0, 0.996D0,         0163.000
     2             0.368D0, 0.763D0, 1.008D0, 1.078D0, 1.226D0,         0164.000
     3             0.414D0, 0.853D0, 1.078D0, 1.124D0, 1.268D0,         0165.000
     4             0.530D0, 0.996D0, 1.226D0, 1.268D0, 1.373D0/         0166.000
C                                                                       0167.000
C     FORCE CONSTANTS FOR ANGLE BENDS (MDYNE-ANGSTROMS/RADIAN**2).      0168.000
C                                                                       0169.000
      DATA AABE/1.0D0/, AABEH/0.64D0/                                   0170.000
C                                                                       0171.000
C     FORCE CONSTANTS FOR TORSIONS (MDYNE-ANGSTROMS/RADIAN**2).         0172.000
C                                                                       0173.000
      DATA AATO1/0.1D0/, AATO2/3.0D0/                                   0174.000
C                                                                       0175.000
 1000 FORMAT ('0IN BERNYH: ATOM',I3,                                    0176.000
     1 ' IS NOT BONDED TO ANY OTHER ATOMS - ',                          0177.000
     2 'REVERT TO DIAGONAL HESSIAN GUESS (IOP(25)=1).'/                 0178.000
     3 ' COVALENT RADIUS MULTIPLYING FACTOR IS',F6.3/                   0179.000
     4 ' WHICH CAN BE MODIFIED USING IOP(27), GUIDED BY THE DISTANCE',  0180.000
     5 ' MATRIX.')                                                      0181.000
 1010 FORMAT ('0TOO MANY INTERNAL COORDINATES FOUND IN BERNYH - ',      0182.000
     1 'REVERT TO DIAGONAL HESSIAN GUESS (IOP(25)=1).'/                 0183.000
     2 ' NUMBER OF BONDS FOUND:',I5/' NUMBER OF BOND ANGLES FOUND:',I5/ 0184.000
     3 ' NUMBER OF DIHEDRAL ANGLES FOUND:',I5/                          0185.000
     4 ' TOTAL NUMBER OF INTERNAL COORDINATES FOUND:',I5)               0186.000
 1020 FORMAT ('0IN BERNYH: A LINEAR BOND WAS FOUND FOR TORSION AROUND', 0187.000
     1 ' THE',I2,'-',I2,' BOND'/                                        0188.000
     2 ' REVERT TO DIAGONAL HESSIAN GUESS (IOP(25)=1).'/                0189.000
     3 '0*** PLEASE REPORT THIS ERROR TO YOUR MONSTERGAUSS SUPPLIER ***'0190.000
     4 )                                                                0191.000
 1030 FORMAT ('0IN BERNYH: NO SUITABLE REFERENCE AXIS COULD BE FOUND',  0192.000
     1 ' THE',I2,2('-',I2),' LINEAR BOND ANGLE'/                        0193.000
     2 ' REVERT TO DIAGONAL HESSIAN GUESS (IOP(25)=1).'/                0194.000
     3 '0*** PLEASE REPORT THIS ERROR TO YOUR MONSTERGAUSS SUPPLIER ***'0195.000
     4 )                                                                0196.000
 1040 FORMAT ('0IN BERNYH: TOO MANY ATTACHED ATOMS FOR TORSION ABOUT',  0197.000
     1 ' THE',I2,'-',I2,' BOND (LIMIT IS',I4,') -'/                     0198.000
     2 ' REVERT TO DIAGONAL HESSIAN GUESS (IOP(25)=1).'/                0199.000
     3 '0*** PLEASE REPORT THIS ERROR TO YOUR MONSTERGAUSS SUPPLIER ***'0200.000
     4 /' *** AND CHANGE PARAMETER MAXIL IN ROUTINE BERNYH ***')        0201.000
C                                                                       0202.000
C     DETERMINE THE COMPLETE SET OF REDUNDANT INTERNAL VALENCE          0203.000
C     COORDINATES FROM THE CARTESIAN COORDINATES.                       0204.000
C                                                                       0205.000
C     THESE ARE ONLY NORMALLY DEFINED FOR BONDED ATOMS.                 0206.000
C     SO WE SHALL FIRST DETERMINE WHICH ATOMS ARE BONDED AND            0207.000
C     SIMULTANEOUSLY GENERATE A DISTANCE MATRIX.                        0208.000
C                                                                       0209.000
      IF (IOP(27) .EQ. 0) THEN                                          0210.000
         BLFACT = ONEP35                                                0211.000
      ELSE                                                              0212.000
         BLFACT = ONE + DFLOAT(IOP(27))*PT1                             0213.000
      END IF                                                            0214.000
      WRITE (6,*) 'BLFACT =', BLFACT                                    0215.000
C                                                                       0216.000
C     DETERMINE IF FCFACT=1.3 REQUIRED FOR MINIMAL BASIS SET.           0217.000
C                                                                       0218.000
      FCFACT = ONE                                                      0219.000
      IF (IOP(6) .EQ. 0) FCFACT = ONEP3 * FCFACT                        0220.000
      WRITE (6,*) 'FCFACT =', FCFACT                                    0221.000
C                                                                       0222.000
C     INITIALIZE THE ENTIRE B-MATRIX TO ZERO: B(NA3,3*NATOMS).          0223.000
C     SINCE MANY MORE INTERNAL COORDINATES MAY BE FOUND THAN REALLY     0224.000
C     EXIST, ALLOW ALL NA3 ROWS OF THE B MATRIX TO BE USED.             0225.000
C                                                                       0226.000
      DO 100 J=1,3*NATOMS                                               0227.000
      DO 100 I=1,NA3                                                    0228.000
  100 B(I,J) = ZERO                                                     0229.000
C                                                                       0230.000
C     INITIALIZE THE FVAL MATRIX (FORCE CONSTANTS FOR THE REDUNDANT     0231.000
C     VALENCE COORDINATES) TO ZERO.                                     0232.000
C                                                                       0233.000
      DO 110 I=1,NA3                                                    0234.000
  110 FVAL(I) = ZERO                                                    0235.000
C                                                                       0236.000
C     KOUNT IS THE NUMBER OF INTERNAL COORDINATES FOUND.                0237.000
C     KOUNTB IS THE NUMBER OF BONDS FOUND.                              0238.000
C     KOUNTA IS THE NUMBER OF BOND ANGLES FOUND.                        0239.000
C     KOUNTD IS THE NUMBER OF DIHEDRAL ANGLES FOUND.                    0240.000
C                                                                       0241.000
      KOUNT = 0                                                         0242.000
      KOUNTB = 0                                                        0243.000
      KOUNTA = 0                                                        0244.000
      KOUNTD = 0                                                        0245.000
C                                                                       0246.000
C     GET THE COVALENT RADII FOR ATOMS I AND J, THEN DETERMINE IF THEY  0247.000
C     MAKE UP AN ACCEPTABLE BOND FOR THE DISTANCE MATRIX.               0248.000
C                                                                       0249.000
      DO 200 I=1,NATOMS                                                 0250.000
      NI = IAN(I)                                                       0251.000
      IF (NI .LE. 0) THEN                                               0252.000
         COVBLA = COVDUM                                                0253.000
      ELSE                                                              0254.000
         COVBLA = COVRAD(NI)                                            0255.000
      END IF                                                            0256.000
C                                                                       0257.000
      DO 200 J=1,I                                                      0258.000
      NJ = IAN(J)                                                       0259.000
      IF (NJ .LE. 0) THEN                                               0260.000
         COVBLB = COVDUM                                                0261.000
      ELSE                                                              0262.000
         COVBLB = COVRAD(NJ)                                            0263.000
      END IF                                                            0264.000
C                                                                       0265.000
      COVTOT = COVBLA + COVBLB                                          0266.000
C                                                                       0267.000
C     CALCULATE THE DISTANCE BETWEEN ATOM I AND J.                      0268.000
C                                                                       0269.000
      DISTAN = DSQRT( (C(I,1)-C(J,1))**2 +                              0270.000
     1                (C(I,2)-C(J,2))**2 +                              0271.000
     2                (C(I,3)-C(J,3))**2 )                              0272.000
C                                                                       0273.000
      IF (DISTAN.LE.PT32 .OR. DISTAN.GT.BLFACT*COVTOT) THEN             0274.000
         DISMAT(I,J) = ZERO                                             0275.000
         DISMAT(J,I) = ZERO                                             0276.000
      ELSE                                                              0277.000
         DISMAT(I,J) = DISTAN                                           0278.000
         DISMAT(J,I) = DISTAN                                           0279.000
C                                                                       0280.000
C     BOND STRETCH.                                                     0281.000
C                                                                       0282.000
         KOUNT = KOUNT + 1                                              0283.000
         KOUNTB = KOUNTB + 1                                            0284.000
         WRITE (6,99990) KOUNT, I, J                                    0285.000
99990    FORMAT (' INTERNAL COORD',I5,' IS BOND LENGTH ',2I3)           0286.000
         IF (KOUNT .GT. NA3) GO TO 900                                  0287.000
         CALL BOST (KOUNT, I, J, C)                                     0288.000
C                                                                       0289.000
C     CALCULATE THE FVAL TERM FOR BOND STRETCHES: F(STR)=AA/(R-BB)**3,  0290.000
C     WHERE BB DEPENDS ON THE ROW OF PERIODIC TABLE.                    0291.000
C     USE F(STR)=10.0 FOR DUMMY/FLOATING ATOMS.                         0292.000
C                                                                       0293.000
         IF (NI.LE.0 .OR. NJ.LE.0) THEN                                 0294.000
            FVAL(KOUNT) = TEN                                           0295.000
         ELSE                                                           0296.000
            IA = (NI+13) / 8                                            0297.000
            IF (NI .GT. 18) IA = (NI+53) / 18                           0298.000
            IB = (NJ+13) / 8                                            0299.000
            IF (NJ .GT. 18) IB = (NJ+53) / 18                           0300.000
            BB = FCONST(IA,IB)                                          0301.000
            FVAL(KOUNT) = FCFACT * AABL / (DISTAN-BB)**3                0302.000
         END IF                                                         0303.000
C                                                                       0304.000
      END IF                                                            0305.000
C                                                                       0306.000
  200 CONTINUE                                                          0307.000
C                                                                       0308.000
      WRITE (6,*) 'FOUND', KOUNT, ' BONDS'                              0309.000
      WRITE (6,*)                                                       0310.000
      WRITE (6,*) 'DISTANCE MATRIX'                                     0311.000
      WRITE (6,*)                                                       0312.000
      CALL MATOUT (DISMAT, NA, NA, NATOMS, NATOMS)                      0313.000
C                                                                       0314.000
C     CHECK TO MAKE SURE THAT EVERY ATOM IS BONDED AT LEAST ONCE.       0315.000
C                                                                       0316.000
      DO 310 J=1,NATOMS                                                 0317.000
      DO 300 I=1,NATOMS                                                 0318.000
      IF (DISMAT(I,J) .NE. ZERO) GO TO 310                              0319.000
  300 CONTINUE                                                          0320.000
      WRITE (IOUT,1000) J, BLFACT                                       0321.000
      IOP(25) = 1                                                       0322.000
      RETURN                                                            0323.000
  310 CONTINUE                                                          0324.000
C                                                                       0325.000
C     DETERMINE ALL THE POSSIBLE ANGLES FROM DISMAT AND CALCULATE       0326.000
C     THE ROW FOR THE B-MATRIX.                                         0327.000
C                                                                       0328.000
      DO 490 J=1,NATOMS                                                 0329.000
      DO 480 I=1,NATOMS                                                 0330.000
C                                                                       0331.000
C     READ ACROSS THE ENTIRE ROW OF DISMAT SEARCHING FOR ALL            0332.000
C     POSSIBLE I ATOMS.                                                 0333.000
C                                                                       0334.000
      IF (DISMAT(I,J) .EQ. ZERO) GO TO 480                              0335.000
C                                                                       0336.000
C     SEARCH FOR THE NEXT BOND J-K IN THE CURRENT ROW.                  0337.000
C                                                                       0338.000
      DO 470 K=I+1,NATOMS                                               0339.000
      IF (DISMAT(K,J) .EQ. ZERO) GO TO 470                              0340.000
      LINEAR = .FALSE.                                                  0341.000
C                                                                       0342.000
C     DETERMINE ANGLE I-J-K.                                            0343.000
C                                                                       0344.000
      IER = 0                                                           0345.000
      KOUNT = KOUNT + 1                                                 0346.000
      KOUNTA = KOUNTA + 1                                               0347.000
      WRITE (6,99991) KOUNT, I, J, K                                    0348.000
99991 FORMAT (' INTERNAL COORD',I5,' IS BOND ANGLE  ',3I3)              0349.000
      IF (KOUNT .GT. NA3) GO TO 900                                     0350.000
      CALL BEND (KOUNT, I, J, K, C, IER)                                0351.000
C                                                                       0352.000
C     IF BOND ANGLE WAS LINEAR - USE LINEAR BEND WITH X AXIS            0353.000
C     AS A REFERENCE POINT. ONLY KEEP ANGLES NEAR 180 DEGREES.          0354.000
C                                                                       0355.000
      IF (IER .NE. 0) THEN                                              0356.000
         IF (IER .EQ. 1) GO TO 470                                      0357.000
         WRITE (6,*) '*** LINEAR ANGLE', I, J, K                        0358.000
         LINEAR = .TRUE.                                                0359.000
         A(1) = ONE                                                     0360.000
         A(2) = ZERO                                                    0361.000
         A(3) = ZERO                                                    0362.000
C                                                                       0363.000
C     SET L=1, SO THAT THE IN-PLANE COMPONENT IS RETURNED.              0364.000
C                                                                       0365.000
         L = 1                                                          0366.000
C                                                                       0367.000
  450    IER = 0                                                        0368.000
         CALL LIBE (KOUNT, I, J, K, C, A, L, IER)                       0369.000
C                                                                       0370.000
C     IF THIS FAILED WITH THE CURRENT REFERENCE POINT, TRY              0371.000
C     AGAIN WITH THE Y AXIS, THEN THE Z AXIS AS REFERENCE.              0372.000
C     NOTE: AT LEAST ONE OF THESE *MUST* WORK, AS THE BOND              0373.000
C     COULD BE ALIGNED ALONG AT MOST ONE OF THE AXES.                   0374.000
C                                                                       0375.000
         IF (IER .NE. 0) THEN                                           0376.000
            IF (A(1) .NE. ZERO) THEN                                    0377.000
               A(1) = ZERO                                              0378.000
               A(2) = ONE                                               0379.000
               GO TO 450                                                0380.000
            ELSE IF (A(2) .NE. ZERO) THEN                               0381.000
               A(2) = ZERO                                              0382.000
               A(3) = ONE                                               0383.000
               GO TO 450                                                0384.000
            ELSE                                                        0385.000
               WRITE (IOUT,1030) I, J, K                                0386.000
               IOP(25) = 1                                              0387.000
               RETURN                                                   0388.000
            END IF                                                      0389.000
         END IF                                                         0390.000
C                                                                       0391.000
C     SET L=2, SO THAT THE OUT-OF-PLANE (PERPENDICULAR) COMPONENT       0392.000
C     IS RETURNED - NOTE THAT THIS CAN NOT FAIL IF L=1 WORKED.          0393.000
C                                                                       0394.000
         KOUNT = KOUNT + 1                                              0395.000
         KOUNTA = KOUNTA + 1                                            0396.000
         L = 2                                                          0397.000
         CALL LIBE (KOUNT, I, J, K, C, A, L, IER)                       0398.000
C                                                                       0399.000
      END IF                                                            0400.000
C                                                                       0401.000
C     DETERMINE THE FORCE CONSTANT FOR ANGLE BEND: F(BEND)=AA.          0402.000
C                                                                       0403.000
      IF (IAN(I).LE.1 .OR. IAN(K).LE.1) THEN                            0404.000
         AA = AABEH                                                     0405.000
      ELSE                                                              0406.000
         AA = AABE                                                      0407.000
      END IF                                                            0408.000
      FVAL(KOUNT) = FCFACT * AA                                         0409.000
      IF (LINEAR) FVAL(KOUNT-1) = FCFACT * AA                           0410.000
C                                                                       0411.000
  470 CONTINUE                                                          0412.000
  480 CONTINUE                                                          0413.000
  490 CONTINUE                                                          0414.000
C                                                                       0415.000
      WRITE (6,*) 'FOUND', KOUNTA, ' BOND ANGLES'                       0416.000
C                                                                       0417.000
C     DETERMINE ALL POSSIBLE BOND TORSIONS FROM THE DISMAT,             0418.000
C     AND CALCULATE A ROW OF THE B-MATRIX.                              0419.000
C     LOOP OVER ALL PAIRS OF J-K ATOMS (CENTRAL ATOMS).                 0420.000
C                                                                       0421.000
      DO 590 J=1,NATOMS-1                                               0422.000
      NJ = IAN(J)                                                       0423.000
      IF (NJ .LE. 0) THEN                                               0424.000
         COVBLA = COVDUM                                                0425.000
      ELSE                                                              0426.000
         COVBLA = COVRAD(NJ)                                            0427.000
      END IF                                                            0428.000
C                                                                       0429.000
      DO 580 K=J+1,NATOMS                                               0430.000
      IF (DISMAT(K,J) .EQ. ZERO) GO TO 580                              0431.000
      NK = IAN(K)                                                       0432.000
      IF (NK .LE. 0) THEN                                               0433.000
         COVBLB = COVDUM                                                0434.000
      ELSE                                                              0435.000
         COVBLB = COVRAD(NK)                                            0436.000
      END IF                                                            0437.000
      COVTOT = COVBLA + COVBLB                                          0438.000
C                                                                       0439.000
C     TRY TO DEFINE JUST ONE TORSION FOR J-K - THIS CAN BE DONE         0440.000
C     PROVIDED NO OTHER ATOMS ARE 'BONDED' TO BOTH J AND K.             0441.000
C                                                                       0442.000
C     CHECK IF EACH OTHER ATOM IS ATTACHED TO J AND/OR K.               0443.000
C                                                                       0444.000
      NI = 0                                                            0445.000
      NL = 0                                                            0446.000
C                                                                       0447.000
      DO 510 I=1,NATOMS                                                 0448.000
      IF (I.EQ.J .OR. I.EQ.K) GO TO 510                                 0449.000
C                                                                       0450.000
C     LOOK FOR AN 'I'-TYPE ATOM.                                        0451.000
C                                                                       0452.000
      IF (DISMAT(I,J) .EQ. ZERO) GO TO 500                              0453.000
      CALL LINTST (I, J, K, C, IER)                                     0454.000
      IF (IER .NE. 0) GO TO 500                                         0455.000
      NI = NI + 1                                                       0456.000
      IF (NI .GT. MAXIL) THEN                                           0457.000
         WRITE (IOUT,1040) J, K, MAXIL                                  0458.000
         IOP(25) = 1                                                    0459.000
         RETURN                                                         0460.000
      END IF                                                            0461.000
      IATOM(NI) = I                                                     0462.000
C                                                                       0463.000
C     LOOK FOR AN 'L'-TYPE ATOM, BUT IF IT WAS ALSO I-TYPE, QUIT.       0464.000
C                                                                       0465.000
  500 IF (DISMAT(I,K) .EQ. ZERO) GO TO 510                              0466.000
      IF (DISMAT(I,J) .EQ. ZERO) GO TO 550                              0467.000
      CALL LINTST (J, K, I, C, IER)                                     0468.000
      IF (IER .NE. 0) GO TO 510                                         0469.000
      NL = NL + 1                                                       0470.000
      IF (NL .GT. MAXIL) THEN                                           0471.000
         WRITE (IOUT,1040) J, K, MAXIL                                  0472.000
         IOP(25) = 1                                                    0473.000
         RETURN                                                         0474.000
      END IF                                                            0475.000
      LATOM(NL) = I                                                     0476.000
  510 CONTINUE                                                          0477.000
C                                                                       0478.000
      IF (NI.EQ.0 .OR. NL.EQ.0) GO TO 580                               0479.000
C                                                                       0480.000
C     DETERMINE TORSION I-J-K-L.                                        0481.000
C                                                                       0482.000
      KOUNT = KOUNT + 1                                                 0483.000
      KOUNTD = KOUNTD + 1                                               0484.000
      WRITE (6,99992) KOUNT, J, K                                       0485.000
99992 FORMAT (' INTERNAL COORD',I5,' IS BOND TORSION',2I3)              0486.000
      WRITE (6,99993) (IATOM(I),I=1,NI)                                 0487.000
99993 FORMAT (' IATOM(S):',20I3)                                        0488.000
      WRITE (6,99994) (LATOM(I),I=1,NI)                                 0489.000
99994 FORMAT (' LATOM(S):',20I3)                                        0490.000
      IF (KOUNT .GT. NA3) GO TO 900                                     0491.000
      CALL TORS (KOUNT, IATOM, NI, J, K, LATOM, NL, C, IER)             0492.000
C                                                                       0493.000
C     TORSION WAS LINEAR - THIS SHOULD NOT HAPPEN.                      0494.000
C                                                                       0495.000
      IF (IER .NE. 0) THEN                                              0496.000
         WRITE (IOUT,1020) J, K                                         0497.000
         IOP(25) = 1                                                    0498.000
         RETURN                                                         0499.000
      END IF                                                            0500.000
C                                                                       0501.000
C     DETERMINE THE FORCE CONSTANT FOR THE TORSIONAL MODE:              0502.000
C     F(TORS)=(AA-BB(R-RCOV)).                                          0503.000
C                                                                       0504.000
      IF (DISMAT(K,J) .LT. COVTOT) THEN                                 0505.000
         T = AATO2*(DISMAT(K,J)-COVTOT)                                 0506.000
      ELSE                                                              0507.000
         T = ZERO                                                       0508.000
      END IF                                                            0509.000
      FVAL(KOUNT) = FCFACT * (AATO1-T)                                  0510.000
      GO TO 580                                                         0511.000
C                                                                       0512.000
C     DEFINE EACH POSSIBLE TORSION SEPARATELY FOR J-K.                  0513.000
C                                                                       0514.000
C     SEARCH FOR THE FIRST BOND (I,J).                                  0515.000
C                                                                       0516.000
  550 DO 570 I=1,NATOMS                                                 0517.000
      IF (DISMAT(I,J) .EQ. ZERO) GO TO 570                              0518.000
      IF (I.EQ.J .OR. I.EQ.K) GO TO 570                                 0519.000
      CALL LINTST (I, J, K, C, IER)                                     0520.000
      IF (IER .NE. 0) GO TO 570                                         0521.000
C                                                                       0522.000
C     SEARCH FOR THE LAST BOND (K,L).                                   0523.000
C                                                                       0524.000
      DO 560 L=1,NATOMS                                                 0525.000
      IF (DISMAT(L,K) .EQ. ZERO) GO TO 560                              0526.000
      IF (L.EQ.I .OR. L.EQ.J .OR. L.EQ.K) GO TO 560                     0527.000
      CALL LINTST (J, K, L, C, IER)                                     0528.000
      IF (IER .NE. 0) GO TO 560                                         0529.000
C                                                                       0530.000
C     DETERMINE TORSION I-J-K-L.                                        0531.000
C                                                                       0532.000
      KOUNT = KOUNT + 1                                                 0533.000
      KOUNTD = KOUNTD + 1                                               0534.000
      WRITE (6,99995) KOUNT, I, J, K, L                                 0535.000
99995 FORMAT (' INTERNAL COORD',I5,' IS BOND TORSION',4I3)              0536.000
      IF (KOUNT .GT. NA3) GO TO 900                                     0537.000
      IATOM(1) = I                                                      0538.000
      LATOM(1) = L                                                      0539.000
      CALL TORS (KOUNT, IATOM, 1, J, K, LATOM, 1, C, IER)               0540.000
C                                                                       0541.000
C     TORSION WAS LINEAR - THIS SHOULD NOT HAPPEN.                      0542.000
C                                                                       0543.000
      IF (IER .NE. 0) THEN                                              0544.000
         WRITE (IOUT,1020) J, K                                         0545.000
         IOP(25) = 1                                                    0546.000
         RETURN                                                         0547.000
      END IF                                                            0548.000
C                                                                       0549.000
C     DETERMINE THE FORCE CONSTANT FOR THE TORSIONAL MODE:              0550.000
C     F(TORS)=(AA-BB(R-RCOV)).                                          0551.000
C                                                                       0552.000
      IF (DISMAT(K,J) .LT. COVTOT) THEN                                 0553.000
         T = AATO2*(DISMAT(K,J)-COVTOT)                                 0554.000
      ELSE                                                              0555.000
         T = ZERO                                                       0556.000
      END IF                                                            0557.000
      FVAL(KOUNT) = FCFACT * (AATO1-T)                                  0558.000
C                                                                       0559.000
  560 CONTINUE                                                          0560.000
  570 CONTINUE                                                          0561.000
C                                                                       0562.000
  580 CONTINUE                                                          0563.000
  590 CONTINUE                                                          0564.000
C                                                                       0565.000
      WRITE (6,*) 'FOUND', KOUNTD, ' TORSIONAL ANGLES'                  0566.000
      WRITE (6,*)                                                       0567.000
      WRITE (6,*) 'BMAT B'                                              0568.000
      WRITE (6,*)                                                       0569.000
      CALL MATOUTX (B, NA3, NA3, KOUNT, 3*NATOMS)                       0570.000
      WRITE (6,*)                                                       0571.000
      WRITE (6,*) 'FVAL VECTOR'                                         0572.000
      CALL MATOUT (FVAL, NA3, 1, KOUNT, 1)                              0573.000
      WRITE (6,*)                                                       0574.000
C                                                                       0575.000
C     CALCULATE FCART = B' * FVAL * B.                                  0576.000
C                                                                       0577.000
      DO 630 I=1,NA3                                                    0578.000
C                                                                       0579.000
C     COPY COLUMN I OF B(J,I) TO WB, INCORPORATING FVAL(J).             0580.000
C                                                                       0581.000
      DO 600 J=1,KOUNT                                                  0582.000
  600 WB(J) = B(J,I) * FVAL(J)                                          0583.000
C                                                                       0584.000
C     NOW FORM FCART(K,I) AS WB*B(J,K), AND STORE IN B (LOWER TRIANGLE).0585.000
C                                                                       0586.000
      DO 620 K=I,NA3                                                    0587.000
      T = ZERO                                                          0588.000
      DO 610 J=1,KOUNT                                                  0589.000
  610 T = T + WB(J)*B(J,K)                                              0590.000
  620 B(K,I) = T                                                        0591.000
C                                                                       0592.000
  630 CONTINUE                                                          0593.000
C                                                                       0594.000
C     RESTORE UPPER TRIANGLE OF MATRIX.                                 0595.000
C                                                                       0596.000
      DO 640 I=1,NA3                                                    0597.000
      DO 640 J=1,I                                                      0598.000
  640 B(J,I) = B(I,J)                                                   0599.000
C                                                                       0600.000
      WRITE (6,*) 'FCART'                                               0601.000
      WRITE (6,*)                                                       0602.000
      CALL MATOUTX (B, NA3, NA3, 3*NATOMS, 3*NATOMS)                    0603.000
      RETURN                                                            0604.000
C                                                                       0605.000
C     ERROR EXIT - TOO MANY INTERNAL COORDINATES FOUND.                 0606.000
C                                                                       0607.000
  900 WRITE (IOUT,1010) KOUNTB, KOUNTA, KOUNTD, KOUNT                   0608.000
      IOP(25) = 1                                                       0609.000
      RETURN                                                            0610.000
      END                                                               0611.000
      SUBROUTINE BOST (NOB, I, J, C)                                    0612.000
C                                                                       0613.000
C     THIS SUBROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A BOND STRETCH 0614.000
C     AS DEFINED BY WILSON.                                             0615.000
C                                                                       0616.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               0617.000
C#                                                                      0618.000
C     PARAMETER (NA=#NA)                                                0619.000
C##                                                                     0620.000
      PARAMETER (NA= 30)                                                0621.000
C###                                                                    0622.000
      PARAMETER (NA3=6*NA)                                              0623.000
C                                                                       0624.000
      EXTENDED BLOCK /C711B/ B(NA3,NA3)                                 0625.000
C                                                                       0626.000
      DIMENSION C(NA,3), EIJ(3)                                         0627.000
C                                                                       0628.000
      DATA TENM8/1.0D-8/                                                0629.000
C                                                                       0630.000
      CALL VEC (EIJ, C, J, I)                                           0631.000
      II = 3 * (I-1)                                                    0632.000
      JJ = 3 * (J-1)                                                    0633.000
C                                                                       0634.000
      DO 20 M=1,3                                                       0635.000
      T = EIJ(M)                                                        0636.000
      IF (DABS(T) .LT. TENM8) GO TO 20                                  0637.000
      B(NOB,II+M) = -T                                                  0638.000
      B(NOB,JJ+M) = T                                                   0639.000
   20 CONTINUE                                                          0640.000
C                                                                       0641.000
      RETURN                                                            0642.000
      END                                                               0643.000
      SUBROUTINE BEND (NOB, I, J, K, C, IER)                            0644.000
C                                                                       0645.000
C     THIS SUBROUTINE COMPUTES THE B MATRIX ELEMENTS OF A VALENCE       0646.000
C     ANGLE BENDING COORDINATE AS DEFINED BY WILSON.                    0647.000
C                                                                       0648.000
C     I AND K ARE THE NUMBERS OF THE END ATOMS.                         0649.000
C     J IS THE NUMBER OF THE CENTRAL ATOM.                              0650.000
C                                                                       0651.000
C     IER IS A RETURN FLAG:                                             0652.000
C      0 - B MATRIX ROW WAS COMPUTED.                                   0653.000
C      1 - ANGLE I-J-K TOO CLOSE TO LINEAR (0 DEGREES).                 0654.000
C     -1 - ANGLE I-J-K TOO CLOSE TO LINEAR (180 DEGREES).               0655.000
C     IN THE LATTER CASE, THE USER SHOULD PROBABLY CALL THE             0656.000
C     LINEAR BEND ROUTINE 'LIBE'.                                       0657.000
C                                                                       0658.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               0659.000
C#                                                                      0660.000
C     PARAMETER (NA=#NA)                                                0661.000
C##                                                                     0662.000
      PARAMETER (NA= 30)                                                0663.000
C###                                                                    0664.000
      PARAMETER (NA3=6*NA)                                              0665.000
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)                                 0666.000
C                                                                       0667.000
      EXTENDED BLOCK /C711B/ B(NA3,NA3)                                 0668.000
C                                                                       0669.000
      DIMENSION C(NA,3), EJI(3), EJK(3)                                 0670.000
C                                                                       0671.000
      DATA TENM8/1.0D-8/, PT9999/0.99999D0/                             0672.000
C                                                                       0673.000
      DJISQ = ZERO                                                      0674.000
      DJKSQ = ZERO                                                      0675.000
C                                                                       0676.000
      DO 20 M=1,3                                                       0677.000
      TP = C(J,M)                                                       0678.000
      T = C(I,M) - TP                                                   0679.000
      EJI(M) = T                                                        0680.000
      DJISQ = DJISQ + T*T                                               0681.000
      T = C(K,M) - TP                                                   0682.000
      EJK(M) = T                                                        0683.000
   20 DJKSQ = DJKSQ + T*T                                               0684.000
C                                                                       0685.000
      DJI = DSQRT(DJISQ)                                                0686.000
      DJK = DSQRT(DJKSQ)                                                0687.000
      DOTJ = ZERO                                                       0688.000
C                                                                       0689.000
      DO 30 M=1,3                                                       0690.000
      T = EJI(M) / DJI                                                  0691.000
      EJI(M) = T                                                        0692.000
      TP = EJK(M) / DJK                                                 0693.000
      EJK(M) = TP                                                       0694.000
   30 DOTJ = DOTJ + T*TP                                                0695.000
C                                                                       0696.000
      IF (DABS(DOTJ) .GT. PT9999) GO TO 60                              0697.000
      SINJ = DSQRT(ONE-DOTJ*DOTJ)                                       0698.000
      II = 3 * (I-1)                                                    0699.000
      JJ = 3 * (J-1)                                                    0700.000
      KK = 3 * (K-1)                                                    0701.000
C                                                                       0702.000
      DO 40 M=1,3                                                       0703.000
      SMI = (DOTJ*EJI(M)-EJK(M)) / (DJI*SINJ)                           0704.000
      IF (DABS(SMI) .GE. TENM8) B(NOB,II+M) = SMI                       0705.000
      SMK = (DOTJ*EJK(M)-EJI(M)) / (DJK*SINJ)                           0706.000
      IF (DABS(SMK) .GE. TENM8) B(NOB,KK+M) = SMK                       0707.000
      SUM = SMI+SMK                                                     0708.000
   40 IF (DABS(SUM) .GE. TENM8) B(NOB,JJ+M) = -SUM                      0709.000
C                                                                       0710.000
      RETURN                                                            0711.000
C                                                                       0712.000
C     ERROR: LINEAR BEND AS COS(I-J-K) IS CLOSE TO +/-1.0.              0713.000
C                                                                       0714.000
   60 IF (DOTJ .GT. ZERO) THEN                                          0715.000
         IER = 1                                                        0716.000
      ELSE                                                              0717.000
         IER = -1                                                       0718.000
      END IF                                                            0719.000
      RETURN                                                            0720.000
      END                                                               0721.000
      SUBROUTINE TORS (NOB, IATOM, NI, J, K, LATOM, NL, C, IER)         0722.000
C                                                                       0723.000
C     THIS SUBROUTINE COMPUTES THE B MATRIX ELEMENTS FOR TORSION AS     0724.000
C     DEFINED BY R L HILDERBRANDT IN J MOLEC SPEC, 44, 599 (1972).      0725.000
C                                                                       0726.000
C     CODED BY M PETERSON, DEPT OF CHEMISTRY, UNIV OF TORONTO.          0727.000
C                                                                       0728.000
C     J AND K DEFINE THE BOND UNDER TORSION.                            0729.000
C     IATOM, LATOM: ATOM NUMBERS FOR THE I- AND L-TYPE ATOMS.           0730.000
C     NI AND NL ARE THE NUMBER OF I-TYPE ATOMS ATTACHED TO J AND        0731.000
C     THE NUMBER OF L-TYPE ATOMS ATTACHED TO K RESPECTIVELY.            0732.000
C                                                                       0733.000
C     IER IS A RETURN FLAG:                                             0734.000
C      0 - TORSION COMPUTED AND ADDED TO B MATRIX ROW.                  0735.000
C      1 - ATOMS I-J-K ARE COLINEAR.                                    0736.000
C     -2 - ATOMS J-K-L ARE COLINEAR.                                    0737.000
C     IN EITHER OF THE LATTER CASES, THE USER MAY WANT TO USE THE       0738.000
C     LINEAR BEND ROUTINE 'LIBE'.                                       0739.000
C                                                                       0740.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               0741.000
C#                                                                      0742.000
C     PARAMETER (NA=#NA)                                                0743.000
C##                                                                     0744.000
      PARAMETER (NA= 30)                                                0745.000
C###                                                                    0746.000
      PARAMETER (NA3=6*NA)                                              0747.000
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)                                 0748.000
C                                                                       0749.000
      EXTENDED BLOCK /C711B/ B(NA3,NA3)                                 0750.000
C                                                                       0751.000
      DIMENSION C(NA,3), SJ(3), SK(3), EJK(3), EIJ(3), ELK(3), CR(3)    0752.000
C                                                                       0753.000
      DIMENSION IATOM(NI), LATOM(NL)                                    0754.000
C                                                                       0755.000
      DATA TENM8/1.0D-8/, PT9999/0.99999D0/                             0756.000
C                                                                       0757.000
      DJKSQ = ZERO                                                      0758.000
C                                                                       0759.000
      DO 20 M=1,3                                                       0760.000
      SJ(M) = ZERO                                                      0761.000
      SK(M) = ZERO                                                      0762.000
      T = C(K,M) - C(J,M)                                               0763.000
      EJK(M) = T                                                        0764.000
   20 DJKSQ = DJKSQ + T*T                                               0765.000
C                                                                       0766.000
      DJK = ONE / DSQRT(DJKSQ)                                          0767.000
C                                                                       0768.000
      DO 30 M=1,3                                                       0769.000
   30 EJK(M) = EJK(M) * DJK                                             0770.000
C                                                                       0771.000
      JJ = 3 * (J-1)                                                    0772.000
      KK = 3 * (K-1)                                                    0773.000
C                                                                       0774.000
C     LOOP OVER THE I-TYPE ATOMS.                                       0775.000
C                                                                       0776.000
      DO 60 N=1,NI                                                      0777.000
      I = IATOM(N)                                                      0778.000
      DIJSQ = ZERO                                                      0779.000
C                                                                       0780.000
      DO 40 M=1,3                                                       0781.000
      T = C(J,M) - C(I,M)                                               0782.000
      EIJ(M) = T                                                        0783.000
   40 DIJSQ = DIJSQ + T*T                                               0784.000
C                                                                       0785.000
      DIJ = ONE / DSQRT(DIJSQ)                                          0786.000
      COSJ = ZERO                                                       0787.000
C                                                                       0788.000
      DO 50 M=1,3                                                       0789.000
      T = EIJ(M) * DIJ                                                  0790.000
      EIJ(M) = T                                                        0791.000
   50 COSJ = COSJ - T*EJK(M)                                            0792.000
C                                                                       0793.000
      IF (DABS(COSJ) .GT. PT9999) GO TO 120                             0794.000
      SIN2J = (ONE-COSJ*COSJ) * DFLOAT(NI)                              0795.000
      II = 3 * (I-1)                                                    0796.000
      CALL VPROD (CR, EIJ, EJK)                                         0797.000
C                                                                       0798.000
      DO 60 M=1,3                                                       0799.000
      T = CR(M) / SIN2J                                                 0800.000
      SMI = T * DIJ                                                     0801.000
      IF (DABS(SMI) .GE. TENM8) B(NOB,II+M) = -SMI                      0802.000
      SMK = T * COSJ * DJK                                              0803.000
      SK(M) = SK(M) + SMK                                               0804.000
      SMJ = SMI - SMK                                                   0805.000
   60 SJ(M) = SJ(M) + SMJ                                               0806.000
C                                                                       0807.000
C     LOOP OVER THE L-TYPE ATOMS.                                       0808.000
C                                                                       0809.000
      DO 90 N=1,NL                                                      0810.000
      L = LATOM(N)                                                      0811.000
      DLKSQ = ZERO                                                      0812.000
C                                                                       0813.000
      DO 70 M=1,3                                                       0814.000
      T = C(K,M) - C(L,M)                                               0815.000
      ELK(M) = T                                                        0816.000
   70 DLKSQ = DLKSQ + T*T                                               0817.000
C                                                                       0818.000
      DLK = ONE / DSQRT(DLKSQ)                                          0819.000
      COSK = ZERO                                                       0820.000
C                                                                       0821.000
      DO 80 M=1,3                                                       0822.000
      T = ELK(M) * DLK                                                  0823.000
      ELK(M) = T                                                        0824.000
   80 COSK = COSK + EJK(M)*T                                            0825.000
C                                                                       0826.000
      IF (DABS(COSK) .GT. PT9999) GO TO 110                             0827.000
      SIN2K = (ONE-COSK*COSK) * DFLOAT(NL)                              0828.000
      LL = 3 * (L-1)                                                    0829.000
      CALL VPROD (CR, EJK, ELK)                                         0830.000
C                                                                       0831.000
      DO 90 M=1,3                                                       0832.000
      T = CR(M) / SIN2K                                                 0833.000
      SML = T * DLK                                                     0834.000
      IF (DABS(SML) .GE. TENM8) B(NOB,LL+M) = -SML                      0835.000
      SMJ = T * COSK * DJK                                              0836.000
      SJ(M) = SJ(M) + SMJ                                               0837.000
      SMK = SML - SMJ                                                   0838.000
   90 SK(M) = SK(M) + SMK                                               0839.000
C                                                                       0840.000
      DO 100 M=1,3                                                      0841.000
      SMJ = SJ(M)                                                       0842.000
      IF (DABS(SMJ) .GE. TENM8) B(NOB,JJ+M) = SMJ                       0843.000
      SMK = SK(M)                                                       0844.000
  100 IF (DABS(SMK) .GE. TENM8) B(NOB,KK+M) = SMK                       0845.000
C                                                                       0846.000
      RETURN                                                            0847.000
C                                                                       0848.000
C     FATAL ERROR - ATOMS J-K-L ARE COLINEAR.                           0849.000
C                                                                       0850.000
  110 IER = -2                                                          0851.000
      RETURN                                                            0852.000
C                                                                       0853.000
C     ATOMS I-J-K ARE COLINEAR - USE LINEAR BEND.                       0854.000
C                                                                       0855.000
  120 IER = 1                                                           0856.000
      RETURN                                                            0857.000
      END                                                               0858.000
      SUBROUTINE LIBE (NOB, I, J, K, C, A, IC, IER)                     0859.000
C                                                                       0860.000
C     THIS SUBROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A LINEAR BEND  0861.000
C     OR FOR A PAIR OF PERPENDICULAR LINEAR BENDS.                      0862.000
C     I AND K ARE THE END ATOMS.                                        0863.000
C     J IS THE CENTRAL ATOM.                                            0864.000
C                                                                       0865.000
C     A GIVES THE CARTESIAN COORDINATES OF A POINT IN SPACE, SUCH       0866.000
C     THAT THE VECTOR FROM ATOM J TO POINT A IS PERPENDICULAR TO        0867.000
C     THE LINE I-J-K AND SERVES TO ORIENT THE COORDINATES IN SPACE.     0868.000
C     IF IC IS 1, THE IN-PLANE COMPONENT IS RETURNED; IF IC IS 2, THE   0869.000
C     OUT-OF-PLANE (PERPENDICULAR) COMPONENT IS RETURNED.               0870.000
C                                                                       0871.000
C     IER IS A RETURN FLAG:                                             0872.000
C      0 - ROW OF THE B MATRIX WAS COMPUTED SUCCESSFULLY.               0873.000
C     -2 - THE POINT A IS COLINEAR WITH I-J-K.                          0874.000
C                                                                       0875.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               0876.000
C#                                                                      0877.000
C     PARAMETER (NA=#NA)                                                0878.000
C##                                                                     0879.000
      PARAMETER (NA= 30)                                                0880.000
C###                                                                    0881.000
      PARAMETER (NA3=6*NA)                                              0882.000
      PARAMETER (ZERO=0.0D0)                                            0883.000
C                                                                       0884.000
      EXTENDED BLOCK /C711B/ B(NA3,NA3)                                 0885.000
C                                                                       0886.000
      DIMENSION C(NA,3), A(3), EJK(3), UP(3), UN(3)                     0887.000
C                                                                       0888.000
      DATA TENM8/1.0D-8/                                                0889.000
C                                                                       0890.000
      DJISQ = ZERO                                                      0891.000
      DJKSQ = ZERO                                                      0892.000
      DJASQ = ZERO                                                      0893.000
      DUPSQ = ZERO                                                      0894.000
C                                                                       0895.000
      DO 20 M=1,3                                                       0896.000
      TP = C(J,M)                                                       0897.000
      T = C(I,M) - TP                                                   0898.000
      DJISQ = DJISQ + T*T                                               0899.000
      T = C(K,M) - TP                                                   0900.000
      EJK(M) = T                                                        0901.000
      DJKSQ = DJKSQ + T*T                                               0902.000
   20 UN(M) = A(M) - TP                                                 0903.000
C                                                                       0904.000
      DJI = DSQRT(DJISQ)                                                0905.000
      DJK = DSQRT(DJKSQ)                                                0906.000
C     MAKE JA PERPENDICULAR TO I-J-K.                                   0907.000
      CALL VPROD (UP, EJK, UN)                                          0908.000
      CALL VPROD (UN, UP, EJK)                                          0909.000
C                                                                       0910.000
      DO 10 M=1,3                                                       0911.000
      DJASQ = DJASQ + UN(M)*UN(M)                                       0912.000
   10 DUPSQ = DUPSQ + UP(M)*UP(M)                                       0913.000
C                                                                       0914.000
      IF (DJASQ .LE. TENM8) GO TO 60                                    0915.000
      II = 3 * (I-1)                                                    0916.000
      JJ = 3 * (J-1)                                                    0917.000
      KK = 3 * (K-1)                                                    0918.000
C                                                                       0919.000
C     DO IN-PLANE COMPONENT.                                            0920.000
C                                                                       0921.000
      IF (IC .EQ. 2) GO TO 45                                           0922.000
      DJA = DSQRT(DJASQ)                                                0923.000
C                                                                       0924.000
      DO 40 M=1,3                                                       0925.000
      T = -UN(M) / DJA                                                  0926.000
      IF (DABS(T) .LT. TENM8) GO TO 40                                  0927.000
      SMI = T / DJI                                                     0928.000
      B(NOB,II+M) = SMI                                                 0929.000
      SMK = T / DJK                                                     0930.000
      B(NOB,KK+M) = SMK                                                 0931.000
      B(NOB,JJ+M) = -SMI - SMK                                          0932.000
   40 CONTINUE                                                          0933.000
C                                                                       0934.000
      RETURN                                                            0935.000
C                                                                       0936.000
C     DO PERPENDICULAR (OUT-OF-PLANE) COMPONENT.                        0937.000
C                                                                       0938.000
   45 DUP = DSQRT(DUPSQ)                                                0939.000
C                                                                       0940.000
      DO 50 M=1,3                                                       0941.000
      T = -UP(M) / DUP                                                  0942.000
      IF (DABS(T) .LT. TENM8) GO TO 50                                  0943.000
      SMI = T / DJI                                                     0944.000
      B(NOB,II+M) = SMI                                                 0945.000
      SMK = T / DJK                                                     0946.000
      B(NOB,KK+M) = SMK                                                 0947.000
      B(NOB,JJ+M) = -SMI - SMK                                          0948.000
   50 CONTINUE                                                          0949.000
C                                                                       0950.000
      RETURN                                                            0951.000
C                                                                       0952.000
C     A IS COLINEAR WITH I-J-K.                                         0953.000
C                                                                       0954.000
   60 IER = -2                                                          0955.000
      RETURN                                                            0956.000
      END                                                               0957.000
      SUBROUTINE LINTST (I, J, K, C, IER)                               0958.000
C                                                                       0959.000
C     THIS SUBROUTINE TESTS WHETHER THE ANGLE I-J-K IS NEAR LINEAR.     0960.000
C                                                                       0961.000
C     I AND K ARE THE NUMBERS OF THE END ATOMS.                         0962.000
C     J IS THE NUMBER OF THE CENTRAL ATOM.                              0963.000
C                                                                       0964.000
C     IER IS A RETURN FLAG:                                             0965.000
C      0 - ANGLE IS NOT NEAR LINEAR.                                    0966.000
C      1 - ANGLE I-J-K CLOSE TO LINEAR (0 DEGREES).                     0967.000
C     -1 - ANGLE I-J-K CLOSE TO LINEAR (180 DEGREES).                   0968.000
C                                                                       0969.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               0970.000
C#                                                                      0971.000
C     PARAMETER (NA=#NA)                                                0972.000
C##                                                                     0973.000
      PARAMETER (NA= 30)                                                0974.000
C###                                                                    0975.000
      PARAMETER (ZERO=0.0D0)                                            0976.000
C                                                                       0977.000
      DIMENSION C(NA,3), EJI(3), EJK(3)                                 0978.000
C                                                                       0979.000
      DATA PT9999/0.99999D0/                                            0980.000
C                                                                       0981.000
      DJISQ = ZERO                                                      0982.000
      DJKSQ = ZERO                                                      0983.000
C                                                                       0984.000
      DO 20 M=1,3                                                       0985.000
      TP = C(J,M)                                                       0986.000
      T = C(I,M) - TP                                                   0987.000
      EJI(M) = T                                                        0988.000
      DJISQ = DJISQ + T*T                                               0989.000
      T = C(K,M) - TP                                                   0990.000
      EJK(M) = T                                                        0991.000
   20 DJKSQ = DJKSQ + T*T                                               0992.000
C                                                                       0993.000
      DJI = DSQRT(DJISQ)                                                0994.000
      DJK = DSQRT(DJKSQ)                                                0995.000
      DOTJ = ZERO                                                       0996.000
C                                                                       0997.000
      DO 30 M=1,3                                                       0998.000
      T = EJI(M) / DJI                                                  0999.000
      EJI(M) = T                                                        1000.000
      TP = EJK(M) / DJK                                                 1001.000
      EJK(M) = TP                                                       1002.000
   30 DOTJ = DOTJ + T*TP                                                1003.000
C                                                                       1004.000
      IF (DABS(DOTJ) .GT. PT9999) GO TO 60                              1005.000
      IER = 0                                                           1006.000
C                                                                       1007.000
      RETURN                                                            1008.000
C                                                                       1009.000
C     LINEAR BEND AS COS(I-J-K) IS CLOSE TO +/-1.0.                     1010.000
C                                                                       1011.000
   60 IF (DOTJ .GT. ZERO) THEN                                          1012.000
         IER = 1                                                        1013.000
      ELSE                                                              1014.000
         IER = -1                                                       1015.000
      END IF                                                            1016.000
      RETURN                                                            1017.000
      END                                                               1018.000
      SUBROUTINE MATOUTX(X,M,N,MM,NN)                                   1019.000
C                                                                       1020.000
C     MATRIX PRINT ROUTINE                                              1021.000
C                                                                       1022.000
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               1023.000
C                                                                       1024.000
      COMMON/IO/IN,IOUT,IODUM(215)                                      1025.000
C                                                                       1026.000
      DIMENSION X(M,N)                                                  1027.000
C                                                                       1028.000
      EXTENDED DUMMY X                                                  1029.000
C                                                                       1030.000
 1000 FORMAT (3X,10I12)                                                 1031.000
 1001 FORMAT ('0')                                                      1032.000
 1002 FORMAT (1X,I3,2X,10F12.6)                                         1033.000
C                                                                       1034.000
      IFLG=1                                                            1035.000
      ILOWER=1                                                          1036.000
    1 IUPPER=ILOWER+9                                                   1037.000
      IF (IUPPER-NN) 3,2,2                                              1038.000
    2 IUPPER=NN                                                         1039.000
      IFLG=0                                                            1040.000
    3 WRITE(IOUT,1000)(J,J=ILOWER,IUPPER)                               1041.000
      WRITE(IOUT,1001)                                                  1042.000
      DO 4 I=1,MM                                                       1043.000
    4 WRITE(IOUT,1002)I,(X(I,J),J=ILOWER,IUPPER)                        1044.000
      WRITE(IOUT,1001)                                                  1045.000
      IF (IFLG .EQ. 0) RETURN                                           1046.000
      ILOWER=ILOWER+10                                                  1047.000
      GO TO 1                                                           1048.000
      END                                                               1049.000
