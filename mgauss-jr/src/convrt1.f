C
C     PROGRAM CONVRT1
C
C     PURPOSE
C        PREPARE A JOB TO COPY ALL MONSTERGAUSS FILES FROM
C        DIRECTORY 'MGAUSS' TO TAPE.
C
C     AUTHOR: MIKE PETERSON, U OF TORONTO CHEMISTRY DEPARTMENT.
C     VERSION: 28 OCT 90.
C
C     INPUT PARAMETERS:
C     EXCLUDE - (OPTIONAL) LIST OF FILE NAMES, 1 PER LINE, TO BE
C               OMITTED FROM THE COPY.
C               TERMINATE THE LIST OF EXCLUDED FILES WITH A BLANK
C               LINE (THE BLANK LINE IS ALWAYS REQUIRED).
C
C     PNAMES  - PARAMETER NAMES AND VALUES TO BE PASSED TO CONVRT2.
C
C               EACH LINE IS IN THE FORM:
C
C               #ABCDEFG=12345678
C
C               WITH THE '#' CHARACTER IN COLUMN 1, WHERE 'ABCDEFG' IS
C               THE NAME OF THE PARAMETER (LEFT JUSTIFIED AND BLANK
C               FILLED), AND '12345678' IS THE PARAMETER VALUE
C               (I8 FORMAT). THE REMAINDER OF THE INPUT LINE MAY
C               BE USED FOR COMMENTS.
C
C               THE FIRST LINE SHOULD BE ONE OF:
C               #MACH   =GLD      CONVERT TO GOULD SYSTEM (FORT77N)
C               #MACH   =GBR      CONVERT TO GOULD SYSTEM (FORTX32)
C               #MACH   =IBM      CONVERT TO IBM SYSTEM
C               #MACH   =CDC      CONVERT TO CDC SYSTEM
C               #MACH   =VAX      CONVERT TO VAX SYSTEM
C               #MACH   =UNX      CONVERT TO UNIX-BASED SYSTEM
C
C     FOR EXAMPLE, TO CREATE AN IBM VERSION FOR UP TO 150 ATOMIC
C     ORBITALS, THE FOLLOWING WOULD BE SUITABLE:
C
C     #MACH   =IBM      CONVERT TO IBM SYSTEM
C     #NA     =      30 NUMBER OF ATOMS
C     #NB     =     150 NUMBER OF BASIS FUNCTIONS (A.O.'S)
C     #NS     =      80 NUMBER OF SHELLS
C     #NP     =     250 NUMBER OF PRIMITIVE GAUSSIANS
C     #NR     =      76 NUMBER OF OCCUPIED ORBITALS IN RHF
C     #NF     =     100 NUMBER OF FOCK MATRICES IN RHF
C     #NC     =     100 NUMBER OF CONFIGURATIONS IN RHF
C     #NO     =     125 NUMBER OF OPTIMIZABLE PARAMETERS
C     #NT     =      60 NUMBER OF TRANSFORMED ORBITALS
C     #NCON   =   32000 NUMBER OF CONFIGURATIONS IN CI
C     #MA     =      15 NUMBER OF ACTIVE OCCUPIED MO IN CI
C     #MV     =      40 NUMBER OF ACTIVE SUBSTITUTING MO IN CI
C     #NDHF   =     100 NUMBER OF BASIS FUNCTIONS IN DHF
C     #IRECL  =    6232 RECORD SIZE IN BYTES (GOULD=6144, IBM=6232,
C                       REST=6400)
C
C     ARRAY DIMENSIONS:
C     MAXFIL: FILENAME, ISEQ.
C     MAXEXC: EXCLUDE.
C
C     NOTE: NEXCL MUST BE SET TO THE NUMBER OF EXCLUDED FILE NAMES
C     SET IN THE ARRAY EXCLUDE BY THE DATA STATEMENT.
C
C     IN ADDITION, MANY OTHER FILES ARE COPIED, WITHOUT ANY
C     EDITTING. FOR THE COMPLETE LIST, SEE THE 'GAUSS.DOC' FILE.
C
C     *** THIS PROGRAM IS INTENDED FOR UNIX SYSTEMS ONLY ***
C
      PROGRAM CONVRT1
C
      PARAMETER (MAXFIL=150, MAXEXC=20)
C
      COMMON /CON1/ ISEQ, FILENAME
C?UNIX
      integer chdir, system
c
      character direct*5, homedir*100
c
      character*30 smd16, filename(maxfil), exclude(maxexc), c16
c??
C?GOULD
C     INTEGER*8 IRDBUFF
C
C     CHARACTER DIRECT*6/'^(GAU)'/
C
C     CHARACTER*16 SMD16, FILENAME(MAXFIL), EXCLUDE(MAXEXC), C16
C??
      INTEGER SMD4(4), ISEQ(MAXFIL)
C
      CHARACTER*8 PNAME, PVALUE
      CHARACTER*1 G, L, HASH, IEQ, I1, COMMENTS(55)
C
      EQUIVALENCE (SMD16,SMD4(1))
C
      DATA NEXCL/1/, NFILES/0/
      DATA G/'G'/, L/'L'/, HASH/'#'/, IEQ/'='/
C?UNIX
      DATA DIRECT/'/tmp/'/
      DATA EXCLUDE/'glkedgld.tst.job', '                ',
     1 18*' '/
C??
C?GOULD
C     DATA EXCLUDE/'GLKEDGLD.TST    ', '                ',
C    1 18*' '/
C??
C
      WRITE (6,10)
   10 FORMAT ('1FILES TO BE EXCLUDED:'/)
      IF (NEXCL .EQ. 0) GO TO 100
      DO 20 I=1,NEXCL
   20 WRITE (6,30) I, EXCLUDE(I)
   30 FORMAT (1X,I4,3X,A,'   (FILE ALWAYS EXCLUDED BY CONVRT1)')
C
C     READ FILE NAMES TO EXCLUDE.
C
  100 READ (5,110,END=9060) C16
  110 FORMAT (A)
      IF (C16 .EQ. ' ') GO TO 200
      NEXCL = NEXCL + 1
      IF (NEXCL .GT. MAXEXC) GO TO 9010
      WRITE (6,120) NEXCL, C16
  120 FORMAT (1X,I4,3X,A,'   (FILE EXCLUDED BY USER)')
      EXCLUDE(NEXCL) = C16
      GO TO 100
C
C     LOG ALL FILES STARTING WITH 'GL' IN THE SOURCE DIRECTORY.
C
  200 IF (NEXCL .EQ. 0) WRITE (6,210)
  210 FORMAT (' THERE ARE NO EXCLUDED FILES')
C?UNIX
      call getenv ('HOME', homedir)
      nhome = index (homedir, ' ')
      if (nhome .le. 0) go to 9100
      homedir(nhome:nhome) = '/'
c
      istat = chdir (homedir(1:nhome)//'src')
      if (istat .ne. 0) go to 9080
c
c     Log all src files that start with 'gl'.
c     Automatically omit all 'gl*.o' files.
c
      istat = system ('ls gl* >/tmp/mgaussls')
      if (istat .ne. 0) go to 9090
      open (unit=1, file='/tmp/mgaussls', form='formatted',
     1 status='old', err=9070, iostat=istat)
c
  300 read (1,310,end=1000) smd16
  310 format (a)
      i = index (smd16, '.o')
      if (i .gt. 0) then
	 if (smd16(i:) .eq. '.o') go to 300
      end if
      if (nexcl .eq. 0) go to 350
      do 330 i=1,nexcl
      if (smd16 .eq. exclude(i)) go to 300
  330 continue
  350 nfiles = nfiles + 1
      if (nfiles .gt. maxfil) go to 9030
      filename(nfiles) = 'src/' // smd16
      iseq(nfiles) = nfiles
      go to 300
C
C     No sort required - names are already in alphabetical order.
C
 1000 if (nfiles .eq. 0) go to 9040
      close (unit=1, status='delete')
C??
C?GOULD
C     CALL X_LOGI (DIRECT, IRDBUFF, 1, SMD4, 4, ISTAT)
C     GO TO 310
C
C 300 CALL X_LOGS (IRDBUFF, 1, SMD4, 4, ISTAT)
C 310 IF (ISTAT) 1000,320,9020
C 320 IF (SMD16(1:1).NE.G .OR. SMD16(2:2).NE.L) GO TO 300
C     IF (NEXCL .EQ. 0) GO TO 350
C     DO 330 I=1,NEXCL
C     IF (SMD16 .EQ. EXCLUDE(I)) GO TO 300
C 330 CONTINUE
C 350 NFILES = NFILES + 1
C     IF (NFILES .GT. MAXFIL) GO TO 9030
C     FILENAME(NFILES) = SMD16
C     ISEQ(NFILES) = NFILES
C     GO TO 300
C
C     SORT FILE NAMES.
C
C1000 IF (NFILES .EQ. 0) GO TO 9040
C     CALL C16SORT (FILENAME, ISEQ, 1, NFILES)
C??
      WRITE (6,1020)
 1020 FORMAT ('1MONSTERGAUSS FILES TO COPY TO TAPE:'/)
C
      DO 1040 I=1,NFILES
      C16 = FILENAME(ISEQ(I))
 1040 WRITE (6,1050) I, C16
 1050 FORMAT (1X,I4,3X,A)
C
C     APPEND EXTRA FILES.
C
      LAST = NFILES
      LASTP1 = LAST + 1
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'doc/gauss.doc'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'doc/gauss.out'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'src/convrt1.f'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'src/convrt2.f'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'bin/convrt150.job'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'src/manual.f'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'bin/manual.job'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'doc/manual.out'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'src/gmerge.f'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'src/gsplit.f'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'src/makefile'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'doc/gauss.l'
      ISEQ(NFILES) = NFILES
C
      DO 1065 I=1,30
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      WRITE (C16,1060) I
 1060 FORMAT ('test/test',I2.2,'.in')
      FILENAME(NFILES) = C16
 1065 ISEQ(NFILES) = NFILES
C
      DO 1075 I=1,30
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      WRITE (C16,1070) I
 1070 FORMAT ('test/test',I2.2,'.out')
      FILENAME(NFILES) = C16
 1075 ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'gould/s.lib77n'
      ISEQ(NFILES) = NFILES
C
      NFILES = NFILES + 1
      IF (NFILES .GT. MAXFIL) GO TO 9030
      FILENAME(NFILES) = 'gould/j.cl'
      ISEQ(NFILES) = NFILES
C
      WRITE (6,1050) (ISEQ(I),FILENAME(I),I=LASTP1,NFILES)
C
C     CREATE JCL FILE CONVERT.INPUT.
C?UNIX
      OPEN (UNIT=1, ERR=9000, FILE=direct//'convert.input',
     1 IOSTAT=ISTAT, STATUS='NEW')
C??
C?GOULD
C     OPEN (UNIT=1, BLOCKED=.TRUE., CLEAR=.TRUE., ERR=9000,
C    1 FILE=DIRECT//'CONVERT.INPUT', FILESIZE=50, IOSTAT=ISTAT,
C    2 STATUS='NEW')
C??
C     WRITE OUT FILE CONVERT.INPUT.
C
      WRITE (1,1090)
 1090 FORMAT ('#! /bin/csh'/
     1 '#'/
     2 '# Convert.input - copy Monstergauss source to /tmp/newgauss.'/
     3 '#'/
     4 '# Usage: /tmp/convert.input >&/tmp/convert.out &'/
     5 '#'/
     6 'if (! -d /tmp/newgauss) then'/
     7 '   echo " Creating /tmp/newgauss directory."'/
     8 '   mkdir /tmp/newgauss'/
     9 'endif'/
     A 'if (-f /tmp/tape.job) then'/
     B '   echo " Deleting /tmp/tape.job."'/
     C '   rm -f /tmp/tape.job'/
     D 'endif'/
     E '#'/
     F 'echo " "'/
     G '~mgauss/bin/convrt2 <<''EOF''')
C??
C?GOULD
C1090 FORMAT ('$JOB CONVERT.INPUT GAU,1080'/
C    1 '$NOTE COPY MONSTERGAUSS SOFTWARE SOURCE TO TAPE.'/
C    2 '$IFF PATH J.TAPEOUT NODELETE'/
C    3 '$NOTE DELETE (GAU)J.TAPEOUT.'/
C    4 '$DELETE J.TAPEOUT'/
C    5 '$DEFNAME NODELETE'/
C    6 '$NOTE CREATE FILE D.NEWVERSION (IF NECESSARY).'/
C    7 '$IFT PATH D.NEWVERSION SKIP'/
C    8 '$EXEC VOLMGR CRE D.NEWVERSION SIZE=15000 NOSAVE=T SEGN=5'/
C    9 '$DEFNAME SKIP'/
C    A '$NOTE CONVERT GOULD SOURCE TO EXTERNAL SOURCE CODE.'/
C    B '$AS 2 TO D.NEWVERSION BBUF=10'/
C    C '$RUN GAUCONV2')
C??
C
C     COPY PARAMETER DEFINITION RECORDS TO FILE CONVERT.INPUT.
C
      WRITE (6,1100)
 1100 FORMAT ('1PARAMETER DEFINITION RECORDS:'/)
C
 1110 READ (5,1120,END=1200) PNAME, I1, PVALUE, COMMENTS
 1120 FORMAT (A8,A1,A8,55A1)
      WRITE (6,1130) PNAME, I1, PVALUE, COMMENTS
 1130 FORMAT (1X,A8,A1,A8,55A1)
      WRITE (1,1120) PNAME, I1, PVALUE, COMMENTS
      IF (PNAME(1:1).EQ.HASH .AND. I1.EQ.IEQ) GO TO 1110
      WRITE (6,1150)
 1150 FORMAT (' *** WARNING: ABOVE LINE IS NOT A PROPER PARAMETER ',
     1 'DEFINITION LINE ***')
      GO TO 1110
C
C     WRITE OUT THE FILE NAMES TO FILE CONVERT.INPUT.
C
 1200 WRITE (1,1210)
 1210 FORMAT (1X)
C
      DO 1220 I=1,LAST
      C16 = FILENAME(ISEQ(I))
 1220 WRITE (1,1230) C16
 1230 FORMAT (A)
C
C     APPEND EXTRA FILES.
C
      DO 1250 I=LASTP1,NFILES
      C16 = FILENAME(ISEQ(I))
 1250 WRITE (1,1260) C16
 1260 FORMAT (A,T67,'NOSCAN')
C
C     CLOSE FILE CONVERT.INPUT.
C
      WRITE (1,1810)
C?UNIX
 1810 FORMAT ('''EOF'''/
     1 '#'/
     2 'if ($status != 0) then'/
     3 '   exit(1)'/
     4 'endif'/
     5 'echo " "'/
     6 'echo " Run file /tmp/tape.job to continue processing."'/
     7 'chmod +x /tmp/tape.job')
C??
C?GOULD
C1810 FORMAT ('$IFT ABORT EXIT'/
C    1 '$NOTE'/
C    2 '$NOTE RUN FILE (GAU)J.TAPEOUT TO CONTINUE PROCESSING.'/
C    3 '$NOTE'/
C    4 '$EOJ'/
C    5 '$$')
C??
      ENDFILE 1
      CLOSE (UNIT=1)
      WRITE (6,1900) DIRECT
 1900 FORMAT (/'0File ',A,'convert.input ready for ',
     1 'editting/submission.')
      STOP
C
C     ERROR SECTION.
C
 9000 WRITE (6,9001) ISTAT
 9001 FORMAT ('0UNABLE TO CREATE/OPEN FILE CONVERT.INPUT, ISTAT =',I12)
      STOP
C
 9010 WRITE (6,9011) MAXEXC
 9011 FORMAT ('0TOO MANY EXCLUDE FILES, LIMIT =',I5)
      STOP
C?GOULD
C9020 WRITE (6,9021) ISTAT
C9021 FORMAT ('0ERROR DURING LOG OF DIRECTORY GAU, ISTAT =',I5)
C     STOP
C??
 9030 WRITE (6,9031) MAXFIL
 9031 FORMAT ('0TOO MANY SOURCE FILES, LIMIT =',I5)
      STOP
C
 9040 WRITE (6,9041)
 9041 FORMAT ('0THERE ARE NO SOURCE FILES SATISFYING THE MASK ',
     1 'GL*.')
      STOP
C
 9060 WRITE (6,9061)
 9061 FORMAT ('0END-OF-FILE ON INPUT FILE BEFORE PARAMETERS FOUND.')
      STOP
C?UNIX
 9070 WRITE (6,9071) ISTAT
 9071 FORMAT ('0Unable to open file /tmp/mgaussls, istat =', i4)
      STOP
C
 9080 WRITE (6,9081) ISTAT
 9081 FORMAT ('0Unable to chdir to src, istat =', i4)
      STOP
C
 9090 WRITE (6,9091) ISTAT
 9091 FORMAT ('0Unable to call system to do the ls, istat =', i4)
      STOP
C
 9100 WRITE (6,9101) HOMEDIR
 9101 FORMAT ('0Home directory name is too long:'/1x,A)
      STOP
C??
      END
