      PROGRAM DOC
C
C     PURPOSE:
C       TO GENERATE THE MONSTERGAUSS DOCUMENTATION MANUAL FROM
C       MONSTERGAUSS FORTRAN SOURCE FILES.
C
C     NOTES:
C
C       THE DOCUMENTATION MANUAL WILL BE WRITTEN INTO FILE
C       'MANUAL.OUT' IN THE 'DOC' DIRECTORY. THIS FILE MUST NOT
C       EXIST BEFORE THE 'MANUAL' PROGRAM IS EXECUTED.
C
C       THIS PROGRAM WILL EXTRACT COMMENT LINES FROM FILES OF
C       THE GENERAL FORMAT 'GLNNNN', WITH 'N' AN INTEGER BETWEEN
C       0 AND 9, AND WILL PRINT COLUMNS 3 TO 72 OF THESE COMMENT
C       LINES IN ORDER TO GENERATE THE MANUAL. FILE NAMES OF THE TYPE
C       'GL*.O' ARE AUTOMATICALLY EXCLUDED. EACH FILE WILL
C       BE LOGGED SEQUENTIALLY AND EXTRACTION OF COMMENT LINES WILL
C       BE INITIATED BY THE CODE 'C=HEADER', WITH 'HEADER' THE
C       CURRENT HEADING TO BE USED IN THE MANUAL.  THIS CODE WILL
C       ALSO FORCE A NEW PAGE. IF THE 'HEADER' IS BLANK, THE CURRENT
C       HEADER IS NOT CHANGED, AND TEXT EXTRACTION RESUMES WITH NO JUMP
C       TO A NEW PAGE.
C
C       ALL SUBSEQUENT COMMENT LINE WILL BE EXTRACTED UNTIL A
C       'C==' CODE IS ENCOUNTERED, AND SEARCH IN THE CURRENT FILE
C       WILL BE TERMINATED UPON A 'C/' CODE (THE 'C/' CODE MEANS THERE
C       IS NO MORE DOCUMENTATION IN THIS SOURCE FILE).
C
C       LINES STARTING WITH 'CN' WHERE 'N' IS AN INTEGER BETWEEN '1'
C       AND '9' ARE IGNORED, AS ARE LINES INSIDE THE 'C##'/'C###'
C       REPLACEMENT SECTION OF 'C#'/'C###' CODE BLOCKS.
C
C       LINES CONTAINING 'C CHAPTER ' STARTING IN COLUMN 1 ARE ADDED
C       TO THE TABLE OF CONTENTS.
C
C       TITLE PAGES ARE PREFACED ONTO THE MANUAL - THEY ARE READ FROM
C       LFC 5 UNTIL AN END OF FILE IS ENCOUNTERED. A LINE WITH 'C='
C       IN COLUMNS 1-2 WILL FORCE A NEW PAGE (THE REST OF THE LINE IS
C       IGNORED). NOTE THAT ONLY COLUMNS 3-72 WILL BE COPIED TO THE
C       DOCUMENTATION MANUAL FILE.
C
C       ON THE GOULD, THE DOCUMENTATION MANUAL WILL BE WRITTEN INTO
C       FILE 'MANUAL.OUT' IN THE 'DOC' DIRECTORY. THIS FILE MUST NOT
C       EXIST BEFORE THE 'MANUAL' PROGRAM IS EXECUTED.
C
C     AUTHORS: MAURICE SYLVAIN AND MIKE PETERSON,
C              UNIVERSITY OF TORONTO CHEMISTRY DEPARTMENT,
C              TORONTO, ONTARIO, CANADA  M5S 1A1.
C
C     VERSION: 07 APR 89.
C
      COMMON /A/ ISEQ
      COMMON /B/ FILENAME
C
      INTEGER*4 ISEQ(75)
C?UNIX
      integer chdir, system
c
      character homedir*100
c??
C?GOULD
C     INTEGER*4 SMD4(4), IABCODE/'CRAP'/
C
C     INTEGER*8 IRDBUFF
C
C     CHARACTER DIRECT*6/'^(GAU)'/
C??
      CHARACTER*16 SMD16, FILENAME(100)
      CHARACTER HEAD*70, INLINE*72, LINE*75
C
C?GOULD
C     EQUIVALENCE (SMD16,SMD4(1))
C??
      DATA MAXFILES/75/, HEAD/' '/
C
C?GOULD
C1000 FORMAT ('<'/' @DOCUMENT(SETUP PORTRAIT2)')
C??
 1010 FORMAT (A72)
 1020 FORMAT (1X,A70)
 1030 FORMAT ('1')
 1040 FORMAT (A75)
C
C
C     CREATE OUTPUT FILE MANUAL.OUT.
C
      NFILES=0
C?UNIX
      call getenv ('HOME', homedir)
      nhome = index (homedir, ' ')
      if (nhome .le. 0) go to 9060
      homedir(nhome:nhome) = '/'
      open (unit=2, file=homedir(1:nhome)//'doc/manual.out',
     1 form='formatted', status='new', err=9010, iostat=Istat)
      open (unit=3, form='formatted',
     1 status='scratch', err=9050, iostat=Istat)
c
      istat = chdir (homedir(1:nhome)//'src')
      if (istat .ne. 0) go to 9080
c
c     Log all src files that start with 'glnnnn' where 'n' is a number
c     between 0 and 9.
c     Eliminate files of the form 'gl*.o' automatically.
c
      istat = system ('ls gl[0-9][0-9][0-9][0-9]* >/tmp/mgaussls')
      if (istat .ne. 0) go to 9090
      open (unit=1, file='/tmp/mgaussls', form='formatted',
     1 status='old', err=9070, iostat=istat)
c
  100 read (1,110,end=200) smd16
  110 format (a)
      i = index (smd16, '.o')
      if (i .gt. 0) then
	 if (smd16(i:) .eq. '.o') go to 100
      end if
      nfiles=nfiles + 1
      if (nfiles.gt.maxfiles) go to 9020
      filename(nfiles)=smd16
      iseq(nfiles)=nfiles
      go to 100
c
c     No sort required - files already in alphabetical order.
c
  200 if (nfiles.eq.0) go to 9030
      close (unit=1, status='delete')
C??
C?GOULD
C     OPEN (UNIT=2, BLOCKED=.TRUE., FILE=DIRECT//'MANUAL.OUT',
C    1 FILESIZE=1000, INCREMENT=100, MININCREMENT=50,
C    2 FORM='FORMATTED', STATUS='NEW', ERR=9010, IOSTAT=ISTAT)
C     OPEN (UNIT=3, BLOCKED=.TRUE., FORM='FORMATTED', FILESIZE=1000,
C    1 STATUS='SCRATCH', ERR=9050, IOSTAT=ISTAT)
C
C     LOG ALL GAU FILES THAT START WITH 'GLNNNN' WHERE 'N' IS A NUMBER
C     BETWEEN 0 AND 9.
C
C     CALL X_LOGI (DIRECT, IRDBUFF, 1, SMD4, 4, ISTAT)
C     GO TO 110
C
C 100 CALL X_LOGS (IRDBUFF, 1, SMD4, 4, ISTAT)
C 110 IF (ISTAT) 200,120,9000
C 120 IF (SMD16(1:1).NE.'G' .OR. SMD16(2:2).NE.'L' .OR.
C    1 SMD16(3:3).LT.'0' .OR. SMD16(3:3).GT.'9') GO TO 100
C     IF (SMD16(4:4).LT.'0' .OR. SMD16(4:4).GT.'9' .OR.
C    1 SMD16(5:5).LT.'0' .OR. SMD16(5:5).GT.'9' .OR.
C    2 SMD16(6:6).LT.'0' .OR. SMD16(6:6).GT.'9') GO TO 100
C     NFILES=NFILES + 1
C     IF (NFILES.GT.MAXFILES) GO TO 9020
C     FILENAME(NFILES)=SMD16
C     ISEQ(NFILES)=NFILES
C     GO TO 100
C
C     SORT FILES INTO ALPHABETICAL ORDER.
C
C 200 IF (NFILES.EQ.0) GO TO 9030
C     CALL C16SORT (FILENAME, ISEQ, 1, NFILES)
C??
C
C     INITIALIZE SOME PARAMETERS.
C
      NPAGE = 0
      HEAD = ' '
C?GOULD
C     WRITE(2,1000)
C??
C     COPY TITLE PAGES FROM LFC 5.
C
  300 READ (5,1010,END=600) INLINE
      IF (INLINE(1:2) .EQ. 'C=') THEN
         WRITE (2,1030)
      ELSE
         WRITE (2,1020) INLINE(3:72)
      END IF
      GO TO 300
C
  600 CLOSE (UNIT=5)
C
C     LOOP OVER EACH FILE.
C
      DO 700 I=1,NFILES
      IQ=ISEQ(I)
      SMD16 = FILENAME(IQ)
C?UNIX
      OPEN (UNIT=1, ERR=9040, FILE=HOMEDIR(1:NHOME)//'src/'//SMD16,
     1 STATUS='OLD', IOSTAT=ISTAT, FORM='FORMATTED')
C??
C?GOULD
C     OPEN (UNIT=1, BLOCKED=.TRUE., ERR=9040, FILE=DIRECT//SMD16,
C    1 OPENMODE='R', STATUS='OLD', IOSTAT=ISTAT, FORM='FORMATTED')
C??
      CALL EXTRACT (HEAD, NPAGE, SMD16, *9900)
 700  CLOSE (UNIT=1)
C
C     COPY TEMPORARY FILE TO 'MANUAL.OUT'.
C
      ENDFILE 3
      REWIND 3
C
  800 READ (3,1040,END=900) LINE
      WRITE (2,1040) LINE
      GO TO 800
C
  900 ENDFILE 2
      CLOSE (UNIT=2)
      STOP
C
C     ERROR EXITS.
C
C?GOULD
C9000 WRITE (6,9001) ISTAT
C9001 FORMAT ('0ERROR DURING LOG OF DIRECTORY GAU, ISTAT =',I5)
C     GO TO 9900
C??
 9010 WRITE (6,9011) ISTAT
 9011 FORMAT ('0UNABLE TO CREATE FILE MANUAL.OUT, ISTAT =',I4)
      GO TO 9900
C
 9020 WRITE (6,9021) MAXFILES
 9021 FORMAT ('0TOO MANY SOURCE FILES, LIMIT = MAXFILES =',I4)
      GO TO 9900
C
 9030 WRITE (6,9031)
 9031 FORMAT ('0THERE ARE NO VALID SOURCE FILE NAMES.')
      GO TO 9900
C
 9040 WRITE (6,9041) SMD16,ISTAT
 9041 FORMAT ('0UNABLE TO OPEN FILE ',A16,', ISTAT =',I4)
      GO TO 9900
C
 9050 WRITE (6,9051) ISTAT
 9051 FORMAT ('0UNABLE TO OPEN SCRATCH FILE, ISTAT =',I4)
      GO TO 9900
C?UNIX
 9060 WRITE (6,9061) HOMEDIR
 9061 FORMAT ('0Home directory name is too long:'/1x,A)
      GO TO 9900
C
 9070 WRITE (6,9071) ISTAT
 9071 FORMAT ('0Unable to open file /tmp/mgaussls, istat =', i4)
      GO TO 9900
C
 9080 WRITE (6,9081) ISTAT
 9081 FORMAT ('0Unable to chdir to src, istat =', i4)
      GO TO 9900
C
 9090 WRITE (6,9091) ISTAT
 9091 FORMAT ('0Unable to call system to do the ls, istat =', i4)
C??
C?UNIX
 9900 CALL EXIT (1)
C??
C?GOULD
C9900 IRDBUFF = 0
C     CALL M:ABORT (IRDBUFF, IABCODE)
C??
      STOP
      END
      SUBROUTINE EXTRACT (HEAD, NPAGE, SMD16, *)
C
C     EXTRACT DOCUMENTATION FROM FILE 'SMD16'.
C
      LOGICAL PAGEEND
C
      CHARACTER SMD16*16, HEAD*70, HOLD*70, CPAGE*3
      CHARACTER*72 LINE, TEXT(8)
C
      DIMENSION IB(8)
C
      DATA PAGEEND/.TRUE./
      DATA MAXLINE/58/
C
 1000 FORMAT (A72)
 1010 FORMAT (' ',A70)
 1020 FORMAT ('1'/1X,A70,I4///)
 1030 FORMAT ('1'/1X,70X,I4///)
 1040 FORMAT ('1'/1X,A3,A67///)
 1050 FORMAT ('1'/1X,A3///)
C
C
      MIN = MAXLINE - 5
C
C     FIND A C= CODE IN THE CURRENT SOURCE FILE.
C
   10 READ (1,1000,END=9000) LINE
      IF (LINE(1:1) .NE. 'C') GO TO 10
      IF (LINE(2:2) .EQ. '/') RETURN
      IF (LINE(2:2) .NE. '=') GO TO 10
      IF (LINE(3:3) .EQ. '=') GO TO 9010
C
C     START A NEW MANUAL PAGE AND GET HEADING.
C
   15 I = 73
      HOLD = HEAD
C
   20 I = I - 1
      IF (I .LT. 3) THEN
         IF (PAGEEND) GO TO 50
         GO TO 100
      END IF
      IF (LINE(I:I) .EQ. ' ') GO TO 20
C
      HEAD = ' '
      IJ = (72-I)/2 + 1
      IF (IJ .LE. 4) GO TO 9020
      HEAD(IJ:70) = LINE(3:I)
      IF (HOLD .EQ. HEAD) GO TO 50
      NPAGE = NPAGE + 1
      IF ((NPAGE/2)*2 .EQ. NPAGE) THEN
         WRITE (3,1030) NPAGE
      ELSE
         CALL LEFTPAGE (NPAGE, CPAGE)
         WRITE (3,1050) CPAGE
      END IF
      GO TO 70
C
   50 NPAGE = NPAGE + 1
      IF ((NPAGE/2)*2 .EQ. NPAGE) THEN
         WRITE (3,1020) HEAD, NPAGE
      ELSE
         CALL LEFTPAGE (NPAGE, CPAGE)
         WRITE (3,1040) CPAGE, HEAD(4:70)
      END IF
C
C     FIND A C= OR A C== CODE IN THE CURRENT SOURCE FILE.
C     IGNORE BLANK AND 'CN' LINES AT THE TOP OF A PAGE.
C
   70 READ (1,1000,END=9030) LINE
      IF (LINE(3:72) .EQ. ' ') GO TO 70
      IF (LINE(1:1).EQ.'C' .AND. LINE(2:2).GE.'0' .AND.
     1 LINE(2:2).LE.'9') GO TO  70
      IF (LINE(1:2) .EQ. 'C/') GO TO 9040
      IFIRST = 1
      PAGEEND = .FALSE.
      IF (LINE(1:3) .EQ. 'C==') GO TO 10
      IF (LINE(1:2) .EQ. 'C=') GO TO 15
      IF (LINE(1:3) .NE. 'C##') GO TO 90
   80 READ (1,1000,END=9050) LINE
      IF (LINE(1:4) .NE. 'C###') GO TO 80
      GO TO 70
C
   90 WRITE (3,1010) LINE(3:72)
      IF (LINE(1:10) .EQ. 'C CHAPTER') CALL CHAPTER (LINE, NPAGE)
      IFIRST = 2
C
C     LOOP OVER LINES ON THE CURRENT PAGE.
C
  100 ISTART = IFIRST
      PAGEEND = .FALSE.
C
      DO 131 I=ISTART,MIN
  110 READ (1,1000,END=9030) LINE
      IF (LINE(1:1).EQ.'C' .AND. LINE(2:2).GE.'0' .AND.
     1 LINE(2:2).LE.'9') GO TO 110
      IF (LINE(1:2) .EQ. 'C/') GO TO 9040
      IFIRST = I
      IF (LINE(1:3) .EQ. 'C==') GO TO 10
      IF (LINE(1:2) .EQ. 'C=') GO TO 15
      IF (LINE(1:3) .NE. 'C##') GO TO 130
  120 READ (1,1000,END=9050) LINE
      IF (LINE(1:4) .NE. 'C###') GO TO 120
      LINE(3:72) = ' '
  130 WRITE (3,1010) LINE(3:72)
  131 IF (LINE(1:10) .EQ. 'C CHAPTER') CALL CHAPTER (LINE, NPAGE)
C
C     HANDLE THE LINES AT THE BOTTOM OF THE PAGE SPECIALLY.
C     READ IN TEXT(I) FOR I = 1 TO 8.
C
      DO 160 I=1,8
  135 READ (1,1000,END=9030) TEXT(I)
      IF (TEXT(I)(1:1).EQ.'C' .AND. TEXT(I)(2:2).GE.'0' .AND.
     1 TEXT(I)(2:2).LE.'9') GO TO 135
      IF (TEXT(I)(1:2) .EQ. 'C/') GO TO 9040
      IF (TEXT(I)(1:3) .EQ. 'C==') GO TO 300
      IF (TEXT(I)(1:2) .EQ. 'C=') GO TO 350
      IF (TEXT(I)(1:3) .NE. 'C##') GO TO 150
  140 READ (1,1000,END=9050) TEXT(I)
      IF (TEXT(I)(1:4) .NE. 'C###') GO TO 140
      TEXT(I)(3:72) = ' '
  150 IB(I) = 1
      IF (TEXT(I)(3:72) .EQ. ' ') IB(I) = 0
  160 CONTINUE
C
C     FIND END OF CURRENT PAGE.
C
      IEND = 5
      IF (IB(5).EQ.0 .OR. IB(6).EQ.0) GO TO  200
      IF (IB(4) .EQ. 0) IEND = 4
      IF (IB(3) .EQ. 0) IEND = 3
      IF (IEND .LT. 5) GO TO 200
      IF (IB(8) .EQ. 0) IEND = 4
      IF (IB(7) .NE. 0) GO TO 170
      IEND = 3
      IF (IB(2) .EQ. 0) IEND = 2
      IF (IB(1) .EQ. 0) IEND = 1
      GO TO 200
C
  170 IF (IEND .EQ. 5) GO TO 200
      IF (IB(2) .EQ. 0) IEND = 2
C
C     WRITE END OF CURRENT PAGE.
C
  200 DO 210 I=1,IEND
  210 WRITE (3,1010) TEXT(I)(3:72)
C
C     START NEW PAGE AND SUPPRESS BLANKS AT TOP OF PAGE.
C
  220 IEND = IEND +1
      IF (IEND .GT. 8) GO TO 50
      IF (IB(IEND) .NE. 0) GO TO 230
      GO TO 220
C
  230 NPAGE = NPAGE + 1
      IF ((NPAGE/2)*2 .EQ. NPAGE) THEN
         WRITE (3,1020) HEAD, NPAGE
      ELSE
         CALL LEFTPAGE (NPAGE, CPAGE)
         WRITE (3,1040) CPAGE, HEAD(4:70)
      END IF
      DO 240 I=IEND,8
  240 WRITE (3,1010) TEXT(I)(3:72)
      IFIRST = 10 - IEND
      GO TO 100
C
C     FOUND A C== CODE IN TEXT(I) AT I.
C
  300 IEND = I - 1
      DO 310 J=1,IEND
  310 WRITE (3,1010) TEXT(J)(3:72)
      PAGEEND = .TRUE.
      GO TO 10
C
C     FOUND A C= CODE.
C
  350 IEND = I - 1
      DO 360 J=1,IEND
  360 WRITE (3,1010) TEXT(J)(3:72)
      LINE = TEXT(I)
      PAGEEND = .TRUE.
      GO TO 15
C
C     ERROR SECTION.
C
 9000 WRITE (6,9001) SMD16
 9001 FORMAT ('0*** WARNING ***   C/ CODE MISSING IN FILE ',A16)
      RETURN
C
 9010 WRITE (6,9011) SMD16
 9011 FORMAT ('0*** ERROR ***   FOUND IN FILE ',A16/
     1 ' A C== CODE FOUND BEFORE A C= CODE')
      RETURN 1
C
 9020 WRITE (6,9021) SMD16, LINE
 9021 FORMAT ('0*** ERROR ***   FOUND IN FILE ',A16/
     1 ' A HEADER LINE (C= LINE) IS TOO LONG:'/1X,A72)
      RETURN 1
C
 9030 WRITE (6,9031) SMD16
 9031 FORMAT ('0*** ERROR ***   FOUND IN FILE ',A16/
     1 ' E.O.F. OCCURRED BEFORE A C/ CODE AND BEFORE A C== CODE')
      RETURN 1
C
 9040 WRITE (6,9041) SMD16
 9041 FORMAT ('0*** ERROR ***   FOUND IN FILE ',A16/
     1 ' FOUND A C/ CODE BEFORE A C== CODE')
      RETURN 1
C
 9050 WRITE (6,9051) SMD16
 9051 FORMAT ('0*** ERROR ***   FOUND IN FILE ',A16/
     1 ' MISSING C### CODE')
      RETURN 1
      END
      SUBROUTINE LEFTPAGE (NPAGE, CPAGE)
C
      CHARACTER*3 CPAGE
C
 1000 FORMAT (I1,2X)
 1010 FORMAT (I2,1X)
 1020 FORMAT (I3)
C
C
      IF (NPAGE .LE. 9) THEN
         WRITE (CPAGE,1000) NPAGE
      ELSE IF (NPAGE .LE. 99) THEN
         WRITE (CPAGE,1010) NPAGE
      ELSE
         WRITE (CPAGE,1020) NPAGE
      END IF
      RETURN
      END
      SUBROUTINE CHAPTER (LINE, NPAGE)
C
C     EXTRACT CHAPTER TITLE AND PLACE IN TABLE OF CONTENTS, WITH
C     THE PROPER PAGE NUMBER.
C
      CHARACTER LINE*72, OUTLINE*72
C
 1000 FORMAT ('0',A72)
 1010 FORMAT (I3)
C
C
C     SCAN BACKWARDS TO FIND THE FIRST NON-BLANK CHARACTER IN THE
C     LINE (ALSO IGNORE A TRAILING '.').
C
      N = 73
C
  100 N = N -1
      IF (LINE(N:N).EQ.' ' .OR. LINE(N:N).EQ.'.') GO TO 100
C
      OUTLINE = ' '
      OUTLINE(1:N-2) = LINE(3:N)
      N = N - 1
C
C     INSERT PAGE NUMBER IN COLUMNS 70-72.
C
      WRITE (OUTLINE(70:72), 1010) NPAGE
C
C     INSERT '.' FROM END OF THE CHAPTER TITLE TO THE PAGE NUMBER.
C
  200 N = N + 1
      IF (OUTLINE(N:N) .NE. ' ') GO TO 300
      OUTLINE(N:N) = '.'
      GO TO 200
C
  300 WRITE (2,1000) OUTLINE
      RETURN
      END
