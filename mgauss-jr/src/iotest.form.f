      COMMON /BUFF/ IBUFF(3200)
C
      OPEN (UNIT=1, BLOCKED=.TRUE., FORM='FORMATTED')
C
C     LOOP OVER BUFFERS.
C
      DO 500 I=1,500
C
C     INITIALIZE THE ARRAY.
C
      DO 100 J=1,3200
  100 IBUFF(J) = I + J
C
C     WRITE OUT THE ARRAY.
C
      WRITE (1,200,END=900) IBUFF
  200 FORMAT (16I5)
C
  500 CONTINUE
C
      ENDFILE 1
      STOP
C
  900 STOP 900
      END
