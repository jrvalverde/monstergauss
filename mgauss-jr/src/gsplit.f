c
c     Split long lines of files into 2 lines for easy mailing over
c     BITNET.
c
c     This program acts like a UNIX filter, reading standard input
c     from FORTRAN unit 5, and writing standard output to FORTRAN
c     unit 6.
c
c     Define a string to hold the input line and a work space.
c
      character line*150, work*150
c
      data lline/150/
c
c
  100 read (5,'(a)',end=900) line
      do 110 l=lline,1,-1
      if (line(l:l) .ne. ' ') go to 200
  110 continue
      l = 1
c
c     Handle lines up to 80 characters without change.
c
  200 if (l .le. 80) then
         write (6,'(a)') line(1:l)
         go to 100
      end if
c
c     Handle long lines.
c
      work = line(79:l)
      line(79:80) = '@@'
      write (6,'(a)') line(1:80)
      write (6,'(a)') work(1:l-78)
      go to 100
c
  900 stop
      end
