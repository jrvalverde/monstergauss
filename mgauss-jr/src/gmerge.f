c
c     Merge long lines of files back together after mailing over
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
      logical gotsplit
c
      data lline/150/
      data gotsplit/.false./
c
c
  100 read (5,'(a)',end=900) line
      do 110 l=lline,1,-1
      if (line(l:l) .ne. ' ') go to 200
  110 continue
      l = 1
c
c     Handle single lines without change.
c
  200 if (gotsplit) go to 300
      if (line(79:80) .ne. '@@') then
         write (6,'(a)') line(1:l)
         go to 100
      end if
c
c     Save first part of the line.
c
      work(1:78) = line(1:78)
      gotsplit = .true.
      go to 100
c
c     Handle long lines.
c
  300 work(79:) = line(1:l)
      write (6,'(a)') work(1:l+78)
      gotsplit = .false.
      go to 100
c
  900 stop
      end
