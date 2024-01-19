          program test_itime
            integer, dimension(3) :: tarray
            call itime(tarray)
            print *, tarray(1)
            print *, tarray(2)
            print *, tarray(3)
          end program test_itime
