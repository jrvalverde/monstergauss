C/    GL0310       16 MAR 87                                         MRP
C?GLD/GBR/IBM/VAX/UNX
      SUBROUTINE DPHNXW
C??
C?CDC
C     PROGRAM DPHNXW
C??
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C#
C     PARAMETER (NA=#NA)
C##
      PARAMETER (NA= 36)
C###
      COMMON /A/ IOP(99),IC1(NA),IC2(7),C(NA,3),FCON(4),IC3(401)
      COMMON /IO/ IN,IOUT,IODUM(15),NFILE(100,2)
C
 1000 FORMAT ('0ROUTINE DPHNXW NOT AVAILABLE (LINK 0310)')
C
      WRITE (IOUT,1000)
      IOP(1) = -2
      RETURN
      END
