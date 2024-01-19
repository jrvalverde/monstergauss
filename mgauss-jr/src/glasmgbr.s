*     GLASMGBR     27 AUG 86                                         MRP
*
         PROGRAM   UNPACK
*
         TITLE     'GOULD BASE REGISTER UNPACK ROUTINE'
         LIST      XREF
*
*        AUTHOR: MIKE PETERSON, UNIV OF TORONTO CHEMISTRY DEPT,
*                TORONTO, ONTARIO, CANADA.
*
*        UNPACK UNPACKS THE INTEGRAL LABEL JA INTO THE FOUR
*        MOLECULAR ORBITAL INDICES I, J, K AND L. ALL VARIABLES ARE
*        INTEGER, AND ARE IN COMMON /PACKED/:
*        COMMON /PACKED/ I, J, K, L, JA
*        EACH INDEX OCCUPIES 7 BITS IN THE LABEL.
*
*        THIS ROUTINE IS USED ONLY IN THE GOULD BASE REGISTER FORTRAN
*        VERSION.
*
UNPACK.C DEFSECT   READ                ! DEFINE CODE SECTION (READ-ONLY)
PACKED   DEFSECT   OVERLAY,WORD        ! DEFINE COMMON BLOCK /PACKED/
*
*        DATA SECTION
*
         USESECT   PACKED              ! CREATE COMMON BLOCK /PACKED/
         DEFDATA   PACKED
*
PACKED   DATA      R(5)0
*
*        CODE SECTION
*
         USESECT   UNPACK.C            ! CODE SECTION
         DEFENTRY  UNPACK
         ALIGN     D
         USEBASE   B1,$                ! REGISTER B1 HAS ENTRY ADDRESS
*
UNPACK   LWBR      B7,C.3              ! GET ADDRESS OF /PACKED/
*
         LB        R2,16*B(B7)         ! GET I
         STW       R2,0*W(B7)          ! STORE I
         LB        R2,17*B(B7)         ! GET J
         STW       R2,1*W(B7)          ! STORE J
         LB        R2,18*B(B7)         ! GET K
         STW       R2,2*W(B7)          ! STORE K
         LB        R2,19*B(B7)         ! GET L
         STW       R2,3*W(B7)          ! STORE L
*
*        USE SHORTENED FORTRAN VERSION OF THE CALL/RETURN STANDARD:
*        RELOAD ONLY REGISTER B1 (THE ONLY ONE OF B0, B1, B2 AND B3
*        THAT WAS MODIFIED IN THIS ROUTINE).
*
         LWBR      B1,12(B2)           ! RELOAD CALLERS BASE ADDRESS
         TRSW      R0                  ! RETURN TO CALLER
*
C.3      DATA      A'PACKED'           ! GET ADDRESS OF /PACKED/
         END
