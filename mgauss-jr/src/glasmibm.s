GLASMIBM TITLE     'IBM ASSEMBLER ROUTINES'
*
*        AUTHOR: MIKE PETERSON, UNIV OF TORONTO CHEMISTRY DEPT,
*                TORONTO, ONTARIO, CANADA.
*
*        UNPACK UNPACKS THE INTEGRAL LABEL JA INTO THE FOUR
*        MOLECULAR ORBITAL INDICES I, J, K AND L. ALL VARIABLES ARE
*        INTEGER*4, AND ARE IN COMMON /PACKED/:
*        COMMON /PACKED/ I, J, K, L, JA
*        EACH INDEX OCCUPIES 7 BITS IN THE LABEL.
*
*        DON'T BOTHER WITH STANDARD SUBROUTINE LINKAGE.
*
UNPACK   START     0
         EXTRN     PACKED
         USING     UNPACK,15
*        SAVE REGISTERS 2 AND 3.
         STM       2,3,28(13)
*        GET ADDRESS OF /PACKED/.
         L         3,=A(PACKED)
*        EXTRACT BYTES TO GET I, J, K AND L.
         SR        2,2
         IC        2,16(3)
         ST        2,0(3)
         IC        2,17(3)
         ST        2,4(3)
         IC        2,18(3)
         ST        2,8(3)
         IC        2,19(3)
         ST        2,12(3)
         LM        2,3,28(13)
         LA        15,0           RETURN CODE OF 0
         BR        14
         END
