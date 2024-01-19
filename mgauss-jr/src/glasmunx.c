/*    GLASMUNX     04 APR 91                                         MRP

      PACK/UNPACK for Unix-based systems - this routine is written in c.

      Author: C. G. Marcellus, Guelph-Waterloo Center for Graduate Work
              in Chemistry (Guelph campus), Guelph, Ontario, Canada.
      Version: November 3, 1986.

      IPACK  - takes 4 unsigned 1-byte integers (range 0-255) and packs
               them into a 32-bit integer.

      UNPACK - reverse of IPACK.

      Note: there is no checking for bad input values.

      Note: Since there is no standard for inter-language calls,
      you may have to play with the declared names of the functions
      to be able to call them from FORTRAN. For example, you may need to
      leave off the trailing underscores on the function names, or
      convert the names to upper case, or both. The 'nm' tool is an easy
      way to see what external symbols are defined/required by a routine.
*/
ipack_(i,j,k,l)
int *i,*j,*k,*l;
{
/*
f77: num = ipack(i,j,k,l)
test data: i=3,j=20,k=19,l=100 ==> num=51647332
*/
return((*i & 0xff) << 24 | (*j & 0xff) << 16 | (*k &0xff) <<
8 | (*l & 0xff));
}

unpack_(num,i,j,k,l)
int *i,*j,*k,*l,*num;
{
/*
f77: call unpack(num,i,j,k,l)
test data: num=21345605 ==> i=1,j=69,k=181,l=69
*/
*i = ((*num >> 24) & 0xff); *j = ((*num >> 16) & 0xff);
*k = ((*num >>  8) & 0xff); *l = *num & 0xff;
return;
}
