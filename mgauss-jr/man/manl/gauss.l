.TH gauss l "November 28, 1991" "" "U/T Chemistry"
.SH NAME
gauss \- Monstergauss ab initio MO program
.SH SYNOPSIS
gauss
.SH DESCRIPTION
This program performs ab initio molecular orbital calculations using
Gaussian basis sets, and features geometry optimization of minima and
transition structures at the SCF level.
Other capabilities include one electron properties, energy decomposition
(Perturbational MO, and Morokuma energy decomposition), direct CI, and
localized orbitals.
.SH ENVIRONMENT VARIABLES
The following environment variables, which should be set to contain the
appropriate file names (which need not exist), are used if set:
.LP
.RS
GAUSS01        C/D input
.br
GAUSS07        "punch" output
.br
GAUSS08        matrix "scratch" file (for restarts)
.br
GAUSS09        two electron integral file
.br
GAUSS10        transformed integrals file
.br
GAUSS11        CI scratch file
.br
GAUSS12        CI Hamiltonian file
.br
GAUSS13        CI integral labels file
.br
GAUSS14        CI coupling constants file
.br
GAUSS15        CI formula tape file
.br
GAUSS16        two electron derivatives file (not used)
.br
GAUSS19        "setup" debugging output
.RE
.SH EXAMPLES
A basic run:
.RS
gauss  <blah.in  >&blah.out
.RE
.LP
A run saving both the "punch" output and the "scratch" file:
.RS
setenv  GAUSS07  blah.z
.br
setenv  GAUSS08  blah.8
.br
gauss  <blah.in  >&blah.out
.br
unsetenv  GAUSS07  GAUSS08
.RE
Note that the "unsetenv" is important, otherwise the same definitions
will be used for every 'gauss' job submitted during the login session.
.SH NOTES
Full details of the program and a documentation manual may be obtained
from Mike Peterson (mail to "mgauss@alchemy.chem.utoronto.ca").
.br
It is considered good behaviour to "renice" your Monstergauss job to
lower priority after it has been submitted, to avoid the wrath of other
users, and probably the system manager too!
.br
Output files should be printed with 'imagenfl2 blah.out' (or 'imagenf
blah.out' for "paper-waster" mode).
.br
You should not leave your "file 8" files lying around unnecessarily -
they are quite large.
.br
The program will work on the following computers: most UNIX variants,
IBM (MVS and VM/CMS), DEC VAX, CDC (NOS and NOS/VE), and Gould MPX (now
Encore). Porting to another system is easy if the system has 32 bit 
integers, and either 32 or 64 bit reals with 15 decimal digits of
accuracy; others have been done, but with more work.
.SH DIAGNOSTICS
The program attempts to trap all causes of failure internally, so
the fact that the program terminates normally as far as the O/S is
concerned doesn't mean much.
.SH FILES
The following files would of use to potential users:
.LP
.nf
.RS 0
~mgauss/doc/manual.out.Z     Monstergauss manual (print with imagenfp2)
.br
~mgauss/doc/gauss.out.Z      Installation manual (print with imagenp2)
.br
~mgauss/test/testnn.in       Input for test case nn (print with imagenp)
.br
~mgauss/test/testnn.out      Output for test case nn (print with imagenfl2)
.RE
.sp
To print the Monstergauss user manual, do 'cd ~mgauss/doc', then
\'zcat manual.out.Z | imagenfp2 - ' (similarly for the installation
manual, but use imagenp2).
.br
Note that the test input/output files are normally kept compressed in
the '~mgauss/test.tar.Z' file, and the names of the available files can
be viewed with 'zcat ~mgauss/test.tar.Z | tar tvf -' or can
be extracted with 'zcat ~mgauss/test.tar.Z | tar xvf -'
.fi
.SH BUGS
None, I hope!
.SH AUTHORS
Mike Peterson, U/Toronto Chemistry Department.
.SH SEE ALSO
g88(local), tx90(local), txinp(local), mopac(local).
