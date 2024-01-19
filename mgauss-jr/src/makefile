#
# Makefile for Monstergauss on Unix systems.
# ------------------------------------------
#
# Version: 18 Jan 1994.
#
# See the Monstergauss installation manual for detailed setup
# information. Also see the end of this file for special compilation
# rules for many systems.
#
# This makefile supports the following commands:
#
# default: build executable "gausstst" from source files.
# gauss:   build executable "~/bin/gauss" from source files.
# gl*.o:   build object "gl*.o" from source file "gl*.f" or "gl*.c".
# install: install "gausstst" as "~/bin/gauss".
# basis:   build executable "~/bin/basis" from source files.
# compar:  build executable "~/bin/compar" from source files.
# convrt:  build executables "~/bin/convrt1" and "~/bin/convrt2" from
#          source files.
# convrt1: build executable "~/bin/convrt1" from source files.
# convrt2: build executable "~/bin/convrt2" from source files.
# gmerge:  build executable "~/bin/gmerge" from source files.
# gsplit:  build executable "~/bin/gsplit" from source files.
# manual:  build executable "~/bin/manual" from source files.
# zmat:    build executable "~/bin/zmat" from source files.
# clean:   clean up "*.o" and "core".
#
# After "gausstst" has been tested successfully, it is normally
# installed as "~/bin/gauss" with 'make install'.
#
# "Basis" is a program to extract basis sets from the basis set
# tables collected by R. Poirier, R. Kari and I.G. Csizmadia.
#
# "Compar" is a program to compare Monstergauss output files, with
# optional neglect of minor formatting differences.
#
# "Convrt1" and "convrt2" are programs to convert Monstergauss from
# one system to another and/or redimension the program.
#
# "Gmerge" and "gsplit" are programs to split long lines of files
# into lines up to 80 characters long for mailing over BITNET,
# and recombine the lines at the other end.
#
# "Manual" is a program to extract the Monstergauss manual from the
# source code.
#
# "Zmat" is a program to convert Z matrices to Cartesian coordinates.
#
#
# Set compiler options:
# ---------------------
#
# NOTE: Check the bottom of this file for special rules for specific
#       links on the system you are using.
#
# HP systems running HP+UX:
#
# If you have the BSD U77 library (which comes with HP-UX 9.01),
# almost everything is fine. If you don't have the real BSD libraries,
# just comment off the 'CALL FLUSH' line in gl0100c.f, and change
# CALL EXIT(1) in convft2.f and manual.f to 'STOP 1'.
# Note that you should use the '+ppu' option at compile time (set in
# your FCOPTS), or add it here, to access the real BSD libraries.
#
# If you don't have them, the dtime/etime/itime calls in gl0100c.f must
# be commented out. You will also need to comment out each call to
# 'GETENV' in gl0100c.f, and insert below it "FILE = ' '"; note that
# this will cripple the capability to save output from one run to another,
# including final Z matrices. At HP-UX 8.07, you can not optimize
# gl0905.f due to a compiler bug.
#
F77=/opt/gcc2/bin/g77
CC=/opt/gcc2/bin/gcc
CFLAGS= -g3
FFLAGS= -g3
FFLAGS2=
#FFLAGS=+ppu +U77 -O +Obb800
#FFLAGS2=+ppu +U77
#LDFLAGS=+U77
LDFLAGS=
FSPLITFLAGS=
CRITIC=CRITIC

.f.o:
	$(F77) $(FFLAGS) -c $<

#
# Apollo systems running SR10 or higher (set FLOAT to cpu type):
#
#F77=f77
#CFLAGS=-O $$FLOAT
#FFLAGS=-O $$FLOAT
#FFLAGS2=$$FLOAT
#LDFLAGS=-lutchem
#FSPLITFLAGS=-e critic
#CRITIC=*
#
# Sun systems (make sure the FLOAT_OPTION is set correctly):
#
#F77=f77
#CFLAGS=-O
#FFLAGS=-O
#FFLAGS2=
#LDFLAGS=
#FSPLITFLAGS=-e critic
#CRITIC=*
#
# Silicon Graphics / MIPS systems:
#
#F77=f77
#CFLAGS=-O2
#FFLAGS=-O2 -Olimit 1500 -nocpp
#FFLAGS2=-O0 -nocpp
#LDFLAGS=
#FSPLITFLAGS=-e critic
#CRITIC=*
#
# IBM RS/6000 systems:
# Note that the dtime/etime/itime calls in gl0100c.f must be commented
# out, since the standard BSD timing routines are not available.
#
#F77=xlf
#CFLAGS=-O
#FFLAGS=-O -NQ240000 -NT240000 -qextname -yn
#FFLAGS2=-g -NQ240000 -NT240000 -qextname -yn
#LDFLAGS=
#FSPLITFLAGS=-e critic
#CRITIC=*
#
# Stardent systems:
# Note that you will need to modify glasmunx.c to change 'ipack_' to
# 'IPACK' and 'unpack_' to 'UNPACK'.
#
#F77=fc
#CFLAGS=-O
#FFLAGS=-O1
#FFLAGS2=-g
#LDFLAGS=
#FSPLITFLAGS=-e critic
#CRITIC=*
#
OBJS=gl0100a.o gl0100b.o gl0100c.o gl0201a.o\
gl0201b.o gl0202a.o gl0202b.o gl0301a.o\
gl0301b.o gl0301c.o gl0301d.o gl0302.o\
gl0303.o gl0306a.o gl0306b.o gl0306c.o\
gl0307a.o gl0307b.o gl0308.o gl0309.o\
gl0310.o gl0311.o gl0401.o gl0402.o\
gl0403.o gl0501.o gl0502.o gl0503a.o\
gl0503b.o gl0505.o gl0506.o gl0507.o\
gl0601.o gl0602.o gl0603.o gl0701.o\
gl0702.o gl0703.o gl0704.o gl0708.o\
gl0709.o gl0710.o gl0711a.o gl0711b.o\
gl0711c.o gl0801.o gl0802a.o gl0802b.o\
gl0803.o gl0901.o gl0902a.o gl0902b.o\
gl0903.o gl0904.o gl0905.o gl0906.o\
gl0907.o gl9990a.o gl9990b.o gl9990c.o\
gl9990d.o gl9991.o gl9999.o glasmunx.o
#
# "default": make the executable (gausstst).
#            Note that this uses the implicit rules for making
#            *.o files from *.f and *.c files.
#	     Leave the executable in the current directory.
#
gausstst: $(OBJS)
	$(F77) -o gausstst $(OBJS) $(LDFLAGS)
#
# "gauss": make the executable (gauss).
#          Note that this uses the implicit rules for making
#          *.o files from *.f and *.c files.
#    	   Install the executable into ~/bin automatically.
#
gauss: $(OBJS)
	$(F77) -o gauss $(OBJS) $(LDFLAGS)
	mv gauss $$HOME/bin/gauss
#
# "gl9990c.o": special rule (don't optimize BLOCK DATA).
#
gl9990c.o: gl9990c.f
	$(F77) -c $(FFLAGS2) gl9990c.f
#
# "install": install the test version "gausstst" into ~/bin
#	     as the production version "gauss".
#	     The executable will be made if necessary.
#
install: gausstst
	mv gausstst $$HOME/bin/gauss
#
# "basis": make the executable (basis).
#    	   Install the executable into ~/bin automatically.
#
basis:
	$(F77) $(FFLAGS) basis.f -o basis $(LDFLAGS)
	mv basis $$HOME/bin/basis
	rm -f basis.o
#
# "compar": make the executable (compar).
#           Install the executable into ~/bin automatically.
#
compar:
	$(F77) $(FFLAGS) compar.f -o compar $(LDFLAGS)
	mv compar $$HOME/bin/compar
	rm -f compar.o
#
# "convrt": make the executables (convrt1 and convrt2).
#    	    Install the executables into ~/bin automatically.
#
convrt:
	make convrt1
	make convrt2
#
# "convrt1": make the executable (convrt1).
#    	     Install the executable into ~/bin automatically.
#
convrt1:
	$(F77) $(FFLAGS) convrt1.f -o convrt1 $(LDFLAGS)
	mv convrt1 $$HOME/bin/convrt1
	rm -f convrt1.o
#
# "convrt2": make the executable (convrt2).
#            Extract routine 'critic' from gl0100c.f into a
#            temporary directory ~/convrt2tmp, which must not
#            exist before this procedure runs.
#    	     Install the executable into ~/bin automatically.
#
convrt2:
	@echo "Creating special $$HOME/convrt2tmp directory."
	mkdir -p $$HOME/convrt2tmp
	@echo "Extracting critic.f from gl0100c.f."
	#cd $$HOME/convrt2tmp; fsplit $(FSPLITFLAGS) $$HOME/src/gl0100c.f
	$(F77) $(FFLAGS) convrt2.f $$HOME/convrt2tmp/$(CRITIC).f -o convrt2 $(LDFLAGS)
	mv convrt2 $$HOME/bin/convrt2
	rm -f convrt2.o critic.f critic.o
	rm -r $$HOME/convrt2tmp
#
# "gmerge": make the executable (gmerge).
#    	    Install the executable into ~/bin automatically.
#
gmerge:
	$(F77) $(FFLAGS) gmerge.f -o gmerge $(LDFLAGS)
	mv gmerge $$HOME/bin/gmerge
	rm -f gmerge.o
#
# "gsplit": make the executable (gsplit).
#    	    Install the executable into ~/bin automatically.
#
gsplit:
	$(F77) $(FFLAGS) gsplit.f -o gsplit $(LDFLAGS)
	mv gsplit $$HOME/bin/gsplit
	rm -f gsplit.o
#
# "manual": make the executable (manual).
#    	    Install the executable into ~/bin automatically.
#
manual:
	$(F77) $(FFLAGS) manual.f -o manual $(LDFLAGS)
	mv manual $$HOME/bin/manual
	rm -f manual.o
#
# "zmat": make the executable (zmat).
#    	  Install the executable into ~/bin automatically.
#
zmat:
	$(F77) $(FFLAGS) zmat.f -o zmat $(LDFLAGS)
	mv zmat $$HOME/bin/zmat
	rm -f zmat.o
#
# "clean": clean up objects and "core" files.
#
clean:
	-rm -f *.o core $$HOME/bin/core $$HOME/core
#
# ********************************************************************
#
# Start special rules for HP-UX systems:
# --------------------------------------
#
# There are no links that require special treatment with the 9.01
# compiler with patches PHSS_2191 (vec_$copy and blas_$dgemm) and
# PHSS_2723 (fixes for f77 FORTRAN compiler).
#
# End special rules for HP-UX systems.
# ------------------------------------
#
# ********************************************************************
#
# Start special rules for Apollo systems:
# ---------------------------------------
#
# There are no links that require special treatment with the ftn 10.7.p
# (version "364") patch compiler. Monstergauss will not work properly
# with the released ftn 10.8.p compiler (version "188"), cause unknown,
# but does work with the ftn 10.8.p (version "200") patch compiler,
# with no links requiring special treatment.
#
# End special rules for Apollo systems.
# -------------------------------------
#
# ********************************************************************
#
# Start special rules for Sun systems:
# ------------------------------------
#
# "gl0501.o": special rule (SunOS 3.x f77 execution-time failure).
#
#gl0501.o: gl0501.f
#	$(F77) -c $(FFLAGS2) gl0501.f
#
# "gl0502.o": special rule (SunOS 3.x f77 execution-time failure).
#
#gl0502.o: gl0502.f
#	$(F77) -c $(FFLAGS2) gl0502.f
#
# "gl0503a.o": special rule (SunOS 3.x f77 execution-time failure).
#
#gl0503a.o: gl0503a.f
#	$(F77) -c $(FFLAGS2) gl0503a.f
#
# "gl0503b.o": special rule (SunOS 3.x f77 execution-time failure).
#
#gl0503b.o: gl0503b.f
#	$(F77) -c $(FFLAGS2) gl0503b.f
#
# "gl0701.o": special rule (SunOS 3.x f77 execution-time failure).
#
#gl0701.o: gl0701.f
#	$(F77) -c $(FFLAGS2) gl0701.f
#
# "gl0901.o": special rule (SunOS 3.x f77 execution-time failure).
#
#gl0901.o: gl0901.f
#	$(F77) -c $(FFLAGS2) gl0901.f
#
# End special rules for Sun systems.
# ----------------------------------
#
# ********************************************************************
#
# Start special rules for Silicon Graphics / MIPS systems:
# --------------------------------------------------------
#
# "gl0802a.o": special rule (MIPS f77 exection-time failure).
#
#gl0802a.o: gl0802a.f
#	$(F77) -c $(FFLAGS2) gl0802a.f
#
# End special rules for Silicon Graphics / MIPS systems.
# ------------------------------------------------------
#
# ********************************************************************
#
# Start special rules for IBM RS/6000 systems:
# --------------------------------------------
#
# There are no links that require special treatment.
# File gl0100c.f must be modified to deactivate the BSD UNIX timing
# routines that are not supplied on RS/6000 systems.
#
# End special rules for IBM RS/6000 systems.
# ------------------------------------------
#
# ********************************************************************
#
# ********************************************************************
#
# Start special rules for Stardent systems:
# -----------------------------------------
#
# There are no links that require special treatment.
# File glasmunx.c must be modified to change the C function names
# to upper case, and remove the trailing underscores.
#
# End special rules for Stardent systems.
# ---------------------------------------
#
# ********************************************************************
