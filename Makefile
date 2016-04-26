#---------------------------------------------------------------------------
#             AUTHOR TRISTAN SALLES -- UNIVERSITY OF SYDNEY
#---------------------------------------------------------------------------
GFORTRAN    = gfortran
LD = -lm

# OPTIMISATION
#FFLAGS = -O0 -g  -Wall -fbacktrace -lstdc++ -cpp -fcheck=bounds -finit-real=nan \
	-ffpe-trap=zero,overflow,invalid -ffree-form -fno-common -Wtabs -Wunused-parameter \
	-Wuninitialized  -ffree-line-length-none -fdump-fortran-optimized -fdump-tree-original

FFLAGS = -O3

.SUFFIXES : .o .f90

all: orographic

orographic : fourrier.o smithmodel.o ororainmain.o
	$(GFORTRAN) $(FFLAGS) -o $@ ororainmain.o smithmodel.o fourrier.o

.f90.o :
	$(GFORTRAN) $(FFLAGS) -c $(*F).f90

clean :
	/bin/rm -f *.o orographic
