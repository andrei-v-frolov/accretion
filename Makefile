################################################################
# Makefile for black hole accretion code (Intel Fortran compiler)
################################################################

# Fortran compiler (adjust for your machine)
FC = ifort
FFLAGS = -O3 -ipo -xHOST -heap-arrays 256 -r8 -pc80 -parallel
#LDFLAGS = -static-intel

# CFITSIO libraries
CFITSIO ?= /opt/healpix
FITINC = -I$(CFITSIO)/include
FITLIB = -L$(CFITSIO)/lib -lcfitsio

# LAPACK libraries (use MKL if compiling with Intel Fortran)
MKLROOT ?= /opt/intel/mkl
LAPINC = -I$(MKLROOT)/include
LAPLIB = -L$(MKLROOT)/lib -lmkl_rt

# Intel's dynamic libraries location
#LDFLAGS += -Wl,-rpath,/opt/intel/lib
LDFLAGS += -Wl,-rpath,$(MKLROOT)/lib

INCS += $(FITINC) $(LAPINC)
LIBS += $(FITLIB) $(LAPLIB)


################### Top-level targets ##########################

all: spectral

clean:
	rm -f spectral `find . -name "*.o" -or -name "*.mod" -or -name "*.py[cod]"`
	rm -f *.aux *.blg *.log *.out *.synctex.gz *Notes.bib

################### Binaries & Dependencies ####################

# binaries
spectral: fitsio.o massive.o starobinsky.o hu-sawicki.o

# generic rules
%: %.f90
	$(FC) $(FFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)
%.o %.mod: %.f90
	$(FC) $(FFLAGS) $(INCS) $< -c
