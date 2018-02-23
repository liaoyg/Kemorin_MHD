#
#  Makefile for Kemo's Dynamo toolkit
#    Written by H. Matsui
#
SHELL           = /bin/sh
#
#  directories of Kemo's Dynamo toolkit
#
SRCDIR = /Users/liaoyangguang/GitHub/Kemorin_MHD
INSTDIR= /Users/liaoyangguang/GitHub/Kemorin_MHD
#
MHDDIR = $(SRCDIR)/MHD
MAKEDIR= $(SRCDIR)/work
BUILDDIR= $(SRCDIR)/bin
#
# MPI settings
#
MPICHDIR =    
MPICHLIBDIR = $(MPICHDIR)/lib
MPICHBINDIR = $(MPICHDIR)/bin
MPICHINCDIR = 
MPILIBS =      
#
#    compilers
#
F90 = gfortran
CC =  gcc-5
CXX = g++-5

F90_LOCAL = $(F90)

MPIF90 = mpif90
MPICC =  mpicc

GMAKE =  make
AR =     ar
RANLIB = ranlib -c

#
# OPENMP_FCFLAGS OprnMP option for fortran compiler
# FC_NAME_MANGLE definition for Fortran-C connection
#    C routtine name mangle flag for Fortran connection
#     -DFC_NAME_LOWER_USCORE, -DFC_NAME_UPPER, or -DFC_NAME_LOWER
# FC_NAME_MANGLE -DWITH_FLAG is defined to use cocoa Framework
# XL_FORTRAN   Set "-WF," option for cpp preprocessor
#
OPENMP_FCFLAGS = -fopenmp
FC_NAME_MANGLE = -DFC_NAME_LOWER_USCORE
COCOA_FLAG =     -D__APPLE__
XL_FORTRAN =     
#
#  optimization
#
F90OPTFLAGS= -O3 -Wall -g -ffpe-trap=invalid,zero,overflow,underflow,denormal -ffpe-summary=invalid,zero,overflow,underflow,denormal   -fopenmp
OPTFLAGS =   -O3 -Wall -g    -fopenmp
CXXFLAGS =   -g -O2  -fopenmp
#
#   BLAS settings
#
BLAS_LIBS = -lblas
#
#   libpng and zlib settings
#
ZLIB_CFLAGS = -I/opt/local/include
ZLIB_LIBS =   -L/opt/local/lib -lz
PNG_CFLAGS =  -I/usr/local/Cellar/libpng/1.6.21/include/libpng16
PNG_LIBS =    -L/usr/local/Cellar/libpng/1.6.21/lib -lpng16
#
#    FFTW3 settings
#
FFTW3_CFLAGS= 
FFTW3_LIBS=   
#
#    HDF5 fortran wrapper settings
#
HDF5_FFLAGS=    
HDF5_LDFLAGS=   
HDF5_FLIBS=     
#
#   GLUT settings (Requirement: OpenGL, GLUT, libpng)
#
OPENGL_INC = -I/usr/X11/include
OPENGL_LIBS= -framework GLUT -framework OpenGL  -framework Cocoa
#
#   PGPLOT settings (Requirement: pgplot, libpng)
#
PGPLOT_LIBS= 
#
#   GTK+ settings (Requirement: OpenGL, GLUT, libpng)
GTK2_CFLAGS= 
GTK2_LIBS= 
#
#   GLUI settings (Requirement: OpenGL, GLUT, libpng)
#     It Turns on if C++ compiler is defined
#
#  --- Please do not chenge the following ---
#
MAKE_MOD_DEP= $(BUILDDIR)/make_f90depends
#

all: makemake
	cd $(MAKEDIR); $(GMAKE)

$(MAKE_MOD_DEP): $(MHDDIR)/module_dependency/make_module_dependency.f90
	if [ ! -d $(MAKEDIR) ]; then \
		mkdir $(MAKEDIR); \
	fi
	if [ ! -d $(BUILDDIR) ]; then \
		mkdir $(BUILDDIR); \
	fi
	$(F90_LOCAL) -c -o $(MAKEDIR)/make_module_dependency.o $<
	$(F90_LOCAL) -o $@ $(MAKEDIR)/make_module_dependency.o

makemake: $(MAKE_MOD_DEP)
	echo "# Construct Makefile"; \
	cd $(MHDDIR) ; \
		$(GMAKE) \
		GMAKE="$(GMAKE)" \
		SRCDIR="$(SRCDIR)" \
		INSTDIR="$(INSTDIR)" \
		MAKEDIR="$(MAKEDIR)" \
		BUILDDIR="$(BUILDDIR)" \
		MHDDIR="$(MHDDIR)" \
		MPICHDIR="$(MPICHDIR)" \
		MPICHLIBDIR="$(MPICHLIBDIR)" \
		MPILIBS="$(MPILIBS)" \
		MPICHBINDIR="$(MPICHBINDIR)" \
		MPICHINCDIR="$(MPICHINCDIR)" \
		BLAS_LIBS="$(BLAS_LIBS)" \
		ZLIB_CFLAGS="$(ZLIB_CFLAGS)" \
		ZLIB_LIBS="$(ZLIB_LIBS)"     \
		PNG_CFLAGS="$(PNG_CFLAGS)"   \
		PNG_LIBS="$(PNG_LIBS)"       \
		FFTW3_CFLAGS="$(FFTW3_CFLAGS)" \
		FFTW3_LIBS="$(FFTW3_LIBS)" \
		OPENGL_INC="$(OPENGL_INC)" \
		OPENGL_LIBS="$(OPENGL_LIBS)" \
		PGPLOT_LIBS="$(PGPLOT_LIBS)" \
		GTK2_CFLAGS="$(GTK2_CFLAGS)" \
		GTK2_LIBS="$(GTK2_LIBS)" \
		HDF5_FFLAGS="$(HDF5_FFLAGS)" \
		HDF5_LDFLAGS="$(HDF5_LDFLAGS)" \
		HDF5_FLIBS="$(HDF5_FLIBS)" \
		FC_NAME_MANGLE="$(FC_NAME_MANGLE)" \
		COCOA_FLAG="$(COCOA_FLAG)" \
		OPTFLAGS="$(OPTFLAGS)" \
		F90OPTFLAGS="$(F90OPTFLAGS)" \
		CXXFLAGS="$(CXXFLAGS)" \
		F90LIB="$(F90LIB)" \
		CC="$(CC)" \
		CXX="$(CXX)" \
		MPIF90="$(MPIF90)" \
		MPICC="$(MPICC)"   \
		AR="$(AR)" \
		RANLIB="$(RANLIB)" \
		MAKE_MOD_DEP="$(MAKE_MOD_DEP)" \
        XL_FORTRAN="$(XL_FORTRAN)" \
		makemake


install:
	cd $(MAKEDIR) ; \
		$(GMAKE) install

clean:
	for dir in $(MAKEDIR) $(MHDDIR) ; do \
	echo "# cleaning $${dir} directory..."; \
		( cd $${dir}; \
		$(GMAKE) clean )\
	done; \
	rm -f mpif.h *.o *.mod *~ *.par *.diag *.a *.f90

distclean:
	echo "# Back to initial package"; \
	rm -fr $(MAKEDIR) $(BUILDDIR) Makefile config.log
