# make.sys.  Generated from make.sys.in by configure.

# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

# most fortran compilers can directly preprocess c-like directives: use
# 	$(MPIF90) $(F90FLAGS) -c $<
# if explicit preprocessing by the C preprocessor is needed, use:
# 	$(CPP) $(CPPFLAGS) $< -o $*.F90 
#	$(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o
# remember the tabulator in the first column !!!

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<

# .f.o and .c.o: do not modify

.f.o:
	$(F77) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS)  -c $<


F90 = ftn
F77 = ftn
CC = cc
MPIF90 = ftn
MPIF77 = ftn
MPICC = cc
CPP = pgCC
LD = ftn -fopenmp
AR = ar
RANLIB = ranlib
 
DFLAGS = -D__GNU -D__PARA -D__MPI -D__FFTW3 -DEXX -D__OPENMP
#DFLAGS = -D__GNU -D__PARA -D__MPI -D__FFTW3 -DEXX
FDFLAGS = $(DFLAGS)
IFLAGS = -I../include -I/opt/intel/composer_xe_2013.1.117/mkl/include/fftw/
MOD_FLAG = -I
 
CFLAGS = -O3 $(DFLAGS) $(IFLAGS)
CPPFLAGS = -E $(DFLAGS) $(IFLAGS)
FFLAGS = -O3 -ffast-math -fdefault-real-8 -fdefault-double-8 -fopenmp
FFLAGS_NOOPT = -O0
F90FLAGS = -cpp -O3 -ffast-math -funroll-loops -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 \
-fopenmp $(FDFLAGS) $(IFLAGS) $(MODFLAGS)
 
LDFLAGS =
ARFLAGS = ruv
 
BLAS_LIBS =
BLAS_LIBS_SWITCH = external
LAPACK_LIBS =
LAPACK_LIBS_SWITCH = external
SCALAPACK_LIBS =
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
FFT_LIBS = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a \
$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a \
$(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -ldl -lpthread -lm
MPI_LIBS =
MASS_LIBS =
PGPLOT_LIBS =
LD_LIBS =
 
FLIB_TARGETS = all
 
LIBOBJS = ../flib/ptools.a ../flib/flib.a ../clib/clib.a \
../iotk/src/libiotk.a
LIBS = $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FFT_LIBS) \
$(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) \
$(PGPLOT_LIBS) $(LD_LIBS)
 
ELPA_LIBS_SWITCH = disabled
 
TOPDIR = /global/homes/b/bsantra/QE-5.0.2-hopper-2013-05-01-GNU
