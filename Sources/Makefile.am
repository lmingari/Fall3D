#
# This file is processed by automake (produces file Makefile.in)
#

# Define the target

# Compiler flags
AM_FCFLAGS = $(NC_INC)
#
# Linker flags
#AM_LDFLAGS =

# Modules must be compiled before all
BUILT_SOURCES = KindType.o Numeric.o Master.o InpOut.o Parallel_ser.o \
	Coordinates.o Domain.o

# Serial
bin_PROGRAMS = Fall3d_ser
Fall3d_ser_LINK=$(FCLINK)

# Switch for paraller execution
if PARALLEL
  BUILT_SOURCES += Parallel_par.o
  bin_PROGRAMS += Fall3d_par
  PARCOMPILE=$(MPIF90) $(AM_FCFLAGS) $(FCFLAGS) -DWITH_MPI
  Fall3d_par_LINK=$(MPIF90) $(AM_FCFLAGS) $(FCFLAGS) -o $@ $(AM_LDFLAGS) \
		$(LDFLAGS)
endif

# Distribute Parallel.F90
EXTRA_DIST = Parallel.F90

# Parallel.F90 needs a special rule (it will recompile if Makefile changes)
Parallel_ser.F90: Parallel.F90
	cp Parallel.F90 $@
Parallel_ser.o: Parallel_ser.F90
	$(FCCOMPILE) -c -o $@ $<
#
Parallel_par.F90: Parallel.F90
	cp Parallel.F90 $@
Parallel_par.o: Parallel_par.F90
	$(PARCOMPILE) -c -o $@ $<
#
# Additional files to be cleaned with 'make clean'
CLEANFILES = *.mod *~ Parallel_ser.F90 Parallel_par.F90 *__genmod.f90

# List of source files
COMMON_SOURCES = KindType.f90 Numeric.f90 Master.f90 InpOut.f90 \
	Coordinates.f90 Domain.f90 Fall3d.f90 accum.f90 addradialwind.f90 \
	addtime.f90 advctx.f90 advcty.f90 advctzc.f90 agrsrc.f90 chkdat.f90 \
	cmass2d.f90 cmass3d.f90 corvver.f90 diffx.f90 diffy.f90 diffz.f90 \
	divcorr.f90 endstep.f90 finddt.f90 inidat.f90 kappa3.f90 meteo.f90 \
	printres_nc.f90 readat.f90 reagrd.f90 reagrn.f90 reamet.f90 \
	reaout.f90 rearst.f90 reasrc.f90 runend.f90 setbcc.f90 setcut.f90 \
	setdivu.f90 setpsi.f90 setsrc.f90 setup.f90 setvdrydep.f90 \
	setvset.f90 sizloop.f90 source.f90 ter2asl.f90 wetdep.f90 \
	wriagr.f90 wridat.f90 wrirst.f90 writim.f90 writps.f90 writps_air.f90

Fall3d_ser_SOURCES = $(COMMON_SOURCES)
nodist_Fall3d_ser_SOURCES = Parallel_ser.F90

Fall3d_par_SOURCES = $(COMMON_SOURCES)
nodist_Fall3d_par_SOURCES = Parallel_par.F90

Fall3d_ser_LDADD = $(top_builddir)/libMaster/libMaster.a $(NC_LIB)
Fall3d_par_LDADD = $(top_builddir)/libMaster/libMaster.a $(NC_LIB)
