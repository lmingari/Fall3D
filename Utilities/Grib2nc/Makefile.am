#
# This file is processed by automake (produces file Makefile.in)
#

# Define the target
bin_PROGRAMS = Grib2nc

BUILT_SOURCES = KindType.o Master.o InpOut.o

EXTRA_DIST = config

# Compiler flags
AM_FCFLAGS = $(NC_INC)
#
# Linker flags
# AM_LDFLAGS =
# Grib2nc_LDFLAGS = $(AM_LDFLAGS)
#

# Force F90 linker (so it does not get confused with F77)
Grib2nc_LINK=$(FCLINK)

#
CLEANFILES = *.mod *~ *__genmod.f90

# List of source files
Grib2nc_SOURCES = KindType.f90 Master.f90 InpOut.f90 Grib2nc.f90 inival.f90 \
		openinp.f90 runend.f90 write_nc_grid.f90 write_nc_var2d.f90 \
		write_nc_var3d.f90

Grib2nc_LDADD = $(top_builddir)/libMaster/libMaster.a $(NC_LIB)
