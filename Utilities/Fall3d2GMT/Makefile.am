#
# This file is processed by automake (produces file Makefile.in)
#

# Define the target
bin_PROGRAMS = Fall3d2GMT

BUILT_SOURCES = KindType.o Master.o InpOut.o TimeFun.o

# Compiler flags
AM_FCFLAGS = $(NC_INC)

#
# Linker flags
# AM_LDFLAGS =
# SetSrc_LDFLAGS = $(AM_LDFLAGS)
#

# Force F90 linker
# Fall3d2GMT_LINK=$(FCLINK)

#
CLEANFILES = *.mod *~ *__genmod.f90

# List of source files
Fall3d2GMT_SOURCES = KindType.f90 Master.f90 InpOut.f90 TimeFun.f90 \
	Fall3d2GMT.f90 openinp.f90 readinp.f90 runend.f90 write_script.f90

Fall3d2GMT_LDADD = $(top_builddir)/libMaster/libMaster.a $(NC_LIB)
