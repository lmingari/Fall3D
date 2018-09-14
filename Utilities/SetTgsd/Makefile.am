#
# This file is processed by automake (produces file Makefile.in)
#

# Define the target
bin_PROGRAMS = SetTgsd

BUILT_SOURCES = KindType.o InpOut.o Master.o

# Compiler flags
# AM_FCFLAGS =
#
# Linker flags
# AM_LDFLAGS =
# SetTgsd_LDFLAGS = $(AM_LDFLAGS)
#

#
CLEANFILES = *.mod *~ *__genmod.f90

# List of source files
SetTgsd_SOURCES = InpOut.f90 KindType.f90 Master.f90 SetTgsd.f90 openinp.f90 \
	readinp.f90 runend.f90 setfrac.f90 wrigrn.f90

SetTgsd_LDADD = $(top_builddir)/libMaster/libMaster.a
