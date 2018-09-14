#
# This file is processed by automake (produces file Makefile.in)
#

# Define the target
bin_PROGRAMS = SetSrc

BUILT_SOURCES = KindType.o Master.o InpOut.o Coordinates.o Plume.o Deposit.o


# Compiler flags
AM_FCFLAGS = $(NC_INC)
#
# Linker flags
# AM_LDFLAGS =
# SetSrc_LDFLAGS = $(AM_LDFLAGS)
#

# Force F90 linker (so it does not get confused with F77)
SetSrc_LINK=$(FCLINK)

#
CLEANFILES = *.mod *~ *__genmod.f90

# List of source files
SetSrc_SOURCES = KindType.f90 Master.f90 InpOut.f90 Coordinates.f90 Plume.f90 \
	Air.f90 Deposit.f90 SetSrc.f90 costa.f90 getnsrc.f90 getsrc.f90 \
	merwind.f90 odepack.f90 openinp.f90 openplumef.f90 opkda1.f opkda2.f \
	opkdmain.f readdbs.f90 readgrn_erup.f90 readgrn_resu.f90 readinp.f90 \
	readinp_erup.f90 readinp_resu.f90 runend.f90 setpsi.f90 solvedust.f90 \
	solveplume.f90 solvepoint.f90 solvesuzuki.f90 wriplumefmass.f90 \
	wriplumeheight.f90 wriplumeprop.f90 wriplumetem.f90 wriresu.f90 \
	wrisrc_mesh.f90 wrisrc_nomesh.f90 wrisrc_resu.f90 wrisrc_void.f90

SetSrc_LDADD = $(top_builddir)/libMaster/libMaster.a $(NC_LIB)
