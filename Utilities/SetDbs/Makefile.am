#
# This file is processed by automake (produces file Makefile.in)
#

# Define the target
bin_PROGRAMS = SetDbs

# These are modules
BUILT_SOURCES = KindType.o InpOut.o Master.o Coordinates.o \
	TimeFun.o ARPA_nc.o CAL_nc.o DBS_nc.o ECMWF_nc.o \
	ETA_nc.o GFS_nc.o  Wind_Rotation.o MathFun.o NCEP1_nc.o NCEP2_nc.o \
	PROF_nc.o WRF_nc.o

# Compiler flags
AM_FCFLAGS = $(NC_INC)

#
# Linker flags
# AM_LDFLAGS=
# SetDbs_LDFLAGS = $(AM_LDFLAGS)

#
CLEANFILES = *.mod *~ *__genmod.f90

# List of source files
SetDbs_SOURCES = KindType.f90 InpOut.f90 Master.f90 Coordinates.f90 \
	TimeFun.f90 ARPA_nc.f90 CAL_nc.f90 DBS_nc.f90 ECMWF_nc.f90 \
	ETA_nc.f90 GFS_nc.f90 Wind_Rotation.f90 MathFun.f90 NCEP1_nc.f90 \
	NCEP2_nc.f90 PROF_nc.f90 WRF_nc.f90 SetDbs.f90 checktime_ARPA.f90 \
	checktime_CAL.f90 checktime_ECMWF.f90 checktime_ETA.f90 \
	checktime_GFS.f90 checktime_NCEP1.f90 checktime_NCEP2.f90 \
	checktime_PROF.f90 checktime_WRF.f90 get_par_ABL.f90 openinp.f90 \
	read_ARPA_data.f90 read_ARPA_grid.f90 read_CAL_data.f90 \
	read_CAL_grid.f90 read_DEBUG_data.f90 read_ECMWF_data.f90 \
	read_ECMWF_grid.f90 read_ETA_data.f90 read_ETA_grid.f90 \
	read_GFS_data.f90 read_GFS_grid.f90 read_NCEP1_data.f90 \
	read_NCEP1_grid.f90 read_NCEP2_data.f90 read_NCEP2_grid.f90 \
	read_PROF_data.f90 read_PROF_grid.f90 read_WRF_data.f90 \
	read_WRF_grid.f90 readinp.f90 readtop.f90 runend.f90 \
	write_DBS_data.f90 write_DBS_grid.f90 writim.f90

SetDbs_LDADD = $(top_builddir)/libMaster/libMaster.a $(NC_LIB)
