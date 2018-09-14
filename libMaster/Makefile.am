#
# This file is processed by automake (produces file Makefile.in)
#

# Define the target
noinst_LIBRARIES = libMaster.a

# Fortran flags
AM_FCFLAGS = $(NC_INC)

# Compiler flags
AM_LDFLAGS=$(NC_LIB)

# These files need to be generated before all
BUILT_SOURCES = KindType.o Decode.o

CLEANFILES = *.mod *~ *__genmod.f90

# List of source files
libMaster_a_SOURCES = KindType.f90 Decode.f90 get_areas_coordinates.f90 \
	get_dbs_dimension.f90 get_dbs_property_cha.f90 \
	get_dbs_property_int.f90 get_dbs_property_rea.f90 \
	get_dbs_value_dimension.f90 get_dbs_value_plane.f90 \
	get_dbs_value_point.f90 get_dbs_value_volume.f90 \
	get_granulometry_name.f90 get_granulometry_nclass.f90 \
	get_granulometry_value.f90 get_input_cha.f90 get_input_int.f90 \
	get_input_npar.f90 get_input_rea.f90 get_points_coordinates.f90 \
	get_points_npts.f90 get_source_coordinates.f90 get_source_mesh.f90 \
	get_source_nclass.f90 get_source_nsrc.f90 get_source_value.f90
