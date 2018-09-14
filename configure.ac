#
# This is the autoconf configuration file for Fall3d
# Process this file with autoconf to produce a configure script.
#
# To produce files: configure, Makefiles and Scripts
#    $ automake
#    $ autoconf
#
# Eventually, if you receive errors, you need to regenerate
# file aclocal.m4 with the command:
#    $ autoreconf
#  or
#    $ aclocal
#    $ autoconf
#
#
# Force an error if autoconf version is earlier than 2.63 (leave commented)
# AC_PREREQ([2.63])
#
# Autoconf initialization [Package name], [version], [email for problems]
AC_INIT([Fall3d], [7.0])
#
# Allows modification of the program file names (keep commented)
#AC_CANONICAL_TARGET
#AC_ARG_PROGRAM
#
# Automake initialization (designed for GNU)
AM_INIT_AUTOMAKE([gnu -Wall -Werror])
#
#
# Set language for configuration checks
AC_LANG(Fortran)
#
# Set the build variables: build_cpu, build_vendor and build_os
# It needs files config.guess and config.sub
AC_CANONICAL_BUILD
#
# Set the AS_PROG_MKDIR_P (usually mkdir -p)
AC_PROG_MKDIR_P
#
# Set default prefix (where directory bin is created)
# This is the default top directory of the distribution (where the
# directory bin is created).
AC_PREFIX_DEFAULT(`pwd`)
#
# Scan for possible compilers (used if environment variable FC is not set)
AC_PROG_FC([ifort,gfortran,xlf90])
# Guess default Fortran linker flags
AC_FC_LIBRARY_LDFLAGS
#
# Assume that F90 is compatible with F77
# Note: F77 is used only in SetSrc to compile odepack
AC_MSG_CHECKING([Fortran 77 compiler])
AC_ARG_VAR([F77],[Fortran 77 compiler command])
AS_IF([test "x$F77" = "x"], [F77=$FC])
AC_MSG_RESULT([$F77])
#
# Define the default location of Run directory
# It defaults to $HOME/Runs (expand the $HOME)
AC_ARG_VAR([RUNDIR],[directory containing Runs])
AS_IF([test "x$RUNDIR" = "x"], [RUNDIR=$HOME/Runs])
AC_MSG_NOTICE([Setting RUNDIR=$RUNDIR])
AC_MSG_NOTICE([Creating RUNDIR directory (if necessary)])
AS_MKDIR_P($RUNDIR)  # Create directory if necessary
#
# Define the default location of data directory
# It defaults to $HOME/Data
AC_ARG_VAR([DATADIR],[directory containing Data])
AS_IF([test "x$DATADIR" = "x"], [DATADIR=$HOME/Data])
AC_MSG_NOTICE([Setting DATADIR=$DATADIR])
AC_MSG_NOTICE([Creating DATADIR directory (if necessary)])
AS_MKDIR_P($DATADIR)  # Create directory if necessary
#
# Location of the scripts (read-only variable)
AC_SUBST([SCRIPTDIR],[$srcdir/Scripts])
#
# *** Check for netcdf ***
#
#  Search sequence:
#  1. Check if configure NETCDF=<directory> is given
#  2. Look for environment variable NETCDF
#  3. Look for nc-config
#  4. Guess NETCDF from the location of ncdump
#
#  If the environmental variable NETCDF is set, NETCDF_INC is set to
#  $NETCDF/include and NETCDF_LIB to $NETCDF/lib.
#  Moreover, both NETCDF_LIB and NETCDF_INC can be set as environmental
#  variables or passed as arguments fo configure; in this case their
#  value overrides the evaluation of $NETCDF/lib or $NETCDF/include
#
#  The compiler and linker flags NC_INC and NC_LIB are set automatically
#  or obtained from nc-config if available, or can be passed to configure.
#
AC_ARG_VAR([NETCDF],[location of the netcdf lib])
AC_ARG_VAR([NC_INC],[compiler flags for netcdf])
AC_ARG_VAR([NC_LIB],[linker flags for netcdf])
if test "X$NETCDF" = "X" ; then
  # NETCDF is not set: look for nc-config
  AC_PATH_PROG([NC_CONF], [nc-config])
  if test -x "$NC_CONF"; then
    AC_MSG_NOTICE([using nc-config])
    AC_MSG_CHECKING(for netCDF includes)
    AS_IF([test "x$NC_INC" = "x"],[NC_INC=`$NC_CONF --fflags`])
    AC_MSG_RESULT($NC_INC)     
    AC_MSG_CHECKING(for netCDF libs)
    AS_IF([test "x$NC_LIB" = "x"],[NC_LIB=`$NC_CONF --flibs`])
    AC_MSG_RESULT($NC_LIB)
    NETCDF=`$NC_CONF --prefix`
    NC_VERSION=`$NC_CONF --version`
    AC_MSG_NOTICE([found $NC_VERSION])
  else
    # Guess location of libnetcdf.a from the directory of ncdump
    AC_PATH_PROG([NETCDF], [ncdump])
    NETCDF=`dirname $NETCDF`   # Strip /ncdump
    NETCDF=`dirname $NETCDF`;  # Strip /bin
    AS_IF([test "x$NETCDF" != "x"],
    [AC_MSG_NOTICE([obtained NETCDF=$NETCDF (from the path of ncdump)])])
  fi
fi
AC_MSG_NOTICE([setting netcdf root directory NETCDF=$NETCDF])
AC_MSG_CHECKING([for the existence of the netcdf root directory])
AS_IF([test -d $NETCDF],[AC_MSG_RESULT([$NETCDF])],
    [AC_MSG_ERROR([directory $NETCDF not found])])
#
# Set the netcdf compiler flags
AC_MSG_CHECKING([for netCDF include directory])
AS_IF([test "x$NETCDF_INC" = "x"],[NETCDF_INC=$NETCDF/include])
AC_MSG_RESULT($NETCDF_INC)
AC_MSG_CHECKING([for the existence of the netcdf include directory])
AS_IF([test -d $NETCDF_INC],[AC_MSG_RESULT([ok])],
    [AC_MSG_ERROR([directory $NETCDF_INC not found])])
AC_MSG_CHECKING([for file netcdf.mod in the include directory])
AS_IF([test -f $NETCDF_INC/netcdf.mod],
   [AC_MSG_RESULT(ok)],
   [AC_MSG_ERROR([netcdf.mod can not be found in $NETCDF_INC])])
AS_IF([ test "x$NC_INC" = "x"],[NC_INC=-I$NETCDF_INC])
#
# Set the netcdf linker flags
AC_MSG_CHECKING([for netCDF library directory])
AS_IF([test "X$NETCDF_LIB" = "X"],[NETCDF_LIB=$NETCDF/lib])
AC_MSG_RESULT([$NETCDF_LIB])
AS_IF([test ! -d $NETCDF_LIB],[AC_MSG_ERROR([No such directory: $NETCDF_LIB])])
AS_IF([test "x$NC_LIB" = "x"],
  [AC_MSG_CHECKING([for -lnetcdff]);
   AS_IF([test -f $NETCDF_LIB/libnetcdff.a],
     [NC_LIB="-L$NETCDF_LIB -lnetcdff -lnetcdf";AC_MSG_RESULT([yes])],
     [NC_LIB="-L$NETCDF_LIB -lnetcdf";AC_MSG_RESULT([no])])])
#
# NC_INC and NC_LIB contain the flags for compilation and linking and are
# exported to Makefiles
AC_SUBST(NC_INC)
AC_SUBST(NC_LIB)
#
#
# *** Check for MPI (if flag --enable-parallel is set)
# Check option "--enable-parallel"
MSG_PARALLEL="(serial)"
AC_ARG_ENABLE([parallel],
  [AS_HELP_STRING([--enable-parallel],[turn on parallel execution (with MPI)])],
  [case "${enableval}" in
    yes) parallel=true ;;
     no) parallel=false ;;
      *) AC_MSG_ERROR([bad value ${enableval} for --enable-parallel]) ;;
   esac], [parallel=false])
AM_CONDITIONAL([PARALLEL], [test "x$parallel" = "xtrue"])
#
# Check the existence of mpif90
AC_ARG_VAR([MPIF90],[parallel compiler (used with --enable-parallel)])
AS_IF([test "x$parallel" = "xtrue"],
  AS_IF([test "x$MPIF90" = "x"],
    [AC_CHECK_PROG([MPIF90], [mpif90], [mpif90])
       AS_IF([test "x$MPIF90" = "x"],
         [AC_MSG_ERROR([cannot find compiler mpif90])])])
  MSG_PARALLEL="(serial and parallel)")
#
# Check the existence of mpiexec
AC_ARG_VAR([MPIEXEC],[run an MPI program (used with --enable-parallel)])
AS_IF([test "x$parallel" = "xtrue"],
  AS_IF([test "x$MPIEXEC" = "x"],
    [AC_CHECK_PROG([MPIEXEC], [mpiexec], [mpiexec])
       AS_IF([test "x$MPIEXEC" = "x"],[MPIEXEC='$MPIEXEC'])]))
#
# Checks for wgrib and wgrib2 (needed by Grib2nc). Warn if not found
AC_ARG_VAR([WGRIB],[Grib files reader (used by the utility Grib2nc)])
AS_IF([test "x$WGRIB" = "x"],
  [AC_CHECK_PROG([WGRIB], [wgrib], [wgrib], [NOTFOUND])
     AS_IF([test "x$WGRIB" = "xNOTFOUND"],
	    [AC_MSG_RESULT([cannot find program wgrib])])])
AC_ARG_VAR([WGRIB2],[Grib2 files reader (used by the utility Grib2nc)])
AS_IF([test "x$WGRIB2" = "x"],
  AC_CHECK_PROG([WGRIB2], [wgrib2], [wgrib2], [NOTFOUND])
     [AS_IF([test "x$WGRIB2" = "xNOTFOUND"],
	    [AC_MSG_RESULT([cannot find program wgrib2])])])
#
# Directory where the Grib2nc configure files are stored
AC_ARG_VAR([GRIBCONFDIR],[Directory where the Grib2nc config files are stored])
AS_IF([test "x$GRIBCONFDIR" = "x"],
	    [GRIBCONFDIR=$ac_pwd/Utilities/Grib2nc/config])
#
# Generate instructions for ranlib
AC_PROG_RANLIB

# List of Makefiles to be processed
AC_CONFIG_FILES([Makefile
		libMaster/Makefile
		Sources/Makefile
		Utilities/SetDbs/Makefile
                Utilities/SetSrc/Makefile
		Utilities/SetTgsd/Makefile
		Utilities/Fall3d2GMT/Makefile
		Utilities/Grib2nc/Makefile
		Utilities/Scripts/Makefile
])

#
# Message to be inserted in the output scripts
AC_SUBST([SCRIPT_WARNING_MESSAGE],
["DO NOT MODIFY (it is generated automatically by configure)"])
#
# List of Scripts to be processed
AC_CONFIG_FILES([
    Scripts/Script-Fall3d_ser:Utilities/Scripts/Script-Fall3d_ser.in 
    Scripts/Script-Fall3d_par:Utilities/Scripts/Script-Fall3d_par.in 
    Scripts/Script-SetDbs:Utilities/Scripts/Script-SetDbs.in 
    Scripts/Script-SetSrc:Utilities/Scripts/Script-SetSrc.in 
    Scripts/Script-SetTgsd:Utilities/Scripts/Script-SetTgsd.in 
    Scripts/Script-Fall3d2GMT:Utilities/Scripts/Script-Fall3d2GMT.in 
    Scripts/Script-arpa-to-nc:Utilities/Scripts/Script-arpa-to-nc.in 
    Scripts/Script-era40-to-nc:Utilities/Scripts/Script-era40-to-nc.in 
    Scripts/Script-eraIn-to-nc:Utilities/Scripts/Script-eraIn-to-nc.in 
    Scripts/Script-eta-to-nc:Utilities/Scripts/Script-eta-to-nc.in 
    Scripts/Script-gfs05deg-to-nc:Utilities/Scripts/Script-gfs05deg-to-nc.in 
    Scripts/Script-gfs1deg-to-nc:Utilities/Scripts/Script-gfs1deg-to-nc.in 
    Scripts/Script-ncep1-to-nc:Utilities/Scripts/Script-ncep1-to-nc.in 
    Scripts/Script-ncep2-to-nc:Utilities/Scripts/Script-ncep2-to-nc.in 
    Scripts/Script-ncepFNL-to-nc:Utilities/Scripts/Script-ncepFNL-to-nc.in 
],[chmod +x Scripts/*])
#
# Produce all the output
AC_OUTPUT
#
# Write configuration on the screen
AC_MSG_NOTICE([---------------------------------------------------------])
AC_MSG_NOTICE([Configuration complete - $PACKAGE_NAME-$PACKAGE_VERSION $MSG_PARALLEL])
AC_MSG_NOTICE([])
AS_IF([test "x$NC_VERSION" != "x"],AC_MSG_NOTICE([Using $NC_VERSION]))
AC_MSG_NOTICE([Fortran 90 compiler:        FC=$FC])
#
AS_IF([test "$parallel" = "true"],
  AC_MSG_NOTICE([Enable parallel version:    --enable-parallel=yes])
  AC_MSG_NOTICE([Fortran 90 PAR compiler:    MPIF90=$MPIF90])
  AC_MSG_NOTICE([Launcher of MPI programs:   MPIEXEC=$MPIEXEC]),
  AC_MSG_NOTICE([Enable parallel version:    --enable-parallel=no])
  AC_MSG_NOTICE([Fortran 90 PAR compiler:    MPIF90=(not used)])
  AC_MSG_NOTICE([Launcher of MPI programs:   MPIEXEC=(not used)]))
#
AC_MSG_NOTICE([Fortran flags:              FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([Fortran 77 compiler:        F77=$F77])
AC_MSG_NOTICE([Root directory of netcdf:   NETCDF=$NETCDF])
AC_MSG_NOTICE([Compiler flags for netcdf:  NC_INC=$NC_INC])
AC_MSG_NOTICE([Linker flags for netcdf:    NC_LIB=$NC_LIB])
AC_MSG_NOTICE([Grib files reader           WGRIB=$WGRIB])
AC_MSG_NOTICE([Grib2 files reader          WGRIB2=$WGRIB2])
AC_MSG_NOTICE([Grib2nc config directory:   GRIBCONFDIR=$GRIBCONFDIR])
AC_MSG_NOTICE([Install prefix:             --prefix=$prefix])
AC_MSG_NOTICE([Executables install prefix: --exec_prefix=$exec_prefix])
AC_MSG_NOTICE([Binary directory:           --bindir=$bindir])
AC_MSG_NOTICE([Run directory:              RUNDIR=$RUNDIR])
AC_MSG_NOTICE([Data directory:             DATADIR=$DATADIR])
AC_MSG_NOTICE([Scripts directory (static): SCRIPTDIR=$SCRIPTDIR])

AC_MSG_NOTICE([---------------------------------------------------------])
