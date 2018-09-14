MODULE parallel
  !***********************************************************************
  !*
  !*    This module is used as an interface between the F90 code and MPI.
  !*    The code is compiled in parallel if the macro WITH_MPI is defined and
  !*    in serial otherwise.
  !*
  !*    LIST OF MODULE VARIABLES
  !*
  !*       nproc    : number of processors            (input)
  !*       ngrp     : number of groups of processors  (input)
  !*
  !*       mpime    : processor number (my task id)   (output)
  !*       root     : root      number (root id)      (output)
  !*       grpnproc : number of processors per group  (output)
  !*       mygrp    : my group number                 (output)
  !*       grpme    : my task id, within MYGRP group  (output)
  !*       grproot  : root group number               (output)
  !*
  !*       group    : father group communicator (MPI_COMM_WORLD)
  !*       grp_intra: child intra group communicator
  !*       grp_inter: child inter group communicator
  !*
  !*    LIST OF ROUTINES IN THE MODULE
  !*
  !*    parallel_startup                        : starts MPI
  !*    parallel_setgroups (ngrp)               : sets the groups of processors
  !*    parallel_hangup                         : closes MPI
  !*    mp_wait( ihand )                        : waits for request
  !*
  !*    bcast        ( array, n, who )                                      : broadcast of an array                                (interfaced)
  !*    parallel_sum ( array, what_group )                                  : performs an array reduction among processors         (interfaced)
  !*    isend        ( array(i), n, idest, itag, ihand, what_group)         : I send    n elements to idest starting from array(i) (interfaced)
  !*    irecv        ( array(i), n, idest, itag, ihand, what_group)         : I recieve n elements to idest starting from array(i) (interfaced)
  !*
  !************************************************************************
  USE KindType
  IMPLICIT NONE
  SAVE
  INTEGER(ip) :: nproc, ngrp
  INTEGER(ip) :: mpime, root, grpnproc, mygrp, grpme, grproot
  INTEGER(ip) :: group, grp_intra, grp_inter
  !
  INTERFACE bcast
     MODULE PROCEDURE bcast_real,bcast_real_1d,bcast_real_2d,bcast_real_3d, &
          bcast_real_4d,bcast_integer,bcast_integer_1d,bcast_integer_2d, &
          bcast_character,bcast_logical
  END INTERFACE
  PRIVATE :: bcast_real,bcast_real_1d,bcast_real_2d,bcast_real_3d, &
       bcast_real_4d, bcast_integer,bcast_integer_1d,bcast_integer_2d, &
       bcast_character,bcast_logical
  !
  INTERFACE parallel_sum
     MODULE PROCEDURE parallel_sum_r1d, parallel_sum_r2d,parallel_sum_r3d, &
          parallel_sum_r4d,parallel_sum_i1d, parallel_sum_i2d, &
          parallel_sum_i3d, parallel_sum_i4d
  END INTERFACE
  PRIVATE ::  parallel_sum_r1d, parallel_sum_r2d,parallel_sum_r3d, &
       parallel_sum_r4d,parallel_sum_i1d, parallel_sum_i2d,parallel_sum_i3d, &
       parallel_sum_i4d
  !
  INTERFACE isend
     MODULE PROCEDURE isend_real,isend_integer
  END INTERFACE
  PRIVATE ::         isend_real,isend_integer
  !
  INTERFACE irecv
     MODULE PROCEDURE irecv_real,irecv_integer
  END INTERFACE
  PRIVATE ::         irecv_real,irecv_integer
  !
  !***********************************************************************
  !*
  !*    Definition of MPI interface routines
  !*
  !***********************************************************************
  !
CONTAINS
  !
  !
  !
  SUBROUTINE parallel_startup
    IMPLICIT NONE
#if defined WITH_MPI
    INCLUDE 'mpif.h'
#endif
    INTEGER(ip) :: ERR
#if defined WITH_MPI
    CALL MPI_INIT(ERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIME,ERR)
    ROOT = 0
    GROUP = MPI_COMM_WORLD
#else
    NPROC = 1
    ROOT  = 0
    MPIME = 0
    GROUP = 0
#endif
    ngrp = 1
    grp_intra = group
    grp_inter = group
    grpme   = mpime
    grproot = root
    grpnproc = nproc / ngrp
    RETURN
  END SUBROUTINE parallel_startup
  !
  !
  !
  SUBROUTINE parallel_setgroups( ngroups )
    IMPLICIT NONE
    INTEGER(ip) :: ngroups
    INTEGER(ip) :: COLOR, KEY, IERR
#if defined WITH_MPI
    INCLUDE 'mpif.h'
#endif
    NGRP = ngroups
    grpnproc = nproc / ngroups
#if defined WITH_MPI
    COLOR = MPIME / GRPNPROC   !  all procs with the same COLOR belong to the same NEWGROUP
    KEY   = MPIME
    CALL MPI_COMM_SPLIT( GROUP, COLOR, KEY, grp_intra, IERR )
    MYGRP = COLOR
    GRPME = MOD( MPIME, GRPNPROC )
    CALL MPI_COMM_SPLIT( GROUP, GRPME, KEY, grp_inter, IERR )
#else
    GRPME = 0
    MYGRP = 0
    grp_intra = GROUP
    grp_inter = GROUP
#endif
    grproot = 0
    RETURN
  END SUBROUTINE parallel_setgroups
  !
  !
  !
  SUBROUTINE MP_WAIT( ihand )
    IMPLICIT NONE
    INTEGER(ip)  :: ihand
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ierr
    INTEGER(ip) :: istatus( MPI_STATUS_SIZE )
    CALL MPI_WAIT( ihand, istatus, ierr )
#endif
    RETURN
  END SUBROUTINE MP_WAIT
  !
  !
  !
  SUBROUTINE parallel_HANGUP
    IMPLICIT NONE
    INTEGER(ip) :: IERR
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    CALL MPI_FINALIZE(IERR)
#endif
    RETURN
  END SUBROUTINE parallel_HANGUP
  !
  !
  !
  SUBROUTINE GLOBALINDEX(LIND,GDIM,NP,ME, gindex)

    !   LIND = INDEX OF AN ELEMENT OF SHEARED ARRAY IN THE LOCAL MEMORY
    !   GDIM = GLOBAL DIMENSION OF AN ARRAY TO BE SHEARED AMONG PROCESSORS
    !   NP   = NUMBER OF PROCESSOR
    !   ME   = INDEX OF THE CALLING PROCESSOR (starting from ZERO)
    !
    !   THIS FUNCTION RETURN THE ABSOLUTE INDEX OF THE LOCAL
    !   ARRAY ELEMENT "LIND"

    IMPLICIT NONE
    INTEGER(ip) :: LIND,GDIM,NP,ME,R,Q, gindex
    Q = INT(GDIM/NP)
    R = MOD(GDIM,NP)
    IF((ME+1).LE.R) THEN
       gindex = (Q+1)*ME + LIND
    ELSE
       gindex = Q*ME + R + LIND
    END IF
    RETURN
  END SUBROUTINE GLOBALINDEX
  !
  !
  !
  SUBROUTINE localdim(gdim,np,me, ldim)

    !   gdim = global dimension of an array to be sheared among processors
    !   np   = number of processor
    !   me   = index of the calling processor (starting from 0)
    !
    !   this function return the number of elements of the array stored
    !   in the local memory of the processor "me"

    IMPLICIT NONE
    INTEGER(ip) :: gdim, np, me, r, q, ldim

    q = INT(gdim / np)
    r = MOD(gdim, np)

    IF( me .LT. r ) THEN
       ldim = q+1
    ELSE
       ldim = q
    END IF

    RETURN
  END SUBROUTINE localdim
  !
  !***********************************************************************
  !*
  !*    LIST OF ROUTINES INTERFACED (module private)
  !*
  !***********************************************************************
  !
  SUBROUTINE BCAST_REAL( ARRAY, N, who )
    IMPLICIT NONE
    REAL(rp)    ::  array
    INTEGER(ip) :: n, who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_REAL
  !
  !
  !
  SUBROUTINE BCAST_REAL_1d( ARRAY, N, who )
    IMPLICIT NONE
    REAL(rp)    ::  array(*)
    INTEGER(ip) :: n, who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_REAL_1d
  !
  !
  !
  SUBROUTINE BCAST_REAL_2d( ARRAY, N, who )
    IMPLICIT NONE
    REAL(rp)    ::  array(:,:)
    INTEGER(ip) :: n, who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_REAL_2d
  !
  !
  !
  SUBROUTINE BCAST_REAL_3d( ARRAY, N, who )
    IMPLICIT NONE
    REAL(rp)    ::  array(:,:,:)
    INTEGER(ip) :: n, who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_REAL_3d
  !
  !
  !
  SUBROUTINE BCAST_REAL_4d( ARRAY, N, who )
    IMPLICIT NONE
    REAL(rp)    ::  array(:,:,:,:)
    INTEGER(ip) :: n, who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_REAL_4d
  !
  !
  !
  SUBROUTINE BCAST_INTEGER(ARRAY,N,who)
    IMPLICIT NONE
    INTEGER(ip) :: ARRAY
    INTEGER(ip) :: N,who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_INTEGER, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_INTEGER
  !
  !
  !
  SUBROUTINE BCAST_INTEGER_1D(ARRAY,N,who)
    IMPLICIT NONE
    INTEGER(ip) :: ARRAY(*)
    INTEGER(ip) :: N,who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_INTEGER, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_INTEGER_1D
  !
  !
  !
  SUBROUTINE BCAST_INTEGER_2D(ARRAY,N,who)
    IMPLICIT NONE
    INTEGER(ip) :: ARRAY(:,:)
    INTEGER(ip) :: N,who
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    CALL MPI_BCAST( array, n, MPI_INTEGER, who, MPI_COMM_WORLD, err)
#endif
    RETURN
  END SUBROUTINE BCAST_INTEGER_2D
  !
  !
  !
  SUBROUTINE BCAST_CHARACTER( ARRAY, N, who )
    IMPLICIT NONE
    CHARACTER(LEN=*) :: ARRAY
    INTEGER(ip)      :: N, who
#if defined WITH_MPI
    INTEGER(ip) :: IERR, I, nn
    INTEGER(ip), ALLOCATABLE :: IARRAY(:)
    nn = LEN( array )
    ALLOCATE( iarray( nn ) )
    DO I=1,nn
       IARRAY(I) = ICHAR( array( i:i ) )
    END DO
    CALL bcast_integer_1d( iarray, nn, who )
    DO I=1,nn
       ARRAY(i:i) = CHAR( iarray( i ) )
    END DO
    DEALLOCATE( iarray )
#endif
    RETURN
  END  SUBROUTINE BCAST_CHARACTER
  !
  !
  !
  SUBROUTINE BCAST_LOGICAL(ARRAY,N,who)
    IMPLICIT NONE
    LOGICAL ARRAY(*)
    INTEGER(ip) :: N
    INTEGER(ip) :: who
#if defined WITH_MPI
    INTEGER(ip) :: IERR, I
    INTEGER(ip), ALLOCATABLE :: IARRAY(:)
    ALLOCATE( iarray( n ) )
    DO I=1,N
       IF(ARRAY(I)) THEN
          IARRAY(I) = 1
       ELSE
          IARRAY(I) = 0
       END IF
    END DO
    CALL bcast_integer_1d( iarray, n, who )
    DO I=1,N
       IF(IARRAY(I).EQ.1) THEN
          ARRAY(I) = .TRUE.
       ELSE
          ARRAY(I) = .FALSE.
       END IF
    END DO
    DEALLOCATE( iarray )
#endif
    RETURN
  END SUBROUTINE BCAST_LOGICAL
  !
  !
  !
  SUBROUTINE ISEND_REAL( array, n, idest, itag, ihand, what_group )
    IMPLICIT NONE
    REAL(rp)    :: array !array(*)
    INTEGER(ip) :: n, idest, itag, ihand, what_group
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ierr
    CALL MPI_ISEND( array, n, MPI_DOUBLE_PRECISION, idest, itag, what_group, ihand, ierr )
#endif
    RETURN
  END SUBROUTINE ISEND_REAL
  !
  !
  !
  SUBROUTINE ISEND_INTEGER( array, n, idest, itag, ihand, what_group )
    IMPLICIT NONE
    INTEGER(ip) :: array !array(*)
    INTEGER(ip) :: n, idest, itag, ihand, what_group
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ierr
    CALL MPI_ISEND( array, n, MPI_INTEGER, idest, itag, what_group, ihand, ierr )
#endif
    RETURN
  END SUBROUTINE ISEND_INTEGER
  !
  !
  !
  SUBROUTINE IRECV_REAL( array, n, isour, itag, ihand, what_group )
    IMPLICIT NONE
    REAL(rp)    :: array !array(*)
    INTEGER(ip) :: n, isour, itag, ihand, what_group
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ierr
    CALL MPI_IRECV( array, n, MPI_DOUBLE_PRECISION, isour, itag, what_group, ihand, ierr )
#endif
    RETURN
  END SUBROUTINE IRECV_REAL
  !
  !
  !
  SUBROUTINE IRECV_INTEGER( array, n, isour, itag, ihand, what_group )
    IMPLICIT NONE
    INTEGER(ip) :: array  !array(*)
    INTEGER(ip) :: n, isour, itag, ihand, what_group
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ierr
    CALL MPI_IRECV( array, n, MPI_INTEGER, isour, itag, what_group, ihand, ierr )
#endif
    RETURN
  END SUBROUTINE IRECV_INTEGER
  !
  !
  !
  SUBROUTINE PARALLEL_SUM_REAL( ARRAY, N, what_grp )
    IMPLICIT NONE
    INTEGER(ip) :: N, what_grp
    REAL(rp)    :: ARRAY(N)
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    REAL(rp), ALLOCATABLE :: ARRAY_TMP(:)
    ALLOCATE( ARRAY_TMP( n ) )
    CALL MPI_ALLREDUCE( ARRAY, ARRAY_TMP, N, MPI_DOUBLE_PRECISION, MPI_SUM,  what_grp, ERR )
    ARRAY(1:n) = ARRAY_TMP(1:n)
    DEALLOCATE( ARRAY_TMP )
#endif
    RETURN
  END SUBROUTINE PARALLEL_SUM_REAL
  !
  !
  !
  SUBROUTINE PARALLEL_SUM_INTEGER(ARRAY, N, what_grp)
    IMPLICIT NONE
    INTEGER N, what_grp
    INTEGER ARRAY(N)
#if defined WITH_MPI
    INCLUDE 'mpif.h'
    INTEGER(ip) :: ib, nb, ipp, nn, err
    INTEGER(ip), ALLOCATABLE :: ARRAY_TMP(:)
    ALLOCATE( ARRAY_TMP( n ) )
    CALL MPI_ALLREDUCE(ARRAY, ARRAY_TMP, N, MPI_INTEGER, MPI_SUM, what_grp, ERR)
    ARRAY(1:n) = ARRAY_TMP(1:n)
    DEALLOCATE( ARRAY_TMP )
#endif
    RETURN
  END SUBROUTINE PARALLEL_SUM_INTEGER
  !
  !
  !
  SUBROUTINE parallel_sum_r1d( array, what_grp )
    REAL   (rp) :: array( : )
    INTEGER(ip) :: what_grp
    CALL parallel_sum_real( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_r1d
  !
  SUBROUTINE parallel_sum_r2d( array, what_grp )
    REAL   (rp) :: array( :, : )
    INTEGER(ip) :: what_grp
    CALL parallel_sum_real( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_r2d
  !
  SUBROUTINE parallel_sum_r3d( array, what_grp )
    REAL   (rp) :: array( :, :, : )
    INTEGER(ip) :: what_grp
    CALL parallel_sum_real( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_r3d
  !
  SUBROUTINE parallel_sum_r4d( array, what_grp )
    REAL   (rp) :: array( :, :, :, : )
    INTEGER(ip) :: what_grp
    CALL parallel_sum_real( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_r4d
  !
  SUBROUTINE parallel_sum_i1d( array, what_grp )
    INTEGER(ip) :: array( : )
    INTEGER(ip) :: what_grp
    CALL parallel_sum_integer( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_i1d
  !
  SUBROUTINE parallel_sum_i2d( array, what_grp )
    INTEGER(ip) :: array( :, : )
    INTEGER(ip) :: what_grp
    CALL parallel_sum_integer( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_i2d
  !
  SUBROUTINE parallel_sum_i3d( array, what_grp )
    INTEGER(ip) :: array( :, :, : )
    INTEGER(ip) :: what_grp
    CALL parallel_sum_integer( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_i3d
  !
  SUBROUTINE parallel_sum_i4d( array, what_grp )
    INTEGER(ip) :: array( :, : ,:, :)
    INTEGER(ip) :: what_grp
    CALL parallel_sum_integer( array, SIZE( array ), what_grp )
    RETURN
  END SUBROUTINE parallel_sum_i4d
  !
  !
END MODULE parallel
