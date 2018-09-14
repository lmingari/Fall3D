MODULE Domain
  !**********************************************************************
  !*
  !*    This module contains all the procedures used for the domain
  !*    decomposition and interprocessor data exchange
  !*
  !**********************************************************************
  USE KindType
  USE Parallel
  IMPLICIT NONE
  SAVE
  INTEGER(ip) :: nz_loc  ! local number of z planes
  INTEGER(ip) :: iz_loc  ! global index of the first local plane
  INTEGER(ip) :: nc_loc  ! local number of particle classes
  INTEGER(ip) :: ic_loc  ! global index of the first particles classes
  INTEGER(ip), ALLOCATABLE :: nz_loc_pe(:)
  INTEGER(ip), ALLOCATABLE :: iz_loc_pe(:)
  INTEGER(ip), ALLOCATABLE :: nc_loc_pe(:)
  INTEGER(ip), ALLOCATABLE :: ic_loc_pe(:)
  !
CONTAINS
  !
  !
  !
  SUBROUTINE decomp( nc, nz )
    USE InpOut, ONLY: nlst
    INTEGER(ip) :: nc, nz
    INTEGER(ip) :: iproc
    INTEGER(ip) :: grpme_pe(0:nproc-1),mygrp_pe(0:nproc-1)
    !
    call localdim( nc , ngrp, mygrp, nc_loc )
    call globalindex( 1, nc , ngrp, mygrp, ic_loc )
    call localdim( nz , grpnproc, grpme, nz_loc )
    call globalindex( 1, nz , grpnproc, grpme, iz_loc )
    !
    ! Each processor must have at least two z-layers
    !
    if(nz_loc < 2) call runend('decomp: a minimum of 2 z-layer per processator required')
    !
    ALLOCATE( nz_loc_pe( 0:nproc-1 ) )
    ALLOCATE( iz_loc_pe( 0:nproc-1 ) )
    ALLOCATE( nc_loc_pe( 0:nproc-1 ) )
    ALLOCATE( ic_loc_pe( 0:nproc-1 ) )
    !
    nz_loc_pe = 0
    iz_loc_pe = 0
    nc_loc_pe = 0
    ic_loc_pe = 0
    !
    nz_loc_pe( mpime ) = nz_loc
    iz_loc_pe( mpime ) = iz_loc
    nc_loc_pe( mpime ) = nc_loc
    ic_loc_pe( mpime ) = ic_loc
    !
    CALL parallel_sum( nz_loc_pe, group )
    CALL parallel_sum( iz_loc_pe, group )
    CALL parallel_sum( nc_loc_pe, group )
    CALL parallel_sum( ic_loc_pe, group )
    !
    grpme_pe = 0
    mygrp_pe = 0
    grpme_pe( mpime ) = grpme
    mygrp_pe( mpime ) = mygrp
    CALL parallel_sum( grpme_pe, group )
    CALL parallel_sum( mygrp_pe, group )
    !
    !***  Writes to the log file
    !
    if( mpime == root ) then
       write(nlst,10)
10     format(&
            '---------------------------------------------------',/,  &
            '                                                   ',/,  &
            '         DISTRIBUTION AMONG PROCESSORS             ',/,  &
            '                                                   ',/,  &
            '---------------------------------------------------')
       write(nlst,15)
15     format(/, &
            '  Proc.  Group  Particle-Class  Vertical-Layers    ')
       do iproc = 0, nproc - 1
          write( nlst, 30 ) iproc, mygrp_pe(iproc), ic_loc_pe(iproc), &
               ic_loc_pe(iproc)+nc_loc_pe(iproc)-1,      &
               iz_loc_pe(iproc), iz_loc_pe(iproc)+nz_loc_pe(iproc)-1
       end do
30     format(3x,i4,3x,i4,3x,i4,' to',i4,3x,i4,' to',i4)
    end if
    !
    RETURN
  END SUBROUTINE decomp
  !
  !
  !
  SUBROUTINE collect_cz( c, nz )
    REAL   (rp) :: c( 0:, 0:, 0: )
    INTEGER(ip) :: k, nz
    !
    IF( SIZE( c, 3 ) /= nz+2 ) THEN
       CALL runend('collect_cz : wrong z dimension ' )
    END IF
    do k = 1, nz
       if( ( k < iz_loc ) .or. ( k > (iz_loc + nz_loc - 1) ) ) then
          c( :, :, k ) = 0.0_rp
       end if
    end do
    if( .not. (   iz_loc                ==  1 ) ) c( :, :, 0 )    = 0.0_rp
    if( .not. ( ( iz_loc + nz_loc - 1 ) == nz ) ) c( :, :, nz+1 ) = 0.0_rp
    CALL parallel_sum( c, grp_intra )
    RETURN
  END SUBROUTINE collect_cz
  !
  !
  !
  SUBROUTINE exchange_1layer( c, nz )
    REAL   (rp) :: c( 0:, 0:, 0: )
    INTEGER(ip) :: isour, idest, ishand, irhand
    INTEGER(ip) :: dim
    INTEGER(ip) :: nz
    !
    dim = SIZE( c, 1 ) * SIZE( c, 2 )
    !
    !  Send data to the next processor
    !
    isour = MOD(grpme - 1 + grpnproc, grpnproc)
    idest = MOD(grpme + 1           , grpnproc)

    IF( grpme > 0 ) THEN
       CALL irecv( c(0,0,iz_loc-1), dim, isour, 1, irhand, grp_intra )
    END IF

    IF( grpme < (grpnproc-1) ) THEN
       CALL isend( c(0,0,iz_loc+nz_loc-1), dim, idest, 1, ishand, grp_intra )
    END IF

    IF( grpme > 0 ) THEN
       CALL mp_wait( irhand )
    END IF

    IF( grpme < (grpnproc-1) ) THEN
       CALL mp_wait( ishand )
    END IF

    !
    !  Send data to the previous processor
    !

    isour = MOD(grpme + 1 + grpnproc, grpnproc)
    idest = MOD(grpme - 1           , grpnproc)

    IF( grpme < (grpnproc-1) ) THEN
       CALL irecv( c(0,0,iz_loc+nz_loc), dim, isour, 2, irhand, grp_intra )
    END IF

    IF( grpme > 0 ) THEN
       CALL isend( c(0,0,iz_loc), dim, idest, 2, ishand, grp_intra )
    END IF

    IF( grpme < (grpnproc-1) ) THEN
       CALL mp_wait( irhand )
    END IF

    IF( grpme > 0 ) THEN
       CALL mp_wait( ishand )
    END IF

    RETURN
  END SUBROUTINE exchange_1layer



  SUBROUTINE exchange_2layer( c, nz )
    REAL   (rp) :: c( 0:, 0:, 0: )
    INTEGER(ip) :: isour, idest, ishand, irhand
    INTEGER(ip) :: dim
    INTEGER(ip) :: nz
    !
    dim = SIZE( c, 1 ) * SIZE( c, 2 )
    !
    !  Send data to the next processor
    !
    isour = MOD(grpme - 1 + grpnproc, grpnproc)
    idest = MOD(grpme + 1           , grpnproc)

    IF( grpme > 0 ) THEN
       CALL irecv( c(0,0,iz_loc-2), 2*dim, isour, 1, irhand, grp_intra )
    END IF

    IF( grpme < (grpnproc-1) ) THEN
       CALL isend( c(0,0,iz_loc+nz_loc-2), 2*dim, idest, 1, ishand, grp_intra )
    END IF

    IF( grpme > 0 ) THEN
       CALL mp_wait( irhand )
    END IF

    IF( grpme < (grpnproc-1) ) THEN
       CALL mp_wait( ishand )
    END IF

    !
    !  Send data to the previous processor
    !

    isour = MOD(grpme + 1 + grpnproc, grpnproc)
    idest = MOD(grpme - 1           , grpnproc)

    IF( grpme < (grpnproc-1) ) THEN
       CALL irecv( c(0,0,iz_loc+nz_loc), 2*dim, isour, 2, irhand, grp_intra )
    END IF

    IF( grpme > 0 ) THEN
       CALL isend( c(0,0,iz_loc), 2*dim, idest, 2, ishand, grp_intra )
    END IF

    IF( grpme < (grpnproc-1) ) THEN
       CALL mp_wait( irhand )
    END IF

    IF( grpme > 0 ) THEN
       CALL mp_wait( ishand )
    END IF

    RETURN
  END SUBROUTINE exchange_2layer
  !
  !
  !
  SUBROUTINE collect_c( c, nc )
    REAL   (rp) :: c( :, :, :, : )
    INTEGER(ip) :: ic, nc
    do ic = 1, nc
       if( (ic < ic_loc) .or. ( ic > (ic_loc + nc_loc - 1) ) ) then
          c( :, :, :, ic ) = 0.0_rp
       end if
    end do
    CALL parallel_sum( c, grp_inter )
    RETURN
  END SUBROUTINE collect_c
  !
  !
  !
  SUBROUTINE collect_cload( c, nc ) ! warning. Modify for z decomposition
    REAL   (rp) :: c( :, :, 0: )
    INTEGER(ip) :: ic, nc
    do ic = 1, nc
       if( (ic < ic_loc) .or. ( ic > (ic_loc + nc_loc - 1) ) ) then
          c( :, :, ic ) = 0.0_rp
       end if
    end do
    if( ic_loc > 1 ) then
       c( :, :, 0 ) = 0.0_rp
    end if
    CALL parallel_sum( c, grp_inter )
    RETURN
  END SUBROUTINE collect_cload
  !
  !
  !
END MODULE Domain
