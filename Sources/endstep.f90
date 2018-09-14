  subroutine endstep
  !***************************************************************************
  !*
  !*    Ends a time step
  !*
  !***************************************************************************
  use KindType
  use InpOut
  use Master
  use Parallel, only: parallel_sum,grp_intra
  use Numeric,  only: nx,ny,nz
  use Domain ,  only: ic_loc, nc_loc, nz_loc, iz_loc, &
                      collect_c,collect_cload,collect_cz,exchange_2layer
  implicit none
  !
  logical          :: tcollect
  integer(ip)      :: ic
  !
  !*** Initializations
  !
  tcollect = .false.
  !
  !*** Prints results
  !
  if( (beg_time + ippfile*print_time) <= time )  then
     ippfile = ippfile+1
     if( .not. tcollect ) then
        do ic = ic_loc, ic_loc + nc_loc - 1
           call collect_cz( c( :, :, :, ic ), nz )  !  return to a serial behaviour
        end do
        call collect_c( c(0:,0:,0:,:), nc )
        call collect_cload( cdry(1:,1:,0:), nc )
        call collect_cload( cwet(1:,1:,0:), nc )
        call getc0
        cload(1:nx,1:ny,0:nc) = cdry(1:nx,1:ny,0:nc) + cwet(1:nx,1:ny,0:nc)
        tcollect = .true.
        !
        if(aggregation) call parallel_sum( Maggr, grp_intra )
        !
     end if
     !
     call printres_nc
     if(track_points) then
        call writps
        call writps_air
     end if
  end if
  !
  !***  Writes log file every 10 minutes
  !
  if(iiter == 1 .or. mod(time,600.0_rp) <= dt) then
     if( .not. tcollect ) then
        do ic = ic_loc, ic_loc + nc_loc - 1
           call collect_cz( c( :, :, :, ic ), nz )  !  return to a serial behaviour
        end do
        call collect_c( c(0:,0:,0:,:), nc )
        call collect_cload( cdry(1:,1:,0:), nc )
        call collect_cload( cwet(1:,1:,0:), nc )
        call getc0
        cload(1:nx,1:ny,0:nc) = cdry(1:nx,1:ny,0:nc) + cwet(1:nx,1:ny,0:nc)
        tcollect = .true.
     end if
     !
     call writim
     !
     !*** Add mass due to the non-null divergence velocity field. This is
     !*** done here and not for each time step in order to avoid excessive
     !*** communication in the parallel version (the correction is based on
     !*** the total mass balance)
     !
     if(gravity_current) call divcorr
     !
  end if
  !
  return
  !
  end subroutine endstep


  subroutine getc0
  !*********************************************************
  !*
  !*  Computes cdry/cwet(i,j,0). Only contributions from
  !*  particles to load.
  !*
  !*********************************************************
  use KindType
  use Master
  use Numeric, ONLY: nx,ny
  implicit none
  !
  integer(ip) :: i,j,ic
  !
  do i = 1,nx
     do j = 1,ny
        cdry(i,j,0) = 0.0_rp
        cwet(i,j,0) = 0.0_rp
        do ic = 1,np  ! particles only
           cdry(i,j,0) = cdry(i,j,0) + cdry(i,j,ic)
           cwet(i,j,0) = cwet(i,j,0) + cwet(i,j,ic)
        end do
     end do
  end do
  !
  return
  end subroutine getc0
