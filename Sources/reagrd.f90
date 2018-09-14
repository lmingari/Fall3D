subroutine reagrd
  !****************************************************************
  !*
  !*    Reads grid data from the dbs file
  !*
  !****************************************************************
  use KindType
  use InpOut
  use Numeric
  use Master
  use Parallel
  implicit none
  !
  integer(ip)           :: istat,ivoid,j,k
  real   (rp)           :: lat,colat
  character(len=10    ) :: cvoid
  character(len=s_mess) :: message
  !
  !***  Reads system of coordinates (and checks for consistency with that of the dbs)
  !
  coord_sys(:) = ' '
  if( mpime == root ) then
     call get_dbs_property_cha(fdbs,'COORDINATES',coord_sys,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( coord_sys, LEN( coord_sys ), root )
  !
  cvoid(:) = ' '
  if( mpime == root ) then
     call get_input_cha (finp,'GRID','COORDINATES',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( cvoid, LEN( cvoid ), root )
  !
  if(TRIM(coord_sys) /= TRIM(cvoid)) &
       call runend('System of coordinates in dbs and input file do not coincide')
  !
  !*** Reads data
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'YEAR',idbsyr,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( idbsyr, 1, root )
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'MONTH',idbsmo,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( idbsmo, 1, root )
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'DAY',idbsdy,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( idbsdy, 1, root )
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'START_TIME',ibdbs,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( ibdbs, 1, root )
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'END_TIME',iedbs,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( iedbs, 1, root )
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'NX',nx,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast(    nx, 1, root )
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'NY',ny,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast(    ny, 1, root )
  !
  if( mpime == root ) then
     call get_dbs_property_int(fdbs,'NZ',nz,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast(    nz, 1, root )
  !
  !*** Calculates the time_lag (if any) between dbs and input data
  !
  call get_time_lag
  !
  !***  Allocates memory for arrays related to nx,ny,nz
  !
  call getmem1
  !
  !***  Horizontal dimensions and scaling factors
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'LONMIN',lonmin,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(lonmin, 1, root )
     !
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'LATMIN',latmin,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(latmin, 1, root )
     !
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'LONMAX',lonmax,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(lonmax, 1, root )
     !
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'LATMAX',latmax,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(latmax, 1, root )
     !
     xorigr = lonmin
     yorigr = latmin
     dx = (lonmax-lonmin)/(nx-1)
     dy = (latmax-latmin)/(ny-1)
     !
     dX1 = Rearth*(dx*pi/180.0_rp)  ! Scaled (dx,dy)
     dX2 = Rearth*(dy*pi/180.0_rp)
     !
  case('UTM')
     !
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'XMIN',xmin,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(xmin, 1, root )
     !
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'YMIN',ymin,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(ymin, 1, root )
     !
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'XMAX',xmax,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(xmax, 1, root )
     !
     if( mpime == root ) then
        call get_dbs_property_rea(fdbs,'YMAX',ymax,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast(ymax, 1, root )
     !
     xorigr = xmin
     yorigr = ymin
     dx = (xmax-xmin)/(nx-1)
     dy = (ymax-ymin)/(ny-1)
     !
     dX1 = dx     ! Scaled (dx,dy)
     dX2 = dy
     !
  END SELECT
  !
  !***  Reads zlayer from the database
  !
  if( mpime == root ) then
     call get_dbs_value_dimension(fdbs,'alt',zlayer,nz,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( zlayer, SIZE( zlayer ), root )
  !
  !***  Computes the array dz
  !
  do k=1,nz-1
     dz(k)= zlayer(k+1)-zlayer(k)
  end do
  dz(0 ) = dz(1)               ! Extend to outer frame
  dz(nz) = dz(nz-1)
  !
  !***  Reads the topography
  !
  if( mpime == root ) then
     call get_dbs_value_plane(fdbs,ibdbs,'TOPG',0.0_rp,nx,ny,topg,ivoid,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( topg, SIZE( topg ), root )
  !
  !***  Reads land use
  !
  if( mpime == root ) then
     call get_dbs_value_plane(fdbs,ibdbs,'LDU',0.0_rp,nx,ny,ldu,ivoid,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) then
     call wriwar('reamet : no value for LDU. Assuming no data LDU = -99')
     LDU = -99
  end if
  call bcast( ldu, SIZE( ldu ), root )
  !
  !***  Computes the scaling factors (see Toon et al.1988)
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     !      For spherical coordinates in terrain following
     !      Hm1 = sin (colat)
     !      Hm2 = 1            (not used)
     !      Vm  = 1            (not used)
     !
     do j=1,ny
        lat   = yorigr + (j-1)*dy
        colat = (90.0_rp-lat)*pi/180.0_rp   ! colatitude in rad
        Hm1(:,j) = sin(colat)
     end do
     !
  case('UTM')
     !
     !      For cartesian coordinates in terrain following
     !      Hm1 = 1
     !      Hm2 = 1            (not used)
     !      Vm  = 1            (not used)
     !
     Hm1(1:nx,1:ny) = 1.0_rp
     !
  END SELECT
  !
  return
end subroutine reagrd
!
!
!
subroutine get_time_lag
  !****************************************************************
  !*
  !*    Calculates the time_lag with the dbs (if any)
  !*
  !****************************************************************
  use KindType
  use Master
  implicit none
  !
  logical     :: found
  integer(ip) :: iyr ,imo ,idy ,ihr ,imi ,ise
  !
  !*** Calculates time lag by iteration (that is, the time in seconds between
  !*** the problem origin and the DBS origin).
  !*** Both origins are referred to 0000UTC, but may belong to different days
  !
  found = .false.
  time_lag = 0.0_rp
  do while(.not.found)
     call addtime(idbsyr,idbsmo,idbsdy,0,iyr,imo,idy,ihr,imi,ise,time_lag)
     !
     if( (iyr == ibyr).and.(imo == ibmo).and.(idy == ibdy) ) then
        found = .true.
     else
        time_lag = time_lag + 86400.0_rp
     end if
     if(INT(time_lag) > iedbs) call runend('get_time_lag: Time lag not found. Check interval consistency')
  end do
  !
  return
end subroutine get_time_lag
