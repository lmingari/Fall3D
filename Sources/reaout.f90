  subroutine reaout
  !****************************************************************
  !*
  !*    Reads output strategy from the input file
  !*
  !****************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  implicit none
  !
  character(len=s_mess) :: message,cvoid
  integer  (ip)  :: istat,i
  real     (rp)  :: rvoid
  !
  !*** Output block
  !
  if( mpime == root ) then
     call get_input_rea(finp,'OUTPUT','POSTPROCESS_TIME_INTERVAL_(HOURS)',rvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( rvoid , 1, root )
  print_time = rvoid*3600.0_rp                ! convert to s
  !
  if( mpime == root ) then
     call get_input_cha (finp,'OUTPUT','POSTPROCESS_3D_VARIABLES',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( cvoid, LEN( cvoid ), root )
  if(TRIM(cvoid) == 'YES') then
     print_3d_var = .true.
  else
     print_3d_var = .false.
  end if
  !
  if( mpime == root ) then
     call get_input_cha (finp,'OUTPUT','POSTPROCESS_CLASSES',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( cvoid, LEN( cvoid ), root )
  if(TRIM(cvoid) == 'YES') then
     print_class = .true.
  else
     print_class = .false.
  end if
  !
  if( mpime == root ) then
     call get_input_cha (finp,'OUTPUT','TRACK_POINTS',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( cvoid, LEN( cvoid ), root )
  if(TRIM(cvoid) == 'YES') then
     track_points = .true.
  else
     track_points = .false.
  end if
  !
  !*** Tracking points
  !
  if(track_points) then
     !
     !***    Get the number of points
     !
     if( mpime == root ) then
        call get_points_npts(fpts,npts,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( npts, 1, root )
     !
     !***    Allocates memory for vectors related to npts (time dependent)
     !
     call getmem3
     !
     !***    Loads vectors
     !
     if( mpime == root ) then
        call get_points_coordinates(fpts,name_pts,xpts,ypts,npts,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     do i=1,npts
        call bcast(name_pts(i), LEN( name_pts(i) ), root )
     end do
     call bcast(xpts , SIZE(xpts), root )
     call bcast(ypts , SIZE(ypts), root )
     !
     !***    Averiguates which points lay within the computational domain
     !
     SELECT CASE(coord_sys)
     case('LON-LAT')
        do i = 1,npts
           if( (xpts(i) >= lonmin).and.(xpts(i) <= lonmax).and. &
                (ypts(i) >= latmin).and.(ypts(i) <= latmax)) use_pts(i) = .true.
        end do
     case('UTM')
        do i = 1,npts
           if( (xpts(i) >= xmin).and.(xpts(i) <= xmax).and. &
                (ypts(i) >= ymin).and.(ypts(i) <= ymax)) use_pts(i) = .true.
        end do
     END SELECT
     !
     !***    Get the interpolation factors
     !
     call get_pts_st
     !
     !***    Set file names
     !
     call get_pts_fname
     !
  end if
  !
  return
end subroutine reaout
!
!
!
subroutine get_pts_fname
  !*******************************************************************
  !*
  !*   Gets the tracking points file names
  !*
  !*******************************************************************
  use KindType
  use InpOut
  implicit none
  !
  logical     :: go_on
  integer(ip) :: ilen,i
  !
  !*** Remove extension
  !
  i = LEN(fpts) + 1
  go_on = .true.
  do while(go_on)
     i = i-1
     if(fpts(i:i) == '.' .or. i == 2) then
        go_on = .false.
        ilen = i-1
     end if
  end do
  !
  do i = 1,npts
     name_file_pts(i) = fpts(1:ilen)//'.tps.'//TRIM(name_pts(i))  ! no extension
  end do
  !
  return
end subroutine get_pts_fname
!
!
!
subroutine get_pts_st
  !*******************************************************************
  !*
  !*   Gets the mesh indexes and interpolation factors for the tracking
  !*   points
  !*
  !********************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  implicit none
  !
  logical     :: found
  integer(ip) :: i,ix,iy
  real   (rp) :: dxo,dyo,xo,yo,x,y
  !
  !***  Loop over tracked points
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     dxo = (lonmax-lonmin)/(nx-1)
     dyo = (latmax-latmin)/(ny-1)
     xo = lonmin
     yo = latmin
     !
  case('UTM')
     !
     dxo = (xmax-xmin)/(nx-1)
     dyo = (ymax-ymin)/(ny-1)
     xo = xmin
     yo = ymin
     !
  END SELECT
  !
  do i = 1,npts
     if(use_pts(i)) then
        !
        found = .false.
        do ix = 1,nx-1
           x = xo + (ix-1)*dxo
           if((xpts(i) >= x).and.(xpts(i) <= (x+dxo))) then
              ipts(i) = ix
              spts(i) = (xpts(i)-x)/dxo
              spts(i) = 2.0_rp*spts(i) - 1.0_rp         ! (-1,1)
              found = .true.
           end if
        end do
        if(.not.found) call runend('get_pts_st: value not found')
        !
        found = .false.
        do iy = 1,ny-1
           y = yo + (iy-1)*dyo
           if((ypts(i) >= y).and.(ypts(i) <= (y+dyo))) then
              jpts(i) = iy
              tpts(i) = (ypts(i)-y)/dyo
              tpts(i) = 2.0_rp*tpts(i) - 1.0_rp         ! (-1,1)
              found = .true.
           end if
        end do
        if(.not.found) call runend('get_pts_st: value not found')
        !
     end if
  end do
  !
  return
end subroutine get_pts_st
