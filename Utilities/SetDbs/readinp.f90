subroutine readinp
  !********************************************************
  !*
  !*    Gets properties from the input file
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  use TimeFun
  use Coordinates
  use wind_rotation
  implicit none
  !
  integer(ip) :: istat,iyr,imo,idy,ihr,imi,ise
  integer(ip) :: ix,iy
  real(rp)    :: dlon,dlat,dx,dy
  character(len=s_mess) :: message
  character(len=3) :: utmzone
  character(len=3) :: rotate_wind_answer
  !
  !***  Initializations
  !
  ibhr = 0
  !
  !***  Reads TIME_UTC block
  !
  call get_input_int &
       (luinpname,'TIME_UTC','YEAR',ibyr,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_int &
       (luinpname,'TIME_UTC','MONTH',ibmo,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_int &
       (luinpname,'TIME_UTC','DAY',ibdy,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea &
       (luinpname,'TIME_UTC','BEGIN_METEO_DATA_(HOURS_AFTER_00)',timeb,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  timeb = timeb*3600.                  ! in s
  !
  call get_input_rea &
       (luinpname,'TIME_UTC','END_METEO_DATA_(HOURS_AFTER_00)',timee,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  timee = timee*3600.                  ! in s
  if(timee.lt.timeb) call runend('Incorrect time interval')
  !
  call get_input_rea &
       (luinpname,'TIME_UTC','TIME_STEP_METEO_DATA_(MIN)',dt,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  dt = dt*60.                          ! in s
  nt = INT((timee-timeb)/dt) + 1
  !
  !***  Reads GRID block
  !
  call get_input_cha &
       (luinpname,'GRID','COORDINATES',coord_sys,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  if(TRIM(coord_sys).ne.'UTM'    .and. &
       TRIM(coord_sys).ne.'LON-LAT') call runend('Inocrrect system of coordinates')
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     call get_input_rea &
          (luinpname,'GRID','LONMIN',lonmin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','LATMIN',latmin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','LONMAX',lonmax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','LATMAX',latmax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
  case('UTM')
     !
     call get_input_rea &
          (luinpname,'GRID','XMIN',xmin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','YMIN',ymin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','XMAX',xmax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','YMAX',ymax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_cha &
          (luinpname,'GRID','UTMZONE',utmzone,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
  END SELECT
  !
  call get_input_int &
       (luinpname,'GRID','NX',nx,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_int &
       (luinpname,'GRID','NY',ny,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_npar &
       (luinpname,'GRID','ZLAYER_(M)',nz,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)

  !
  ! Read parameters for wind rotation (HIDDEN TO USERS)
  ! Does not generate error message if ROTATE_WIND is absent
  rotate_wind = .false.   ! Default
  call get_input_cha &
       (luinpname,'GRID','ROTATE_WIND',rotate_wind_answer,1,istat,message)
  if(istat /= 0) then   ! Here error means rotate_wind=.false.
     rotate_wind = .false.   ! Default
  elseif(TRIM(rotate_wind_answer).eq.'YES'.or.TRIM(rotate_wind_answer).eq.'yes') then
     rotate_wind = .true.
     call get_input_rea &
          (luinpname,'GRID','ROTATION_LONPOLE',rotation_lonpole,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','ROTATION_LATPOLE',rotation_latpole,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea &
          (luinpname,'GRID','ROTATION_ANGLE',rotation_angle,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)

  endif
  !
  !***  Allocates memory
  !
  call getmem
  !
  !***  Reads zlayer
  !
  call get_input_rea &
       (luinpname,'GRID','ZLAYER_(M)',zlayer,nz,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  !***  Computes timesec.  Seconds after 0000UTC
  !
  do it=1,nt
     timesec(it) = timeb + (it-1)*dt
  end do
  !
  !***  Computes time in the format YYYYMMDDHHMMSS
  !
  do it = 1,nt
     call addtime(ibyr,ibmo,ibdy,ibhr,iyr,imo,idy,ihr,imi,ise,timesec(it))
     time(it) = 1d10*iyr + 1d8 *imo + 1d6 *idy + 1d4 *ihr + 1d2*imi + ise
  end do
  !
  !***   Computes lonp, latp, xutm and yutm
  !
  SELECT CASE(coord_sys)
  CASE('LON-LAT')
     !
     dlon = (lonmax-lonmin)/(nx-1)
     dlat = (latmax-latmin)/(ny-1)
     do iy = 1,ny
        do ix = 1,nx
           lonp(ix,iy) = lonmin + (ix-1)*dlon
           latp(ix,iy) = latmin + (iy-1)*dlat
        end do
     end do
     !
     do iy = 1,ny
        do ix = 1,nx
           call ll2utm(latp(ix,iy),lonp(ix,iy),zutm(ix,iy),xutm(ix,iy),yutm(ix,iy),WGS_84_DATUM,istat)
           if(istat == -1) call runend('Error converting to UTM')
        end do
     end do
     !
  CASE('UTM')
     !
     dx = (xmax-xmin)/(nx-1)
     dy = (ymax-ymin)/(ny-1)
     do iy = 1,ny
        do ix = 1,nx
           xutm(ix,iy) = xmin + (ix-1)*dx
           yutm(ix,iy) = ymin + (iy-1)*dy
        end do
     end do
     !
     do iy = 1,ny
        do ix = 1,nx
           zutm(ix,iy) = utmzone
           call utm2ll(xutm(ix,iy),yutm(ix,iy),zutm(ix,iy),lonp(ix,iy),latp(ix,iy),WGS_84_DATUM,istat)
           if(istat == -1) call runend('Error converting to LON-LAT')
        end do
     end do
     !
  END SELECT
  !
  return
end subroutine readinp
