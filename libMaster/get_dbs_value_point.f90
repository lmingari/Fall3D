subroutine get_dbs_value_point &
     (fname,timesec,word,x,y,z,value,endsec,istat,message)
  !**************************************************************************
  !*
  !*    Gets a variable value at a point from a netCDF database
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the database file
  !*    integer         timesec  Time (in seconds) at which data is extracted
  !*                             Time origin is date at 0000UTC
  !*    character*(*)   word     Variable to select. Possibilities are
  !*
  !*         Time-independent 2D var: XLON, XLAT, UTM_X, UTM_Y, UTM_ZONE, ZONE
  !*         Time-dependent   2D var: PBLH,  UST,  TOPG
  !*         Time-dependent   3D var:    U,    V,     W,     T,       TV, P
  !*
  !*    real            x,y,z    Point coordinates (z terrain following)
  !*
  !*    OUTPUT:
  !*    real            value    Value
  !*    integer         endsec   Time (in seconds) at until which data is valid
  !*                             Time origin is date at 0000UTC
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !************************************************************************
  use kindtype
  use netcdf
  implicit none
  !
  character(len=*)      :: fname,message,word
  integer(ip)           :: timesec,endsec,istat
  real(rp)              :: value,x,y,z
  !
  character(len=s_mess)   :: mymessage
  character(len=10)       :: coord_sys
  integer(ip)             :: ilen,nx,ny,ix,iy,i
  integer(ip)             :: ncID,lonID,latID
  real   (rp)             :: t,s,sm,sp,tm,tp,myshape(4)
  real   (rp),allocatable :: lon(:),lat(:),var(:,:)
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  coord_sys(:) = ' '
  istat = 0
  !
  !*** Gets database values
  !
  call get_dbs_property_int(fname,'NX',nx,istat,message)
  if(istat /= 0) return
  call get_dbs_property_int(fname,'NY',ny,istat,message)
  if(istat /= 0) return
  call get_dbs_property_cha(fname,'COORDINATES',coord_sys,istat,message)
  if(istat /= 0) return
  !
  !*** Read lat/lon and search if the point lays in the domain
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) goto 101
  !
  allocate(lon(nx))
  if(TRIM(coord_sys).eq.'LON-LAT') then
     if( nf90_inq_varid(ncID,'lon',lonID) /= 0) goto 102
  else if(TRIM(coord_sys).eq.'UTM') then
     if( nf90_inq_varid(ncID,'x',lonID) /= 0) goto 102   ! use lon as auxiliar
  else
     goto 106
  end if
  if( nf90_get_var(ncID,lonID,lon) /= 0) goto 103
  if((x.lt.lon(1)).or.(x.gt.lon(nx))) goto 104
  !
  allocate(lat(ny))
  if(TRIM(coord_sys).eq.'LON-LAT') then
     if( nf90_inq_varid(ncID,'lat',latID) /= 0) goto 102
  else if(TRIM(coord_sys).eq.'UTM') then
     if( nf90_inq_varid(ncID,'y',latID) /= 0) goto 102  ! use lat as auxiliar
  else
     goto 106
  end if
  if( nf90_get_var(ncID,latID,lat) /= 0) goto 103
  if((y.lt.lat(1)).or.(y.gt.lat(ny))) goto 104
  !
  if( nf90_close(ncID) /= 0) goto 105
  !
  !*** Finds the interpolation factors
  !
  ix = nx-1
  do i=1,nx-1                              ! index ix position (left point)
     if( (x.ge.lon(i)).and.(x.le.lon(i+1)) ) ix = i
  end do
  s = (x-lon(ix))/(lon(ix+1)-lon(ix))  ! parameter s in (0,1)
  s = 2.0_rp*s-1.0_rp                          ! parameter s in (-1,1)
  !
  iy = ny-1
  do i=1,ny-1                              ! index ix position (left point)
     if( (y.ge.lat(i)).and.(y.le.lat(i+1)) ) iy = i
  end do
  t = (y-lat(iy))/(lat(iy+1)-lat(iy))  ! parameter s in (0,1)
  t = 2.0_rp*t-1.0_rp                    ! parameter t in (-1,1)
  !
  deallocate(lon)
  deallocate(lat)
  !
  !**  Value in a plane z at the required time
  !
  allocate(var(nx,ny))
  call get_dbs_value_plane(fname,timesec,word,z,nx,ny,var,endsec,istat,message)
  if(istat /= 0) return
  !
  !*** Interpolate
  !
  sm = 0.5_rp*(1.0_rp-s)
  tm = 0.5_rp*(1.0_rp-t)
  sp = 0.5_rp*(1.0_rp+s)
  tp = 0.5_rp*(1.0_rp+t)
  myshape(1) = sm*tm
  myshape(2) = sp*tm
  myshape(3) = sp*tp
  myshape(4) = sm*tp
  !
  value = myshape(1)*var(ix  ,iy  ) + &
       myshape(2)*var(ix+1,iy  ) + &
       myshape(3)*var(ix+1,iy+1) + &
       myshape(4)*var(ix  ,iy+1)
  !
  deallocate(var)
  !
  !***  Successful end
  !
  return
  !
  !***  List of errors
  !
101 istat = -1
  mymessage = 'get_dbs_value_point : Error in nf90_open for file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
102 istat = -1
  mymessage = 'get_dbs_value_point : Error getting netCDF ID for variable ' &
       //TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
103 istat = -1
  mymessage = 'get_dbs_value_point : Error reading variable '//TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
104 istat = -1
  mymessage = 'get_dbs_value_point: point not in the database domain'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
105 istat = -1
  mymessage = 'get_dbs_value_point : Error in nf90_close '
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
106 istat = -1
  mymessage = 'get_dbs_value_point : coordinate system '// &
       TRIM(coord_sys)//' not found '
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_dbs_value_point
