subroutine read_WRF_grid(fname,npoin,lonp,latp)
  !********************************************************************
  !*
  !*   Reads WRF grid in NetCDF format. It also creates all the
  !*   interpolation arrays used later by read_WRF_data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use TimeFun
  use MathFun
  use WRF_nc
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: npoin
  real   (rp)      :: lonp(npoin),latp(npoin)
  !
  logical     :: found
  integer(ip) :: ix,iy
  integer(ip) :: iyr ,imo ,idy ,ihr ,imi ,ise
  integer(ip) :: iyr2,imo2,idy2,ihr2
  !
  !*** Initializations
  !
  nDimMax   = NF90_MAX_VAR_DIMS
  npoin_DAT = npoin
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_WRF : Error in nf90_open')
  !
  !*** Global Inquire. Gets nDim nVar nAttr
  !
  if( nf90_inquire(ncID, nDim, nVar, nAttr) /= 0 ) call runend('read_WRF : Error in nf90_inquire')
  !
  !*** Allocates memeory for time-independent variables
  !
  call WRF_getmem01
  !
  !**  Inquires the names and lengths of the dimensions
  !
  do iDim = 1,nDim
     if( nf90_inquire_dimension(ncID,iDim,dimName(iDim),dimLen(iDim)) /= 0 ) &
          call runend('read_WRF : Error in nf90_inquire_dimension')
  end do
  !
  !*** Inquires the names and dimensions of variables
  !
  do iVar = 1,nVar
     if( nf90_inquire_variable(ncID,iVar,varName(iVar),varType(iVar),varDims(iVar),ivoid) /= 0 ) &
          call runend('read_WRF : Error in nf90_inquire_variable')
     varDimID(iVar,1:nDimMax) = ivoid(1:nDimMax)
  end do
  !
  !*** Reads global attributes needed later
  !
  if( nf90_get_att(ncID, NF90_GLOBAL, WRF_map_proj_name, WRF_map_proj) /= 0 ) &
       call runend('read_WRF : Error in f90_get_att for '//TRIM(WRF_map_proj_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, WRF_cen_lon_name, WRF_cen_lon) /= 0 ) &
       call runend('read_WRF : Error in f90_get_att for '//TRIM(WRF_cen_lon_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, WRF_truelat1_name, WRF_truelat1) /= 0 ) &
       call runend('read_WRF : Error in f90_get_att for '//TRIM(WRF_truelat1_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, WRF_truelat2_name, WRF_truelat2) /= 0 ) &
       call runend('read_WRF : Error in f90_get_att for '//TRIM(WRF_truelat2_name))
  !
  !*** Determines the WFR grid dimensions in 4D
  !
  do iDim = 1,nDim
     if(TRIM(dimName(iDim)) == TRIM(nt_WRF_name)) then
        nt_WRF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nx_WRF_name)) then
        nx_WRF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nxstag_WRF_name)) then
        nxstag_WRF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(ny_WRF_name)) then
        ny_WRF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nystag_WRF_name)) then
        nystag_WRF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nz_WRF_name)) then
        nz_WRF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nzstag_WRF_name)) then
        nzstag_WRF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nsoil_WRF_name)) then
        nsoil_WRF = dimLen(iDim)
     end if
  end do
  !
  if (nt_WRF     == 0) call runend('read_WRF : Error in dimensions.'//TRIM(nt_WRF_name)//'=0')
  if (nx_WRF     == 0) call runend('read_WRF : Error in dimensions.'//TRIM(nx_WRF_name)//'=0')
  if (ny_WRF     == 0) call runend('read_WRF : Error in dimensions.'//TRIM(ny_WRF_name)//'=0')
  if (nz_WRF     == 0) call runend('read_WRF : Error in dimensions.'//TRIM(nz_WRF_name)//'=0')
  if (nxstag_WRF == 0) call runend('read_WRF : Error in dimensions.'//TRIM(nxstag_WRF_name)//'=0')
  if (nystag_WRF == 0) call runend('read_WRF : Error in dimensions.'//TRIM(nystag_WRF_name)//'=0')
  if (nzstag_WRF == 0) call runend('read_WRF : Error in dimensions.'//TRIM(nzstag_WRF_name)//'=0')
  if (nsoil_WRF  == 0) call runend('read_WRF : Error in dimensions.'//TRIM(nsoil_WRF_name)//'=0')
  !
  !*** Determinates the nodal conectivities of the WRF grid. This is necessary
  !*** in order to interpolate data onto a GENERIC set of points (which do not
  !*** necessary form a structurated mesh!). WFR results are transported to the
  !*** non staggered grid (box center)
  !
  nelem_WRF = (nx_WRF - 1)*(ny_WRF - 1)
  !
  !*** Allocates memory for more time-independent variables
  !
  call WRF_getmem02
  !
  !*** Reads WRF time instants in format YYYY-MM-DD_HH:MM:SS
  !
  found = .false.
  iVar = 0
  do while(.not.found)
     iVar = iVar + 1
     if( iVar == nVar+1 ) &
          call runend('read_WRF : variable '//TRIM(time_WRF_name)//' not found')
     if( TRIM(varName(iVar)) == TRIM(time_WRF_name) ) then
        if( nf90_get_var(ncID,iVar,time_WRF_string)  /= 0) &
             call runend('read_WRF : error reading variable'//TRIM(time_WRF_name))
        found = .true.
     end if
  end do
  !
  !*** WRF begin time data and checks
  !
  ibyr_WRF = stoi1(time_WRF_string(1)(1 :1 ))*1000 + &
       stoi1(time_WRF_string(1)(2 :2 ))*100  + &
       stoi1(time_WRF_string(1)(3 :3 ))*10   + &
       stoi1(time_WRF_string(1)(4 :4 ))
  ibmo_WRF = stoi1(time_WRF_string(1)(6 :6 ))*10   + &
       stoi1(time_WRF_string(1)(7 :7 ))
  ibdy_WRF = stoi1(time_WRF_string(1)(9 :9 ))*10   + &
       stoi1(time_WRF_string(1)(10:10))
  ibhr_WRF = stoi1(time_WRF_string(1)(12:12))*10   + &
       stoi1(time_WRF_string(1)(13:13))
  ibmi_WRF = stoi1(time_WRF_string(1)(15:15))*10   + &
       stoi1(time_WRF_string(1)(16:16))
  ibse_WRF = stoi1(time_WRF_string(1)(18:18))*10   + &
       stoi1(time_WRF_string(1)(19:19))
  !
  !*** Averiguates the WRF time step by guess (1h increment iteration) and
  !*** comparing with the second value
  !
  if(nt_WRF.gt.1) then
     iyr2 = stoi1(time_WRF_string(2)(1 :1 ))*1000 + &
          stoi1(time_WRF_string(2)(2 :2 ))*100  + &
          stoi1(time_WRF_string(2)(3 :3 ))*10   + &
          stoi1(time_WRF_string(2)(4 :4 ))
     imo2 = stoi1(time_WRF_string(2)(6 :6 ))*10   + &
          stoi1(time_WRF_string(2)(7 :7 ))
     idy2 = stoi1(time_WRF_string(2)(9 :9 ))*10   + &
          stoi1(time_WRF_string(2)(10:10))
     ihr2 = stoi1(time_WRF_string(2)(12:12))*10   + &
          stoi1(time_WRF_string(2)(13:13))
     !
     dt_WRF = 0.
     found = .false.
     do while(.not.found)
        dt_WRF = dt_WRF + 3600.
        call addtime(ibyr_WRF,ibmo_WRF,ibdy_WRF,ibhr_WRF, &
             iyr,imo,idy,ihr,imi,ise,dt_WRF)
        if( (iyr.eq.iyr2).and.(imo.eq.imo2).and.(idy.eq.idy2).and.(ihr.eq.ihr2) ) found = .true.
        if( dt_WRF.gt.(24.*3600.) ) call runend('read_WRF_grid : unable to find dt_WRF')
     end do
  else
     dt_WRF = 0.
  end if
  !
  !*** Calculates timesec_WRF(:) in sec after 0000UTC
  !
  timesec_WRF(1) = 3600.*ibhr_WRF + 60.*ibmi_WRF + ibse_WRF
  do it_WRF = 2,nt_WRF
     timesec_WRF(it_WRF) = timesec_WRF(1) + (it_WRF-1)*dt_WRF
  end do
  !
  !***  Calculates time_WRF(:)    in format YYYYMMDDHHMMSS
  !
  do it_WRF = 1,nt_WRF
     time_WRF(it_WRF) = stoi1(time_WRF_string(it_WRF)(1 :1 ))*(1d13) + &
          stoi1(time_WRF_string(it_WRF)(2 :2 ))*(1d12) + &
          stoi1(time_WRF_string(it_WRF)(3 :3 ))*(1d11) + &
          stoi1(time_WRF_string(it_WRF)(4 :4 ))*(1d10) + &
          stoi1(time_WRF_string(it_WRF)(6 :6 ))*(1d9 ) + &
          stoi1(time_WRF_string(it_WRF)(7 :7 ))*(1d8 ) + &
          stoi1(time_WRF_string(it_WRF)(9 :9 ))*(1d7 ) + &
          stoi1(time_WRF_string(it_WRF)(10:10))*(1d6 ) + &
          stoi1(time_WRF_string(it_WRF)(12:12))*(1d5 ) + &
          stoi1(time_WRF_string(it_WRF)(13:13))*(1d4 ) + &
          stoi1(time_WRF_string(it_WRF)(15:15))*(1d3 ) + &
          stoi1(time_WRF_string(it_WRF)(16:16))*(1d2 ) + &
          stoi1(time_WRF_string(it_WRF)(18:18))*(1d1 ) + &
          stoi1(time_WRF_string(it_WRF)(19:19))*(1d0 )
  end do
  !
  !*** Reads lat for WRF meshes at box center
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(lat_WRF_name) ) then
        if( nf90_get_var(ncID,iVar,lat_WRF) /= 0) call runend('read_WRF : Error in nf90_get_var')
     end if
  end do
  !
  !*** Reads lon for WRF meshes at box center
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(lon_WRF_name) ) then
        if( nf90_get_var(ncID,iVar,lon_WRF) /= 0 ) call runend('read_WRF : Error in nf90_get_var')
     end if
     !
     !*** Check if the International Date Line is intersected. Note that WRF gives lon in the
     !*** range (-180,180). In order to have growing values (required by the interpolation) WRF
     !*** longitudes are shifted if necessary
     !
     if( ((lon_WRF(1,1).gt.0.0).and.(lon_WRF(nx_WRF,1).lt.0.0)).or. &
         ((lon_WRF(1,1).gt.0.0).and.(lon_WRF(1,ny_WRF).lt.0.0)) ) then
        do iy = 1,ny_WRF
        do ix = 1,nx_WRF
           if(lon_WRF(ix,iy).gt.0.0) lon_WRF(ix,iy) = -180.0-lon_WRF(ix,iy)
        end do
        end do
     end if
  end do
  !
  !*** Builds the nodal connectivities
  !
  call set_lnods(nx_WRF,ny_WRF,nelem_WRF,lnods_WRF)
  !
  !*** Builds elements to points connectivities and computes the interpolation
  !*** factors in 2D. It also verifies that the points lay within the WRF domain.
  !
  call set_lelpo(nelem_WRF,npoin_DAT,lon_WRF,lat_WRF,lonp,latp,lnods_WRF, &
       lelpo_WRF,s_po_WRF,t_po_WRF)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_WRF : Error in nf90_close')
  !
  return
end subroutine read_WRF_grid
