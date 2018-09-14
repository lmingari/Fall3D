subroutine read_NCEP2_grid(fname,npoin,lonp,latp)
  !********************************************************************
  !*
  !*   Reads NCEP2 grid in NetCDF format. It also creates all the
  !*   interpolation arrays used later by read_NCEP2_data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use TimeFun
  use MathFun
  use NCEP2_nc
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: npoin
  real   (rp)      :: lonp(npoin),latp(npoin)
  !
  logical     :: found
  integer(ip) :: iyr ,imo ,idy ,ihr ,imi ,ise, it, ix,iy
  real   (rp) :: dlat,dlon
  !
  !*** Initializations
  !
  nDimMax   = NF90_MAX_VAR_DIMS
  npoin_DAT = npoin
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_NCEP2 : Error in nf90_open')
  !
  !*** Global Inquire. Gets nDim nVar nAttr
  !
  if( nf90_inquire(ncID, nDim, nVar, nAttr) /= 0 ) call runend('read_NCEP2 : Error in nf90_inquire')
  !
  !*** Allocates memeory for time-independent variables
  !
  call NCEP2_getmem01
  !
  !**  Inquires the names and lengths of the dimensions
  !
  do iDim = 1,nDim
     if( nf90_inquire_dimension(ncID,iDim,dimName(iDim),dimLen(iDim)) /= 0 ) &
          call runend('read_NCEP2 : Error in nf90_inquire_dimension')
  end do
  !
  !*** Inquires the names and dimensions of variables
  !
  do iVar = 1,nVar
     if( nf90_inquire_variable(ncID,iVar,varName(iVar),varType(iVar),varDims(iVar),ivoid) /= 0 ) &
          call runend('read_NCEP2 : Error in nf90_inquire_variable')
     varDimID(iVar,1:nDimMax) = ivoid(1:nDimMax)
  end do
  !
  !*** Reads global attributes
  !
  ibmi_NCEP2 = 0
  ibse_NCEP2 = 0
  if( nf90_get_att(ncID, NF90_GLOBAL, ibyr_NCEP2_name, ibyr_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(ibyr_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibmo_NCEP2_name, ibmo_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(ibmo_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibdy_NCEP2_name, ibdy_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(ibdy_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibhr_NCEP2_name, ibhr_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(ibhr_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, dt_NCEP2_name, dt_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(dt_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmin_NCEP2_name, lonmin_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(lonmin_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmax_NCEP2_name, lonmax_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(lonmax_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmin_NCEP2_name, latmin_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(latmin_NCEP2_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmax_NCEP2_name, latmax_NCEP2) /= 0 ) &
       call runend('read_NCEP2 : Error in f90_get_att for '//TRIM(latmax_NCEP2_name))
  !
  !*** Determines the NCEP2 grid dimensions in 4D
  !
  do iDim = 1,nDim
     if(TRIM(dimName(iDim)) == TRIM(nt_NCEP2_name)) then
        nt_NCEP2 = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nx_NCEP2_name)) then
        nx_NCEP2 = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(ny_NCEP2_name)) then
        ny_NCEP2 = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nz_NCEP2_name)) then
        nz_NCEP2 = dimLen(iDim)
     end if
  end do
  !
  if (nt_NCEP2     == 0) call runend('read_NCEP2 : Error in dimensions.'//TRIM(nt_NCEP2_name)//'=0')
  if (nx_NCEP2     == 0) call runend('read_NCEP2 : Error in dimensions.'//TRIM(nx_NCEP2_name)//'=0')
  if (ny_NCEP2     == 0) call runend('read_NCEP2 : Error in dimensions.'//TRIM(ny_NCEP2_name)//'=0')
  if (nz_NCEP2     == 0) call runend('read_NCEP2 : Error in dimensions.'//TRIM(nz_NCEP2_name)//'=0')
  !
  !*** Determinates the nodal conectivities of the NCEP2 grid. This is necessary
  !*** in order to interpolate data onto a GENERIC set of points (which do not
  !*** necessary form a structurated mesh!).
  !
  nelem_NCEP2 = (nx_NCEP2 - 1)*(ny_NCEP2 - 1)
  !
  !*** Allocates memory for more time-independent variables
  !
  call NCEP2_getmem02
  !
  !*** Reads timesec_NCEP2(:) in sec after 0000UTC
  !
  found = .false.
  iVar = 0
  do while(.not.found)
     iVar = iVar + 1
     if( iVar == nVar+1 ) &
          call runend('read_NCEP2 : variable '//TRIM(time_NCEP2_name)//' not found')
     if( TRIM(varName(iVar)) == TRIM(time_NCEP2_name) ) then
        if( nf90_get_var(ncID,iVar,timesec_NCEP2)  /= 0) &
             call runend('read_NCEP2 : error reading variable'//TRIM(time_NCEP2_name))
        found = .true.
     end if
  end do
  !
  !***  Calculates time_NCEP2(:) in format YYYYMMDDHHMMSS
  !
  do it = 1,nt_NCEP2
     call addtime(ibyr_NCEP2,ibmo_NCEP2,ibdy_NCEP2,0,iyr,imo,idy,ihr,imi,ise,timesec_NCEP2(it))
     time_NCEP2(it) = 1d10*iyr + 1d8 *imo + 1d6 *idy + 1d4 *ihr + 1d2*imi + ise
  end do
  !
  !*** Computes lat and lon for NCEP2 mesh (GLOBAL!). Origin is at
  !*** (lonmin_NCEP2,latmin_NCEP2)=(-180,-90).
  !
  dlon = (lonmax_NCEP2-lonmin_NCEP2)/(nx_NCEP2-1)
  dlat = (latmax_NCEP2-latmin_NCEP2)/(ny_NCEP2-1)
  do iy = 1,ny_NCEP2
     do ix = 1,nx_NCEP2
        lon_NCEP2(ix,iy) = lonmin_NCEP2 + (ix-1)*dlon
        lat_NCEP2(ix,iy) = latmin_NCEP2 + (iy-1)*dlat
     end do
  end do
  !
  !*** Reads p levels for NCEP2 mesh
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(p_NCEP2_name) ) then
        if( nf90_get_var(ncID,iVar,p_NCEP2) /= 0 ) call runend('read_NCEP2 : Error in nf90_get_var')
        p_NCEP2 = p_NCEP2*101.3_rp  ! (converted from mb to Pa)
     end if
  end do
  !
  !*** Builds the nodal connectivities
  !
  call set_lnods(nx_NCEP2,ny_NCEP2,nelem_NCEP2,lnods_NCEP2)
  !
  !*** Builds elements to points connectivities and computes the interpolation
  !*** factors in 2D. It also verifies that the points lay within the NCEP2 domain.
  !
  call set_lelpo(nelem_NCEP2,npoin_DAT,lon_NCEP2,lat_NCEP2,lonp,latp,lnods_NCEP2, &
       lelpo_NCEP2,s_po_NCEP2,t_po_NCEP2)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_NCEP2 : Error in nf90_close')
  !
  return
end subroutine read_NCEP2_grid
