subroutine read_ETA_grid(fname,npoin,lonp,latp)
  !********************************************************************
  !*
  !*   Reads ETA grid in NetCDF format. It also creates all the
  !*   interpolation arrays used later by read_ETA_data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use TimeFun
  use MathFun
  use ETA_nc
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
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_ETA : Error in nf90_open')
  !
  !*** Global Inquire. Gets nDim nVar nAttr
  !
  if( nf90_inquire(ncID, nDim, nVar, nAttr) /= 0 ) call runend('read_ETA : Error in nf90_inquire')
  !
  !*** Allocates memeory for time-independent variables
  !
  call ETA_getmem01
  !
  !**  Inquires the names and lengths of the dimensions
  !
  do iDim = 1,nDim
     if( nf90_inquire_dimension(ncID,iDim,dimName(iDim),dimLen(iDim)) /= 0 ) &
          call runend('read_ETA : Error in nf90_inquire_dimension')
  end do
  !
  !*** Inquires the names and dimensions of variables
  !
  do iVar = 1,nVar
     if( nf90_inquire_variable(ncID,iVar,varName(iVar),varType(iVar),varDims(iVar),ivoid) /= 0 ) &
          call runend('read_ETA : Error in nf90_inquire_variable')
     varDimID(iVar,1:nDimMax) = ivoid(1:nDimMax)
  end do
  !
  !*** Reads global attributes
  !
  ibmi_ETA = 0
  ibse_ETA = 0
  if( nf90_get_att(ncID, NF90_GLOBAL, ibyr_ETA_name, ibyr_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(ibyr_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibmo_ETA_name, ibmo_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(ibmo_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibdy_ETA_name, ibdy_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(ibdy_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibhr_ETA_name, ibhr_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(ibhr_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, dt_ETA_name, dt_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(dt_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmin_ETA_name, lonmin_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(lonmin_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmax_ETA_name, lonmax_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(lonmax_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmin_ETA_name, latmin_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(latmin_ETA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmax_ETA_name, latmax_ETA) /= 0 ) &
       call runend('read_ETA : Error in f90_get_att for '//TRIM(latmax_ETA_name))
  !
  !*** Determines the ETA grid dimensions in 4D
  !
  do iDim = 1,nDim
     if(TRIM(dimName(iDim)) == TRIM(nt_ETA_name)) then
        nt_ETA = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nx_ETA_name)) then
        nx_ETA = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(ny_ETA_name)) then
        ny_ETA = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nz_ETA_name)) then
        nz_ETA = dimLen(iDim)
     end if
  end do
  !
  if (nt_ETA     == 0) call runend('read_ETA : Error in dimensions.'//TRIM(nt_ETA_name)//'=0')
  if (nx_ETA     == 0) call runend('read_ETA : Error in dimensions.'//TRIM(nx_ETA_name)//'=0')
  if (ny_ETA     == 0) call runend('read_ETA : Error in dimensions.'//TRIM(ny_ETA_name)//'=0')
  if (nz_ETA     == 0) call runend('read_ETA : Error in dimensions.'//TRIM(nz_ETA_name)//'=0')
  !
  !*** Determinates the nodal conectivities of the ETA grid. This is necessary
  !*** in order to interpolate data onto a GENERIC set of points (which do not
  !*** necessary form a structurated mesh!).
  !
  nelem_ETA = (nx_ETA - 1)*(ny_ETA - 1)
  !
  !*** Allocates memory for more time-independent variables
  !
  call ETA_getmem02
  !
  !*** Reads timesec_ETA(:) in sec after 0000UTC
  !
  found = .false.
  iVar = 0
  do while(.not.found)
     iVar = iVar + 1
     if( iVar == nVar+1 ) &
          call runend('read_ETA : variable '//TRIM(time_ETA_name)//' not found')
     if( TRIM(varName(iVar)) == TRIM(time_ETA_name) ) then
        if( nf90_get_var(ncID,iVar,timesec_ETA)  /= 0) &
             call runend('read_ETA : error reading variable'//TRIM(time_ETA_name))
        found = .true.
     end if
  end do
  !
  !***  Calculates time_ETA(:) in format YYYYMMDDHHMMSS
  !
  do it = 1,nt_ETA
     call addtime(ibyr_ETA,ibmo_ETA,ibdy_ETA,0,iyr,imo,idy,ihr,imi,ise,timesec_ETA(it))
     time_ETA(it) = 1d10*iyr + 1d8 *imo + 1d6 *idy + 1d4 *ihr + 1d2*imi + ise
  end do
  !
  !*** Computes lat and lon for ETA mesh.
  !
  dlon = (lonmax_ETA-lonmin_ETA)/(nx_ETA-1)
  dlat = (latmax_ETA-latmin_ETA)/(ny_ETA-1)
  do iy = 1,ny_ETA
     do ix = 1,nx_ETA
        lon_ETA(ix,iy) = lonmin_ETA + (ix-1)*dlon
        lat_ETA(ix,iy) = latmin_ETA + (iy-1)*dlat
     end do
  end do
  !
  !*** Reads p levels for ETA mesh
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(p_ETA_name) ) then
        if( nf90_get_var(ncID,iVar,p_ETA) /= 0 ) call runend('read_ETA : Error in nf90_get_var')
        p_ETA = p_ETA*101.3_rp  ! (converted from mb to Pa)
     end if
  end do
  !
  !*** Builds the nodal connectivities
  !
  call set_lnods(nx_ETA,ny_ETA,nelem_ETA,lnods_ETA)
  !
  !*** Builds elements to points connectivities and computes the interpolation
  !*** factors in 2D. It also verifies that the points lay within the ETA domain.
  !
  call set_lelpo(nelem_ETA,npoin_DAT,lon_ETA,lat_ETA,lonp,latp,lnods_ETA, &
       lelpo_ETA,s_po_ETA,t_po_ETA)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_ETA : Error in nf90_close')
  !
  return
end subroutine read_ETA_grid
