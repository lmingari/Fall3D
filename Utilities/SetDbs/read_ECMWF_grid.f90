subroutine read_ECMWF_grid(fname,npoin,lonp,latp)
  !********************************************************************
  !*
  !*   Reads ECMWF grid in NetCDF format. It also creates all the
  !*   interpolation arrays used later by read_ECMWF_data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use TimeFun
  use MathFun
  use ECMWF_nc
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
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_ECMWF : Error in nf90_open')
  !
  !*** Global Inquire. Gets nDim nVar nAttr
  !
  if( nf90_inquire(ncID, nDim, nVar, nAttr) /= 0 ) call runend('read_ECMWF : Error in nf90_inquire')
  !
  !*** Allocates memeory for time-independent variables
  !
  call ECMWF_getmem01
  !
  !**  Inquires the names and lengths of the dimensions
  !
  do iDim = 1,nDim
     if( nf90_inquire_dimension(ncID,iDim,dimName(iDim),dimLen(iDim)) /= 0 ) &
          call runend('read_ECMWF : Error in nf90_inquire_dimension')
  end do
  !
  !*** Inquires the names and dimensions of variables
  !
  do iVar = 1,nVar
     if( nf90_inquire_variable(ncID,iVar,varName(iVar),varType(iVar),varDims(iVar),ivoid) /= 0 ) &
          call runend('read_ECMWF : Error in nf90_inquire_variable')
     varDimID(iVar,1:nDimMax) = ivoid(1:nDimMax)
  end do
  !
  !*** Reads global attributes
  !
  ibmi_ECMWF = 0
  ibse_ECMWF = 0
  if( nf90_get_att(ncID, NF90_GLOBAL, ibyr_ECMWF_name, ibyr_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(ibyr_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibmo_ECMWF_name, ibmo_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(ibmo_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibdy_ECMWF_name, ibdy_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(ibdy_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibhr_ECMWF_name, ibhr_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(ibhr_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, dt_ECMWF_name, dt_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(dt_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmin_ECMWF_name, lonmin_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(lonmin_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmax_ECMWF_name, lonmax_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(lonmax_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmin_ECMWF_name, latmin_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(latmin_ECMWF_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmax_ECMWF_name, latmax_ECMWF) /= 0 ) &
       call runend('read_ECMWF : Error in f90_get_att for '//TRIM(latmax_ECMWF_name))
  !
  !*** Determines the ECMWF grid dimensions in 4D
  !
  do iDim = 1,nDim
     if(TRIM(dimName(iDim)) == TRIM(nt_ECMWF_name)) then
        nt_ECMWF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nx_ECMWF_name)) then
        nx_ECMWF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(ny_ECMWF_name)) then
        ny_ECMWF = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nz_ECMWF_name)) then
        nz_ECMWF = dimLen(iDim)
     end if
  end do
  !
  if (nt_ECMWF     == 0) call runend('read_ECMWF : Error in dimensions.'//TRIM(nt_ECMWF_name)//'=0')
  if (nx_ECMWF     == 0) call runend('read_ECMWF : Error in dimensions.'//TRIM(nx_ECMWF_name)//'=0')
  if (ny_ECMWF     == 0) call runend('read_ECMWF : Error in dimensions.'//TRIM(ny_ECMWF_name)//'=0')
  if (nz_ECMWF     == 0) call runend('read_ECMWF : Error in dimensions.'//TRIM(nz_ECMWF_name)//'=0')
  !
  !*** Determinates the nodal conectivities of the ECMWF grid. This is necessary
  !*** in order to interpolate data onto a GENERIC set of points (which do not
  !*** necessary form a structurated mesh!).
  !
  nelem_ECMWF = (nx_ECMWF - 1)*(ny_ECMWF - 1)
  !
  !*** Allocates memory for more time-independent variables
  !
  call ECMWF_getmem02
  !
  !*** Reads timesec_ECMWF(:) in sec after 0000UTC
  !
  found = .false.
  iVar = 0
  do while(.not.found)
     iVar = iVar + 1
     if( iVar == nVar+1 ) &
          call runend('read_ECMWF : variable '//TRIM(time_ECMWF_name)//' not found')
     if( TRIM(varName(iVar)) == TRIM(time_ECMWF_name) ) then
        if( nf90_get_var(ncID,iVar,timesec_ECMWF)  /= 0) &
             call runend('read_ECMWF : error reading variable'//TRIM(time_ECMWF_name))
        found = .true.
     end if
  end do
  !
  !***  Calculates time_ECMWF(:) in format YYYYMMDDHHMMSS
  !
  do it = 1,nt_ECMWF
     call addtime(ibyr_ECMWF,ibmo_ECMWF,ibdy_ECMWF,0,iyr,imo,idy,ihr,imi,ise,timesec_ECMWF(it))
     time_ECMWF(it) = 1d10*iyr + 1d8 *imo + 1d6 *idy + 1d4 *ihr + 1d2*imi + ise
  end do
  !
  !*** Computes lat and lon for ECMWF mesh.
  !
  dlon = (lonmax_ECMWF-lonmin_ECMWF)/(nx_ECMWF-1)
  dlat = (latmax_ECMWF-latmin_ECMWF)/(ny_ECMWF-1)
  do iy = 1,ny_ECMWF
     do ix = 1,nx_ECMWF
        lon_ECMWF(ix,iy) = lonmin_ECMWF + (ix-1)*dlon
        lat_ECMWF(ix,iy) = latmin_ECMWF + (iy-1)*dlat
     end do
  end do
  !
  !*** Reads p levels for ECMWF mesh
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(p_ECMWF_name) ) then
        if( nf90_get_var(ncID,iVar,p_ECMWF) /= 0 ) call runend('read_ECMWF : Error in nf90_get_var')
        p_ECMWF = p_ECMWF*101.3_rp  ! (converted from mb to Pa)
     end if
  end do
  !
  !*** Builds the nodal connectivities
  !
  call set_lnods(nx_ECMWF,ny_ECMWF,nelem_ECMWF,lnods_ECMWF)
  !
  !*** Builds elements to points connectivities and computes the interpolation
  !*** factors in 2D. It also verifies that the points lay within the ECMWF domain.
  !
  call set_lelpo(nelem_ECMWF,npoin_DAT,lon_ECMWF,lat_ECMWF,lonp,latp,lnods_ECMWF, &
       lelpo_ECMWF,s_po_ECMWF,t_po_ECMWF)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_ECMWF : Error in nf90_close')
  !
  return
end subroutine read_ECMWF_grid
