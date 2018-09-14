subroutine read_GFS_grid(fname,npoin,lonp,latp)
  !********************************************************************
  !*
  !*   Reads GFS grid in NetCDF format. It also creates all the
  !*   interpolation arrays used later by read_GFS_data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use TimeFun
  use MathFun
  use GFS_nc
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
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_GFS : Error in nf90_open')
  !
  !*** Global Inquire. Gets nDim nVar nAttr
  !
  if( nf90_inquire(ncID, nDim, nVar, nAttr) /= 0 ) call runend('read_GFS : Error in nf90_inquire')
  !
  !*** Allocates memeory for time-independent variables
  !
  call GFS_getmem01
  !
  !**  Inquires the names and lengths of the dimensions
  !
  do iDim = 1,nDim
     if( nf90_inquire_dimension(ncID,iDim,dimName(iDim),dimLen(iDim)) /= 0 ) &
          call runend('read_GFS : Error in nf90_inquire_dimension')
  end do
  !
  !*** Inquires the names and dimensions of variables
  !
  do iVar = 1,nVar
     if( nf90_inquire_variable(ncID,iVar,varName(iVar),varType(iVar),varDims(iVar),ivoid) /= 0 ) &
          call runend('read_GFS : Error in nf90_inquire_variable')
     varDimID(iVar,1:nDimMax) = ivoid(1:nDimMax)
  end do
  !
  !*** Reads global attributes
  !
  ibmi_GFS = 0
  ibse_GFS = 0
  if( nf90_get_att(ncID, NF90_GLOBAL, ibyr_GFS_name, ibyr_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(ibyr_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibmo_GFS_name, ibmo_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(ibmo_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibdy_GFS_name, ibdy_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(ibdy_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibhr_GFS_name, ibhr_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(ibhr_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, dt_GFS_name, dt_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(dt_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmin_GFS_name, lonmin_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(lonmin_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, lonmax_GFS_name, lonmax_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(lonmax_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmin_GFS_name, latmin_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(latmin_GFS_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, latmax_GFS_name, latmax_GFS) /= 0 ) &
       call runend('read_GFS : Error in f90_get_att for '//TRIM(latmax_GFS_name))
  !
  !*** Determines the GFS grid dimensions in 4D
  !
  do iDim = 1,nDim
     if(TRIM(dimName(iDim)) == TRIM(nt_GFS_name)) then
        nt_GFS = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nx_GFS_name)) then
        nx_GFS = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(ny_GFS_name)) then
        ny_GFS = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nz_GFS_name)) then
        nz_GFS = dimLen(iDim)
     end if
  end do
  !
  if (nt_GFS     == 0) call runend('read_GFS : Error in dimensions.'//TRIM(nt_GFS_name)//'=0')
  if (nx_GFS     == 0) call runend('read_GFS : Error in dimensions.'//TRIM(nx_GFS_name)//'=0')
  if (ny_GFS     == 0) call runend('read_GFS : Error in dimensions.'//TRIM(ny_GFS_name)//'=0')
  if (nz_GFS     == 0) call runend('read_GFS : Error in dimensions.'//TRIM(nz_GFS_name)//'=0')
  !
  !*** Determinates the nodal conectivities of the GFS grid. This is necessary
  !*** in order to interpolate data onto a GENERIC set of points (which do not
  !*** necessary form a structurated mesh!).
  !
  nelem_GFS = (nx_GFS - 1)*(ny_GFS - 1)
  !
  !*** Allocates memory for more time-independent variables
  !
  call GFS_getmem02
  !
  !*** Reads timesec_GFS(:) in sec after 0000UTC
  !
  found = .false.
  iVar = 0
  do while(.not.found)
     iVar = iVar + 1
     if( iVar == nVar+1 ) &
          call runend('read_GFS : variable '//TRIM(time_GFS_name)//' not found')
     if( TRIM(varName(iVar)) == TRIM(time_GFS_name) ) then
        if( nf90_get_var(ncID,iVar,timesec_GFS)  /= 0) &
             call runend('read_GFS : error reading variable'//TRIM(time_GFS_name))
        found = .true.
     end if
  end do
  !
  !***  Calculates time_GFS(:) in format YYYYMMDDHHMMSS
  !
  do it = 1,nt_GFS
     call addtime(ibyr_GFS,ibmo_GFS,ibdy_GFS,0,iyr,imo,idy,ihr,imi,ise,timesec_GFS(it))
     time_GFS(it) = 1d10*iyr + 1d8 *imo + 1d6 *idy + 1d4 *ihr + 1d2*imi + ise
  end do
  !
  !*** Computes lat and lon for GFS mesh (GLOBAL!). Origin is at
  !*** (lonmin_GFS,latmin_GFS)=(-180,-90).
  !
  dlon = (lonmax_GFS-lonmin_GFS)/(nx_GFS-1)
  dlat = (latmax_GFS-latmin_GFS)/(ny_GFS-1)
  do iy = 1,ny_GFS
     do ix = 1,nx_GFS
        lon_GFS(ix,iy) = lonmin_GFS + (ix-1)*dlon
        lat_GFS(ix,iy) = latmin_GFS + (iy-1)*dlat
     end do
  end do
  !
  !*** Reads p levels for GFS mesh
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(p_GFS_name) ) then
        if( nf90_get_var(ncID,iVar,p_GFS) /= 0 ) call runend('read_GFS : Error in nf90_get_var')
        p_GFS = p_GFS*101.3_rp  ! (converted from mb to Pa)
     end if
  end do
  !
  !*** Builds the nodal connectivities
  !
  call set_lnods(nx_GFS,ny_GFS,nelem_GFS,lnods_GFS)
  !
  !*** Builds elements to points connectivities and computes the interpolation
  !*** factors in 2D. It also verifies that the points lay within the GFS domain.
  !
  call set_lelpo(nelem_GFS,npoin_DAT,lon_GFS,lat_GFS,lonp,latp,lnods_GFS, &
       lelpo_GFS,s_po_GFS,t_po_GFS)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_GFS : Error in nf90_close')
  !
  return
end subroutine read_GFS_grid
