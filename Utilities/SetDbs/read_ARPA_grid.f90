subroutine read_ARPA_grid(fname,npoin,lonp,latp)
  !********************************************************************
  !*
  !*   Reads ARPA grid in NetCDF format. It also creates all the
  !*   interpolation arrays used later by read_ARPA_data
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use TimeFun
  use MathFun
  use ARPA_nc
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: npoin
  real   (rp)      :: lonp(npoin),latp(npoin)
  !
  logical     :: found
  integer(ip) :: iyr ,imo ,idy ,ihr ,imi ,ise, it
  integer :: status
  real(rp),allocatable :: tmpVar2d(:,:)
  !
  !*** Initializations
  !
  nDimMax   = NF90_MAX_VAR_DIMS
  npoin_DAT = npoin
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) call runend('read_ARPA : Error in nf90_open')
  !
  !*** Global Inquire. Gets nDim nVar nAttr
  !
  if( nf90_inquire(ncID, nDim, nVar, nAttr) /= 0 ) call runend('read_ARPA : Error in nf90_inquire')
  !
  !*** Allocates memeory for time-independent variables
  !
  call ARPA_getmem01
  !
  !**  Inquires the names and lengths of the dimensions
  !
  do iDim = 1,nDim
     if( nf90_inquire_dimension(ncID,iDim,dimName(iDim),dimLen(iDim)) /= 0 ) &
          call runend('read_ARPA : Error in nf90_inquire_dimension')
     ! write(*,*) iDim,trim(dimName(iDim)),dimLen(iDim)
  end do
  !
  !*** Inquires the names and dimensions of variables
  !
  do iVar = 1,nVar
     if( nf90_inquire_variable(ncID,iVar,varName(iVar),varType(iVar), &
          varDims(iVar),dimids=varDimID(iVar,:)) /= 0 ) &
          call runend('read_ARPA : Error in nf90_inquire_variable')
     !  varDimID(iVar,1:nDimMax) = ivoid(1:nDimMax)
     !  write(*,*) trim(varName(ivar)),varDimID(iVar,1:varDims(iVar))
  end do
  !
  !*** Reads global attributes
  !
  ibmi_ARPA = 0
  ibse_ARPA = 0
  if( nf90_get_att(ncID, NF90_GLOBAL, ibyr_ARPA_name, ibyr_ARPA) /= 0 ) &
       call runend('read_ARPA : Error in f90_get_att for '//TRIM(ibyr_ARPA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibmo_ARPA_name, ibmo_ARPA) /= 0 ) &
       call runend('read_ARPA : Error in f90_get_att for '//TRIM(ibmo_ARPA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibdy_ARPA_name, ibdy_ARPA) /= 0 ) &
       call runend('read_ARPA : Error in f90_get_att for '//TRIM(ibdy_ARPA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ibhr_ARPA_name, ibhr_ARPA) /= 0 ) &
       call runend('read_ARPA : Error in f90_get_att for '//TRIM(ibhr_ARPA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, dt_ARPA_name, dt_ARPA) /= 0 ) &
       call runend('read_ARPA : Error in f90_get_att for '//TRIM(dt_ARPA_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ARPA_cen_lon_name, ARPA_cen_lon) /= 0 ) &
       call runend('read_ARPA : Error in f90_get_att for '//TRIM(ARPA_cen_lon_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, ARPA_cen_lat_name, ARPA_cen_lat) /= 0 ) &
       call runend('read_ARPA : Error in f90_get_att for '//TRIM(ARPA_cen_lat_name))
  !
  !*** Determines the ARPA grid dimensions in 4D
  !
  do iDim = 1,nDim
     if(TRIM(dimName(iDim)) == TRIM(nt_ARPA_name)) then
        nt_ARPA_dimid = iDim
        nt_ARPA = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nx_ARPA_name)) then
        nx_ARPA_dimid = iDim
        nx_ARPA = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(ny_ARPA_name)) then
        ny_ARPA_dimid = iDim
        ny_ARPA = dimLen(iDim)
     else if(TRIM(dimName(iDim)) == TRIM(nz_ARPA_name)) then
        nz_ARPA_dimid = iDim
        nz_ARPA = dimLen(iDim)
     end if
  end do
  !
  if (nt_ARPA     == 0) call runend('read_ARPA : Error in dimensions.'//TRIM(nt_ARPA_name)//'=0')
  if (nx_ARPA     == 0) call runend('read_ARPA : Error in dimensions.'//TRIM(nx_ARPA_name)//'=0')
  if (ny_ARPA     == 0) call runend('read_ARPA : Error in dimensions.'//TRIM(ny_ARPA_name)//'=0')
  if (nz_ARPA     == 0) call runend('read_ARPA : Error in dimensions.'//TRIM(nz_ARPA_name)//'=0')
  !
  !*** Determinates the nodal conectivities of the ARPA grid. This is necessary
  !*** in order to interpolate data onto a GENERIC set of points (which do not
  !*** necessary form a structurated mesh!).
  !
  nelem_ARPA = (nx_ARPA - 1)*(ny_ARPA - 1)
  !
  !*** Allocates memory for more time-independent variables
  !

  call ARPA_getmem02
  !
  !*** Reads timesec_ARPA(:) in sec after 0000UTC
  !
  found = .false.
  iVar = 0
  do while(.not.found)
     iVar = iVar + 1
     if( iVar == nVar+1 ) &
          call runend('read_ARPA : variable '//TRIM(time_ARPA_name)//' not found')
     if( TRIM(varName(iVar)) == TRIM(time_ARPA_name) ) then
        if( nf90_get_var(ncID,iVar,timesec_ARPA)  /= 0) &
             call runend('read_ARPA : error reading variable'//TRIM(time_ARPA_name))
        found = .true.
     end if
  end do
  !
  !***  Calculates time_ARPA(:) in format YYYYMMDDHHMMSS
  !
  do it = 1,nt_ARPA
     call addtime(ibyr_ARPA,ibmo_ARPA,ibdy_ARPA,0,iyr,imo,idy,ihr,imi,ise,timesec_ARPA(it))
     time_ARPA(it) = 1d10*iyr + 1d8 *imo + 1d6 *idy + 1d4 *ihr + 1d2*imi + ise
  end do
  !
  !*** Reads lat for ARPA mesh
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(lat_ARPA_name) ) then
        if(varDims(iVar) /= 2) call runend('RLAT dimensions /= 2')
        if(varDimID(iVar,1)==nx_ARPA_dimid .and. varDimID(iVar,2)==ny_ARPA_dimid) then
           do_transpose= .false.
        else if(varDimID(iVar,2)==nx_ARPA_dimid.and.varDimID(iVar,1)==ny_ARPA_dimid) then
           do_transpose= .true.
        else
           call runend('RLAT order of dimensions is confusing')
        end if
        if(do_transpose) then
           allocate(tmpVar2d(ny_ARPA,nx_ARPA))
        else
           allocate(tmpVar2d(nx_ARPA,ny_ARPA))
        end if
        status=nf90_get_var(ncID,iVar,tmpVar2d)
        if( status /= 0) then
           call runend('read_ARPA lat mesh: Error in nf90_get_var: ' &
                //TRIM(varName(iVar))//' ('//TRIM(lat_ARPA_name)//') '// &
                nf90_strerror(status))
        end if
        if(do_transpose) then
           lat_ARPA = transpose(tmpVar2d)
        else
           lat_ARPA = tmpVar2d
        end if
        deallocate(tmpVar2d)
        exit  ! Found: exit from loop
     end if
  end do
  !
  !*** Reads lon for ARPA mesh
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(lon_ARPA_name) ) then
        if(varDims(iVar) /= 2) call runend('RLON dimensions /= 2')
        if(varDimID(iVar,1)==nx_ARPA_dimid .and. varDimID(iVar,2)==ny_ARPA_dimid) then
           if(do_transpose) call runend('Inconsistent value of do_transpose')
           do_transpose= .false.
        else if(varDimID(iVar,2)==nx_ARPA_dimid.and.varDimID(iVar,1)==ny_ARPA_dimid) then
           do_transpose= .true.
        else
           call runend('RLON order of dimensions is confusing')
        end if
        if(do_transpose) then
           allocate(tmpVar2d(ny_ARPA,nx_ARPA))
        else
           allocate(tmpVar2d(nx_ARPA,ny_ARPA))
        end if
        status=nf90_get_var(ncID,iVar,tmpVar2d)
        if(status /= 0 ) then
           call runend('read_ARPA lon mesh : Error in nf90_get_var: ' &
                //TRIM(varName(iVar))//' ('//TRIM(lon_ARPA_name)//') '// &
                nf90_strerror(status))
        end if
        if(do_transpose) then
           lon_ARPA = transpose(tmpVar2d)
        else
           lon_ARPA = tmpVar2d
        end if
        deallocate(tmpVar2d)
        exit  ! Found: exit from loop
     end if
  end do
  !
  !*** Reads p levels for ARPA mesh
  !
  do iVar = 1,nVar
     if( TRIM(varName(iVar)) == TRIM(p_ARPA_name) ) then
        if( nf90_get_var(ncID,iVar,p_ARPA) /= 0 ) &
             call runend('read_ARPA p-levels: Error in nf90_get_var: ' &
             //TRIM(varName(iVar))//' ('//TRIM(lat_ARPA_name)//')' )
        p_ARPA = p_ARPA*101.3_rp  ! (converted from mb to Pa)
        exit  ! Found: exit from loop
     end if
  end do
  !
  !*** Builds the nodal connectivities
  !
  call set_lnods(nx_ARPA,ny_ARPA,nelem_ARPA,lnods_ARPA)
  !
  !*** Builds elements to points connectivities and computes the interpolation
  !*** factors in 2D. It also verifies that the points lay within the ARPA domain.
  !
  call set_lelpo(nelem_ARPA,npoin_DAT,lon_ARPA,lat_ARPA,lonp,latp,lnods_ARPA, &
       lelpo_ARPA,s_po_ARPA,t_po_ARPA)
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) call runend('read_ARPA : Error in nf90_close')
  !
  return
end subroutine read_ARPA_grid
