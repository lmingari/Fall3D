  subroutine rearst
  !**********************************************************
  !*
  !*    Reads the NetCDF restart file
  !*
  !**********************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use netcdf
  use Parallel
  implicit none
  !
  character(len=s_name) :: dimName
  integer(ip) :: istat
  integer(ip) :: iDim,nDim,nVar,nAttr,dimLen
  integer(ip) :: nx_rst,ny_rst,nz_rst,nc_rst
  integer(ip) :: ibyr_rst,ibmo_rst,ibdy_rst
  integer(ip) :: ibeg_time_rst,iend_time_rst
  real   (rp) :: lonmin_rst,lonmax_rst,latmin_rst,latmax_rst
  real   (rp), allocatable :: work3(:,:,:),work4(:,:,:,:)
  !
  !*** Open netCDF file and get ncID
  !
  if( mpime == root ) then
     istat = nf90_open(TRIM(frst),NF90_NOWRITE, ncID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in nf90_open')
  !
  !*** Global Inquire. Gets nDim nVar nAttr
  !
  if( mpime == root ) then
     istat = nf90_inquire(ncID, nDim, nVar, nAttr)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in nf90_inquire')
  call bcast( nDim , 1, root )
  call bcast( nVar , 1, root )
  call bcast( nAttr, 1, root )
  !
  !**  Inquires the names and lengths of the dimensions
  !
  do iDim = 1,nDim
     if( mpime == root ) then
        istat = nf90_inquire_dimension(ncID,iDim,dimName,dimLen)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('rearst : Error in nf90_inquire_dimension')
     call bcast( dimName, LEN(dimName), root )
     call bcast  ( dimLen , 1           , root )
     !
     if(TRIM(dimName) == 'nx' ) nx_rst = dimLen
     if(TRIM(dimName) == 'ny' ) ny_rst = dimLen
     if(TRIM(dimName) == 'nz2') nz_rst = dimLen - 2
     if(TRIM(dimName) == 'nc' ) nc_rst = dimLen
  end do
  !
  !*** Reads global attributes
  !
  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'YEAR', ibyr_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( ibyr_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'MONTH', ibmo_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( ibmo_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'DAY', ibdy_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( ibdy_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'RUN_START', ibeg_time_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( ibeg_time_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'RUN_END'  , iend_time_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( iend_time_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'ERUPTED_MASS', erumass)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( erumass, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'OUT_MASS', outmass)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( outmass, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'LONMIN', lonmin_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( lonmin_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'LONMAX', lonmax_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( lonmax_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'LATMIN', latmin_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error in f90_get_att ')
  call bcast( latmin_rst, 1, root )

  if( mpime == root ) then
     istat = nf90_get_att(ncID, NF90_GLOBAL, 'LATMAX', latmax_rst)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('readres : Error in f90_get_att')
  call bcast( latmax_rst, 1, root )
  !
  !*** Checks for consistency and sets beg_time
  !
  if(nx_rst /= nx) call runend('restart inconsistency in nx')
  if(ny_rst /= ny) call runend('restart inconsistency in ny')
  if(nz_rst /= nz) call runend('restart inconsistency in nz')
  if(nc_rst /= nc) call runend('restart inconsistency in nc')
  !
  if(lonmin_rst /= lonmin) call runend('restart inconsistency in lonmin')
  if(lonmax_rst /= lonmax) call runend('restart inconsistency in lonmax')
  if(latmin_rst /= latmin) call runend('restart inconsistency in latmin')
  if(latmax_rst /= latmax) call runend('restart inconsistency in latmax')
  !
  !*** Allocates memory
  !
  allocate(work4(nx+2,ny+2,nz+2,nc  ))
  allocate(work3(nx  ,ny  ,     nc+1))
  !
  !*** Reads values
  !
  if( mpime == root ) then
     istat = nf90_inq_varid(ncID,cload_nc_name,cload_nc_ID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error getting cload_nc_ID')

  if( mpime == root ) then
     istat = nf90_get_var(ncID,cload_nc_ID,work3)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error reading cload')
  call bcast( work3, SIZE(work3), root )
  !
  if( mpime == root ) then
     istat = nf90_inq_varid(ncID,conce_nc_name,conce_nc_ID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error getting conce_nc_ID')

  if( mpime == root ) then
     istat = nf90_get_var  (ncID,conce_nc_ID,work4)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('rearst : Error reading conce')
  call bcast( work4, SIZE(work4), root )
  !
  !*** Assigns and releases memory
  !
  cload(1:nx  ,1:ny         ,0:nc) = work3(:,:,:)
  c    (0:nx+1,0:ny+1,0:nz+1,1:nc) = work4(:,:,:,:)
  !
  deallocate(work3)
  deallocate(work4)
  !
  return
  end subroutine rearst
