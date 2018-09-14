  subroutine wrirst
  !**************************************************************************
  !*
  !*   Writes the restart file in netCDF format
  !*
  !**************************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use netcdf
  use Parallel
  implicit none
  !
  integer(ip) :: istat
  integer(ip) :: nx_ID ,ny_ID        ,nc_ID
  integer(ip) :: nx2_ID,ny2_ID,nz2_ID,nc1_ID
  !
  !*** Open netCDF file (define mode) and get ncID
  !
  if( mpime == root ) then
     istat = nf90_create(TRIM(frst),NF90_CLOBBER, ncID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : Error in nf90_create')
  !
  !**  Define dimensions
  !
  if( mpime == root ) then
     istat = nf90_def_dim(ncID, 'nx', nx, nx_ID )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 )  call runend('wrirst : error in nf90_def_dim for nx')

  if( mpime == root ) then
     istat = nf90_def_dim(ncID, 'ny', ny, ny_ID )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_dim for ny')

  if( mpime == root ) then
     istat = nf90_def_dim(ncID, 'nc', nc, nc_ID )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_dim for nc')

  if( mpime == root ) then
     istat = nf90_def_dim(ncID, 'nc1', nc+1, nc1_ID )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_dim for nc1')

  if( mpime == root ) then
     istat = nf90_def_dim(ncID, 'nx2', nx+2, nx2_ID )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_dim for nx2')

  if( mpime == root ) then
     istat = nf90_def_dim(ncID, 'ny2', ny+2, ny2_ID )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_dim for ny2')

  if( mpime == root ) then
     istat = nf90_def_dim(ncID, 'nz2', nz+2, nz2_ID )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_dim for nz2')
  !
  !*** Define variables
  !
  if( mpime == root ) then
     istat = nf90_def_var(ncID, conce_nc_name, NF90_FLOAT,(/nx2_ID,ny2_ID,nz2_ID,nc_ID/),conce_nc_ID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_var for variable '//TRIM(conce_nc_name))

  if( mpime == root ) then
     istat = nf90_def_var(ncID, cload_nc_name, NF90_FLOAT,(/nx_ID,ny_ID,nc1_ID/)        ,cload_nc_ID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_def_var for variable '//TRIM(cload_nc_name))
  !
  !*** Define attributes
  !
  attr_desc  = 'Scaled ground deposit load'
  attr_units = 'kg/m2'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, cload_nc_ID, 'units', attr_units)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')

  if( mpime == root ) then
     istat = nf90_put_att(ncID, cload_nc_ID, 'description', attr_desc)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  attr_desc  = 'Scaled concentration'
  attr_units = 'gr/m3'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, conce_nc_ID, 'units', attr_units)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')

  if( mpime == root ) then
     istat = nf90_put_att(ncID, conce_nc_ID, 'description', attr_desc)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  !*** Put global attributes
  !
  attr_title = 'Fall3d 6.2 results'
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'TITLE', attr_title)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  attr_title = TRIM(coord_sys)
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'COORDINATES', attr_title)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'YEAR', ibyr )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'MONTH', ibmo )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'DAY', ibdy )
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'RUN_START', INT(beg_time))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'RUN_END', INT(end_time))
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'ERUPTED_MASS', erumass)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  if( mpime == root ) then
     istat = nf90_put_att(ncID, NF90_GLOBAL, 'OUT_MASS', outmass)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LONMIN', lonmin)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LATMIN', latmin)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LONMAX', lonmax)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'LATMAX', latmax)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
  case('UTM')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'XMIN', xmin)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'YMIN', ymin)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'XMAX', xmax)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
     if( mpime == root ) then
        istat = nf90_put_att(ncID, NF90_GLOBAL, 'YMAX', ymax)
     end if
     call bcast( istat, 1, root )
     if( istat /= 0 ) call runend('wrirst : error in nf90_put_att')
     !
  END SELECT
  !
  !*** Leave the define mode
  !
  if( mpime == root ) then
     istat = nf90_enddef(ncID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst: error in nf90_enddef')

  if( mpime == root ) then
     istat = nf90_close (ncID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst: error in nf90_close')
  !
  !***  Start to write results
  !***  Open netCDF file and get ncID
  !
  if( mpime == root ) then
     istat = nf90_open(TRIM(frst),NF90_WRITE, ncID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : Error in nf90_open')
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID,conce_nc_ID,c)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_var for variable c')
  !
  if( mpime == root ) then
     istat = nf90_put_var(ncID,cload_nc_ID,cload)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('wrirst : error in nf90_put_var for variable cload')
  !
  !*** Close the file
  !
  if( mpime == root ) then
     istat = nf90_close(ncID)
  end if
  call bcast( istat, 1, root )
  if( istat /= 0 ) call runend('printres_nc : Error in nf90_close')
  !
  return
end subroutine wrirst
