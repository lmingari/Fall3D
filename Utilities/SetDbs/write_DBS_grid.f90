subroutine write_dbs_grid(fname,nx,ny,nz,nt,lon,lat,zlev,time,timesec, &
     xutm,yutm,zutm,coord_sys)
  !**************************************************************************
  !*
  !*   Writes grid data (including coordinate variables) in netCDF format
  !*
  !**************************************************************************
  use KindType, ONLY :  ip,rp
  use InpOut
  use DBS_nc
  use netcdf
  implicit none
  !
  character(len=*) :: fname,coord_sys
  character(len=3) :: zutm(nx,ny)
  integer(ip)      :: nx,ny,nz,nt,iy,im,id
  real   (rp)      :: lon (nx,ny),lat (nx,ny),zlev(nz),time(nt),timesec(nt)
  real   (rp)      :: xutm(nx,ny),yutm(nx,ny)
  !
  real(rp) :: val
  !
  !*** Open netCDF file (define mode) and get ncID
  !
  !     if( nf90_create(TRIM(fname),NF90_CLOBBER, ncID) /= 0 ) &
  !         call runend('write_dbs_grid : Error in nf90_create')
  ! MN
  if( nf90_create(TRIM(fname),cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncID) /= 0 ) &
       call runend('write_dbs_grid : Error in nf90_create')
  !
  !**  Define dimensions
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( nf90_def_dim(ncID, lon_DBS_name , nx, nx_DBS_ID ) /= 0 ) &
          call runend('write_dbs_grd : error in nf90_def_dim for lon')
     if( nf90_def_dim(ncID, lat_DBS_name , ny, ny_DBS_ID ) /= 0 ) &
          call runend('write_dbs_grd : error in nf90_def_dim for lat')
  case('UTM')
     if( nf90_def_dim(ncID, x_DBS_name, nx, nx_DBS_ID ) /= 0 ) &
          call runend('write_dbs_grd : error in nf90_def_dim for x')
     if( nf90_def_dim(ncID, y_DBS_name, ny, ny_DBS_ID ) /= 0 ) &
          call runend('write_dbs_grd : error in nf90_def_dim for y')
  END SELECT
  if( nf90_def_dim(ncID, alt_DBS_name , nz, nz_DBS_ID ) /= 0 ) &
       call runend('write_dbs_grd : error in nf90_def_dim for alt')
  if( nf90_def_dim(ncID, time_DBS_name, nt, nt_DBS_ID) /= 0 ) &
       call runend('write_dbs_grd : error in nf90_def_dim for time')
  !
  !*** Define coordinate variables
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( nf90_def_var(ncID, lon_DBS_name ,NF90_FLOAT, (/nx_DBS_ID/), lon_DBS_ID) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_def_variable for variable '//TRIM(lon_DBS_name))
     if( nf90_def_var(ncID, lat_DBS_name ,NF90_FLOAT, (/ny_DBS_ID/), lat_DBS_ID) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_def_variable for variable '//TRIM(lat_DBS_name))
  case('UTM')
     if( nf90_def_var(ncID, x_DBS_name ,NF90_FLOAT, (/nx_DBS_ID/), x_DBS_ID) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_def_variable for variable '//TRIM(x_DBS_name))
     if( nf90_def_var(ncID, y_DBS_name ,NF90_FLOAT, (/ny_DBS_ID/), y_DBS_ID) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_def_variable for variable '//TRIM(y_DBS_name))
  END SELECT
  if( nf90_def_var(ncID, alt_DBS_name ,NF90_FLOAT, (/nz_DBS_ID/), alt_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_variable for variable '//TRIM(alt_DBS_name))
  if( nf90_def_var(ncID, time_DBS_name,NF90_DOUBLE, (/nt_DBS_ID/), time_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_variable for variable '//TRIM(time_DBS_name))
  !
  !*** Define other variables
  !
  !
  !***    Coordinate variables and time-independent variables
  !
  if( nf90_def_var(ncID, xlon_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), xlon_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(xlon_DBS_name))
  if( nf90_def_var(ncID, xlat_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), xlat_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(xlat_DBS_name))
  if( nf90_def_var(ncID, xutm_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), xutm_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(xutm_DBS_name))
  if( nf90_def_var(ncID, yutm_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), yutm_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(yutm_DBS_name))
  if( nf90_def_var(ncID, topg_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), topg_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(topg_DBS_name))
  if( nf90_def_var(ncID, ldu_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), ldu_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(ldu_DBS_name))
  if( nf90_def_var(ncID, soil_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), soil_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(soil_DBS_name))
  if( nf90_def_var(ncID, vegfra_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID/), vegfra_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(vegfra_DBS_name))
  !
  !***   2D variables
  !
  if( nf90_def_var(ncID, pblh_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), pblh_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(pblh_DBS_name))
  if( nf90_def_var(ncID, ust_DBS_name , NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), ust_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(ust_DBS_name))
  if( nf90_def_var(ncID, smoi_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), smoi_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(smoi_DBS_name))
  if( nf90_def_var(ncID, rmol_DBS_name , NF90_FLOAT,(/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), rmol_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(rmol_DBS_name))
  if( nf90_def_var(ncID, znt_DBS_name , NF90_FLOAT,(/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), znt_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(znt_DBS_name))
  if( nf90_def_var(ncID, spd10_DBS_name , NF90_FLOAT,(/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), spd10_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(spd10_DBS_name))
  if( nf90_def_var(ncID, L_DBS_name   , NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), L_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(L_DBS_name))
  if( nf90_def_var(ncID, prat_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nt_DBS_ID/), prat_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(prat_DBS_name))
  !
  !***   3D variables
  !
  if( nf90_def_var(ncID, u_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), u_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(u_DBS_name))
  if( nf90_def_var(ncID, v_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), v_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(u_DBS_name))
  if( nf90_def_var(ncID, w_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), w_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(u_DBS_name))
  if( nf90_def_var(ncID, T_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), T_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(T_DBS_name))
  if( nf90_def_var(ncID, Tp_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), Tp_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(Tp_DBS_name))
  if( nf90_def_var(ncID, p_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), p_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(p_DBS_name))
  if( nf90_def_var(ncID, qv_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), qv_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(qv_DBS_name))
  if( nf90_def_var(ncID, ro_DBS_name, NF90_FLOAT, (/nx_DBS_ID,ny_DBS_ID,nz_DBS_ID,nt_DBS_ID/), ro_DBS_ID) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_def_var for variable '//TRIM(ro_DBS_name))
  !
  !*** Define attributes for coordinate variables
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     attr_desc  = 'longitude. East positive'
     attr_units = 'degrees_east'
     if( nf90_put_att(ncID, lon_DBS_ID, 'units', attr_units) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, lon_DBS_ID, 'description', attr_desc) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     !
     attr_desc  = 'latitude. North positive'
     attr_units = 'degrees_north'
     if( nf90_put_att(ncID, lat_DBS_ID, 'units', attr_units) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, lat_DBS_ID, 'description', attr_desc) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     !
  case('UTM')
     attr_desc  = 'West-East UTM distance'
     attr_units = 'm'
     if( nf90_put_att(ncID, x_DBS_ID, 'units', attr_units) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, x_DBS_ID, 'description', attr_desc) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     !
     attr_desc  = 'South-North UTM distance'
     attr_units = 'm'
     if( nf90_put_att(ncID, y_DBS_ID, 'units', attr_units) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, y_DBS_ID, 'description', attr_desc) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     !
  END SELECT
  !
  attr_desc  = 'height. Terrain following levels'
  attr_units = 'm'
  if( nf90_put_att(ncID, alt_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, alt_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'time after 0000UTC'
  attr_units = 's'
  if( nf90_put_att(ncID, time_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, time_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  !*** Define attributes for other variables
  !
  attr_desc  = 'longitude. East positive'
  attr_units = 'degrees_east'
  if( nf90_put_att(ncID, xlon_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, xlon_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'latitude. North positive'
  attr_units = 'degrees_north'
  if( nf90_put_att(ncID, xlat_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, xlat_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'x-UTM coordinate'
  attr_units = 'm'
  if( nf90_put_att(ncID, xutm_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, xutm_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'y-UTM coordinate'
  attr_units = 'm'
  if( nf90_put_att(ncID, yutm_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, yutm_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'UTM zone'
  attr_units = '-'
  !
  attr_desc  = 'Topography'
  attr_units = 'm'
  if( nf90_put_att(ncID, topg_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, topg_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Land use category. USGS 24-class'
  attr_units = '-'
  if( nf90_put_att(ncID, ldu_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, ldu_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Soil category. 19-class'
  attr_units = '-'
  if( nf90_put_att(ncID, soil_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, soil_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Vegetation Fraction'
  attr_units = '-'
  if( nf90_put_att(ncID, vegfra_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, vegfra_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Precipitation rate'
  attr_units = 'mm s-1'
  if( nf90_put_att(ncID, prat_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, prat_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Planetary boundary layer height'
  attr_units = 'm'
  if( nf90_put_att(ncID, pblh_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, pblh_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'u*'
  attr_units = 'm/s'
  if( nf90_put_att(ncID, ust_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, ust_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Soil moisture (first soil layer)'
  attr_units = 'm3/m3'
  if( nf90_put_att(ncID, smoi_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, smoi_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Inv. Monin-Obukhov lenght'
  attr_units = '1/m'
  if( nf90_put_att(ncID, rmol_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, rmol_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Roughness length'
  attr_units = 'm'
  if( nf90_put_att(ncID, znt_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, znt_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Wind speed at 10 M'
  attr_units = 'm/s'
  if( nf90_put_att(ncID, spd10_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, spd10_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Monin-Obukhov lenght'
  attr_units = 'm'
  if( nf90_put_att(ncID, L_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, L_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'u-velocity. Tangent to lon'
  attr_units = 'm/s'
  if( nf90_put_att(ncID, u_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, u_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'v-velocity. Tangent to lat'
  attr_units = 'm/s'
  if( nf90_put_att(ncID, v_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, v_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'vertical velocity'
  attr_units = 'm/s'
  if( nf90_put_att(ncID, w_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, w_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Temperature'
  attr_units = 'K'
  if( nf90_put_att(ncID, T_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, T_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Potential temperature'
  attr_units = 'K'
  if( nf90_put_att(ncID, Tp_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, Tp_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Pressure'
  attr_units = 'Pa'
  if( nf90_put_att(ncID, p_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, p_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Specific humidity'
  attr_units = 'kg/kg'
  if( nf90_put_att(ncID, qv_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, qv_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_desc  = 'Air density'
  attr_units = 'kg/m3'
  if( nf90_put_att(ncID, ro_DBS_ID, 'units', attr_units) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, ro_DBS_ID, 'description', attr_desc) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  !*** Put global attributes
  !
  attr_title = 'Meteo database for Fall3d derived from '//TRIM(TYPE_DATA)
  if( nf90_put_att(ncID, NF90_GLOBAL, 'TITLE', attr_title) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  attr_title = TRIM(coord_sys)
  if( nf90_put_att(ncID, NF90_GLOBAL, 'COORDINATES', attr_title) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  iy = INT(time(1)/1e10_rp)
  im = INT( (time(1)-1e10_rp*iy)/1e8_rp )
  id = INT( (time(1)-1e10_rp*iy-1e8_rp*im)/1e6_rp )
  if( nf90_put_att(ncID, NF90_GLOBAL, 'YEAR', iy ) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'MONTH', im ) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'DAY', id ) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  if( nf90_put_att(ncID, NF90_GLOBAL, 'START_TIME', INT(timesec(1))) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  if( nf90_put_att(ncID, NF90_GLOBAL, 'END_TIME', INT(timesec(nt))) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  if( nf90_put_att(ncID, NF90_GLOBAL, 'TIME_INCREMENT',INT((timesec(nt)-timesec(1))/(nt-1))) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     val =  lon(1,1)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'LONMIN', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     val =  lat(1,1)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'LATMIN', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     val =  lon(nx,1)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'LONMAX', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     val =  lat(1,ny)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'LATMAX', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     !
  case('UTM')
     val =  xutm(1,1)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'XMIN', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     val =  yutm(1,1)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'YMIN', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     val =  xutm(nx,1)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'XMAX', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     val =  yutm(1,ny)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'YMAX', val) /= 0 ) &
          call runend('write_dbs_grid : error in nf90_put_att')
     !
  END SELECT
  !
  val =  0.0
  if( nf90_put_att(ncID, NF90_GLOBAL, 'ZMIN', val) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  val =  zlev(nz)
  if( nf90_put_att(ncID, NF90_GLOBAL, 'ZMAX', val) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NX', nx) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NY', ny) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NZ', nz) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NT', nt) /= 0 ) &
       call runend('write_dbs_grid : error in nf90_put_att')
  !
  !*** Leave the define mode
  !
  if( nf90_enddef(ncID) /= 0 ) call runend('write_dbs_grid: error in nf90_enddef')
  !
  !*** Write coordinate variables
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     if( nf90_put_var(ncID, lon_DBS_ID, lon(1:nx,1   )) /= 0 )  &
          call runend('write_dbs_grid: error in nf90_put_var for varialbe lon')
     !
     allocate(work_DBS(ny))
     work_DBS(1:ny) = lat(1,1:ny)
     if( nf90_put_var(ncID, lat_DBS_ID, work_DBS(1:ny)) /= 0 )  &
          call runend('write_dbs_grid: error in nf90_put_var for varialbe lat')
     deallocate(work_DBS)
     !
  case('UTM')
     if( nf90_put_var(ncID, x_DBS_ID, xutm(1:nx,1   )) /= 0 )  &
          call runend('write_dbs_grid: error in nf90_put_var for varialbe xutm')
     !
     allocate(work_DBS(ny))
     work_DBS(1:ny) = yutm(1,1:ny)
     if( nf90_put_var(ncID, y_DBS_ID, work_DBS(1:ny)) /= 0 )  &
          call runend('write_dbs_grid: error in nf90_put_var for varialbe yutm')
     deallocate(work_DBS)
     !
  END SELECT
  !
  if( nf90_put_var(ncID, alt_DBS_ID, zlev(1:nz)) /= 0 )  &
       call runend('write_dbs_grid: error in nf90_put_var for varialbe alt')
  if( nf90_put_var(ncID, time_DBS_ID, timesec(1:nt)) /= 0 )  &
       call runend('write_dbs_grid: error in nf90_put_var for varialbe time')
  !
  !*** Write other variables (time-independent)
  !
  if( nf90_put_var(ncID, xlon_DBS_ID, lon  ) /= 0 )  &
       call runend('write_dbs_grid: error in nf90_put_var for varialbe xlon')
  if( nf90_put_var(ncID, xlat_DBS_ID, lat  ) /= 0 )  &
       call runend('write_dbs_grid: error in nf90_put_var for varialbe xlat')
  if( nf90_put_var(ncID, xutm_DBS_ID, xutm ) /= 0 )  &
       call runend('write_dbs_grid: error in nf90_put_var for varialbe xutm')
  if( nf90_put_var(ncID, yutm_DBS_ID, yutm ) /= 0 )  &
       call runend('write_dbs_grid: error in nf90_put_var for varialbe yutm')
  !
  !*** Close the file
  !
  if( nf90_close(ncID) /= 0) call runend('write_dbs_grid : Error in nf90_close')
  !
  return
end subroutine write_dbs_grid
