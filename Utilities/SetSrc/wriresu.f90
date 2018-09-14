      subroutine wriresu_grid
!**************************************************************************
!*
!*   Writes source data (including coordinate variables) in netCDF format
!*   This is used for the resuspension case only
!*
!**************************************************************************
     use KindType, ONLY :  ip,rp
     use InpOut
     use Master
     use Deposit
     use netcdf
     implicit none
!
     integer(ip)      :: ic
     character(len=2) :: ext
!
!*** Open netCDF file (define mode) and get ncID
!
     if( nf90_create(TRIM(lusncname),cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncID) /= 0 ) &
         call runend(' wriresu_grid : Error in nf90_create')
!
!**  Define dimensions
!
      if( nf90_def_dim(ncID, lon_RES_name   , nx    , nx_RES_ID ) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_dim for lon')
      if( nf90_def_dim(ncID, lat_RES_name   , ny    , ny_RES_ID ) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_dim for lat')
      if( nf90_def_dim(ncID, nc_RES_name    , nc    , nc_RES_ID ) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_dim for nc')
      if( nf90_def_dim(ncID, nc_dep_RES_name, nc_dep, nc_dep_RES_ID ) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_dim for nc_dep')
      if( nf90_def_dim(ncID, time_RES_name  , ndt   , nt_RES_ID ) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_dim for time')
!
!*** Define coordinate variables
!
      if( nf90_def_var(ncID, lon_RES_name ,NF90_FLOAT, (/nx_RES_ID/), lon_RES_ID) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_variable for variable '//TRIM(lon_RES_name))
      if( nf90_def_var(ncID, lat_RES_name ,NF90_FLOAT, (/ny_RES_ID/), lat_RES_ID) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_variable for variable '//TRIM(lat_RES_name))
      if( nf90_def_var(ncID, time_RES_name ,NF90_INT, (/nt_RES_ID/), time_RES_ID) /= 0 ) &
          call runend('wriresu_grid : error in nf90_def_variable for variable '//TRIM(time_RES_name))
!
!*** Define other variables
!
     if( nf90_def_var(ncID, load_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID/), load_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(load_RES_name))
!
     do ic = 1,nc_dep
         if( nf90_def_var(ncID, clas_RES_name(ic), NF90_FLOAT, (/nx_RES_ID,ny_RES_ID/), clas_RES_ID(ic)) /= 0 ) &
             call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(clas_RES_name(ic)))
     end do
!
     if( nf90_def_var(ncID, soil_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID/), soil_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(soil_RES_name))
!
     if( nf90_def_var(ncID, vegfra_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID/), vegfra_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(vegfra_RES_name))
!
     if( nf90_def_var(ncID, ust_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID,nt_RES_ID/), ust_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(ust_RES_name))
!
     if( nf90_def_var(ncID, smoi_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID,nt_RES_ID/), smoi_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(smoi_RES_name))
!     
     if( nf90_def_var(ncID, rmol_RES_name, NF90_FLOAT,(/nx_RES_ID,ny_RES_ID,nt_RES_ID/), rmol_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(rmol_RES_name))
!
     if( nf90_def_var(ncID, densi_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID,nt_RES_ID/), densi_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(densi_RES_name))
!
     if( nf90_def_var(ncID, tust_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID,nc_dep_RES_ID,nt_RES_ID/), tust_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(tust_RES_name))
!
     if( nf90_def_var(ncID, emis_RES_name, NF90_FLOAT, (/nx_RES_ID,ny_RES_ID,nc_RES_ID,nt_RES_ID/), emis_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(tust_RES_name))
!
     if( nf90_def_var(ncID, diam_RES_name, NF90_FLOAT, (/nc_RES_ID/), diam_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(diam_RES_name))
!
     if( nf90_def_var(ncID, rhop_RES_name, NF90_FLOAT, (/nc_RES_ID/), rhop_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(rhop_RES_name))
!
     if( nf90_def_var(ncID, sphe_RES_name, NF90_FLOAT, (/nc_RES_ID/), sphe_RES_ID) /= 0 ) &
         call runend('wriresu_grid : error in nf90_def_var for variable '//TRIM(sphe_RES_name))
!
!*** Define attributes for coordinate variables
!
      attr_desc  = 'longitude. East positive'
      attr_units = 'degrees_east'
      if( nf90_put_att(ncID, lon_RES_ID, 'units', attr_units) /= 0 ) &
         call runend('wriresu_grid : error in nf90_put_att')
      if( nf90_put_att(ncID, lon_RES_ID, 'description', attr_desc) /= 0 ) &
         call runend('wriresu_grid : error in nf90_put_att')
!
      attr_desc  = 'latitude. North positive'
      attr_units = 'degrees_north'
      if( nf90_put_att(ncID, lat_RES_ID, 'units', attr_units) /= 0 ) &
         call runend('wriresu_grid : error in nf90_put_att')
      if( nf90_put_att(ncID, lat_RES_ID, 'description', attr_desc) /= 0 ) &
         call runend('wriresu_grid : error in nf90_put_att')
!
      attr_desc  = 'time after 0000UTC'
      attr_units = 's'
      if( nf90_put_att(ncID, time_RES_ID, 'units', attr_units) /= 0 ) &
         call runend('wriresu_grid : error in nf90_put_att')
      if( nf90_put_att(ncID, time_RES_ID, 'description', attr_desc) /= 0 ) &
         call runend('wriresu_grid : error in nf90_put_att')
!
!*** Define attributes for other variables
!
     attr_desc  = 'Particle ground deposit load'
     attr_units = 'kg/m2'
     if( nf90_put_att(ncID, load_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, load_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     do ic = 1,nc_dep
        ext  = '00'
        if(ic.lt.10) then
          write(ext(2:2),'(i1.1)') ic
        else if(ic.le.100) then
          write(ext(1:2),'(i2.2)') ic
        end if
        attr_desc  = 'Deposit load for particle class '//TRIM(ext)
        attr_units = 'kg/m2'
        if( nf90_put_att(ncID, clas_RES_ID(ic), 'units', attr_units) /= 0 ) &
            call runend('wriresu_grid : error in nf90_put_att')
        if( nf90_put_att(ncID, clas_RES_ID(ic), 'description', attr_desc) /= 0 ) &
            call runend('wriresu_grid : error in nf90_put_att')
     end do
!
     attr_desc  = 'Soil categories'
     attr_units = '-'
     if( nf90_put_att(ncID, soil_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, soil_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Vegetation fraction'
     attr_units = '-'
     if( nf90_put_att(ncID, vegfra_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, vegfra_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Friction velocity'
     attr_units = 'm/s'
     if( nf90_put_att(ncID, ust_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, ust_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Soil moisture (first soil layer)'
     attr_units = 'm3/m3'
     if( nf90_put_att(ncID, smoi_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, smoi_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Inv. Monin-Obukhov length'
     attr_units = '1/m'
     if( nf90_put_att(ncID, rmol_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, rmol_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Air density at surface'
     attr_units = 'kg/m3'
     if( nf90_put_att(ncID, densi_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, densi_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Threshold friction velocity'
     attr_units = 'm/s'
     if( nf90_put_att(ncID, tust_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, tust_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Emission rate'
     attr_units = 'kg/s'
     if( nf90_put_att(ncID, emis_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, emis_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Particle diameter'
     attr_units = 'microns'
     if( nf90_put_att(ncID, diam_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, diam_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Particle density'
     attr_units = 'kg/m3'
     if( nf90_put_att(ncID, rhop_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, rhop_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_desc  = 'Particle sphericity'
     attr_units = '-'
     if( nf90_put_att(ncID, sphe_RES_ID, 'units', attr_units) /= 0 ) &
       call runend('wriresu_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, sphe_RES_ID, 'description', attr_desc) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
!*** Put global attributes
!
     attr_title = 'Source database for resuspension in Fall3d'
     if( nf90_put_att(ncID, NF90_GLOBAL, 'TITLE', attr_title) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     attr_title = TRIM(emission_scheme)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'EMISSION_SCHEME', attr_title) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if(moisture_correction) then
        attr_title = ' YES'
     else
        attr_title = ' NO'
     end if
     if( nf90_put_att(ncID, NF90_GLOBAL, 'MOISTURE_CORRECTION', attr_title) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'MAX_RESUSPENSION_SIZE_MM', diam_max) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'EMISSION_FACTOR', emission_factor) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if(TRIM(emission_scheme).eq.'WESTPHAL') then
        if( nf90_put_att(ncID, NF90_GLOBAL, 'THRESHOLD_UST',tust_cte) /= 0 ) &
            call runend('wriresu_grid : error in nf90_put_att')
     end if
!
     attr_title = TRIM(coord_sys)
     if( nf90_put_att(ncID, NF90_GLOBAL, 'COORDINATES', attr_title) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'YEAR', idbsyr) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'MONTH', idbsmo) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'DAY', idbsday) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'START_TIME', idbsbeg) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'END_TIME', idbsend) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
     if( nf90_put_att(ncID, NF90_GLOBAL, 'TIME_INCREMENT', idbsincr) /= 0 ) &
        call runend('wriresu_grid : error in nf90_put_att')
!
!*** Leave the define mode
!
     if( nf90_enddef(ncID) /= 0 ) call runend('wriresu_grid : error in nf90_enddef')
!
!*** Write coordinate variables
!
     if( nf90_put_var(ncID, lon_RES_ID, lon(1:nx)) /= 0 )  &
         call runend('wriresu_grid : error in nf90_put_var for varialbe lon')
!
     if( nf90_put_var(ncID, lat_RES_ID, lat(1:ny)) /= 0 )  &
          call runend('wriresu_grid : error in nf90_put_var for varialbe lat')
!
     
     if( nf90_put_var(ncID, time_RES_ID, idt_src(1:ndt)) /= 0 )  &
          call runend('wriresu_grid : error in nf90_put_var for varialbe time')
!
!*** Write other variables (time-independent)
!
     if( nf90_put_var(ncID, load_RES_ID, load ) /= 0 )  &
        call runend('wriresu_grid : error in nf90_put_var for varialbe load')
!
     do ic = 1,nc_dep
        if( nf90_put_var(ncID, clas_RES_ID(ic), clas(1:nx,1:ny,ic)  ) /= 0 )  &
          call runend('wriresu_grid : error in nf90_put_var for varialbe clas')
     end do
!
     if( nf90_put_var(ncID, diam_RES_ID, 1d6 * diam ) /= 0 )  &
        call runend('wriresu_grid : error in nf90_put_var for varialbe diameter')
!
     if( nf90_put_var(ncID, rhop_RES_ID, rhop ) /= 0 )  &
        call runend('wriresu_grid : error in nf90_put_var for varialbe density')
!
     if( nf90_put_var(ncID, sphe_RES_ID, sphe ) /= 0 )  &
        call runend('wriresu_grid : error in nf90_put_var for varialbe sphericity')
!
!*** Close the file
!
     if( nf90_close(ncID) /= 0) call runend('wriresu_grid : Error in nf90_close')
!
     return
     end subroutine wriresu_grid
!
!
!
    subroutine wriresu_res
!**************************************************************************
!*
!*   Writes source data results for the given time step in netCDF format
!*   Must be called after wriresu_grid
!*   This is used for the resuspension case only
!*
!**************************************************************************
    use KindType, ONLY :  ip,rp
    use InpOut
    use Master
    use Deposit
    use netcdf
    implicit none
    !
    !*** Open netCDF file and get ncID
    !
    if( nf90_open(TRIM(lusncname),NF90_WRITE, ncID) /= 0 ) call runend('wriresu_res : Error in nf90_open')
    !
    !*** Put variables for the current time instant (idt)
    !
    if( nf90_inq_varid(ncID,ust_RES_name,ust_RES_ID) /= 0) &
        call runend('wriresu_res : Error getting ust_DBS_ID')
    if( nf90_put_var(ncID, ust_RES_ID, ust, start=(/1,1,idt/), count=(/nx,ny,1/) ) /= 0 ) &
        call wriwar('wriresu_res : error in nf90_put_var for variable ust')
    !
    if( nf90_inq_varid(ncID,smoi_RES_name,smoi_RES_ID) /= 0) &
        call runend('wriresu_res : Error getting smoi_DBS_ID')
    if( nf90_put_var(ncID, smoi_RES_ID, smoi, start=(/1,1,idt/), count=(/nx,ny,1/) ) /= 0 ) &
        call wriwar('wriresu_res : error in nf90_put_var for variable smoi')
    !
    if( nf90_inq_varid(ncID,rmol_RES_name,rmol_RES_ID) /= 0) &
        call runend('wriresu_res : Error getting rmol_DBS_ID')
    if( nf90_put_var(ncID, rmol_RES_ID, rmol, start=(/1,1,idt/),count=(/nx,ny,1/) ) /= 0 ) &
        call wriwar('wriresu_res : error in nf90_put_var for variable rmol')
    !
    if( nf90_inq_varid(ncID,densi_RES_name,densi_RES_ID) /= 0) &
        call runend('wriresu_res : Error getting ust_DBS_ID')
    if( nf90_put_var(ncID, densi_RES_ID, densi, start=(/1,1,idt/), count=(/nx,ny,1/) ) /= 0 ) &
        call wriwar('wriresu_res : error in nf90_put_var for variable densi')
    !
    if( nf90_inq_varid(ncID,tust_RES_name,tust_RES_ID) /= 0) &
        call runend('wriresu_res : Error getting tust_DBS_ID')
    if( nf90_put_var(ncID, tust_RES_ID, tust, start=(/1,1,1,idt/), count=(/nx,ny,nc_dep,1/) ) /= 0 ) &
        call wriwar('wriresu_res : error in nf90_put_var for variable tust')
    !
    if( nf90_inq_varid(ncID,emis_RES_name,emis_RES_ID) /= 0) &
        call runend('wriresu_res : Error getting emis_DBS_ID')
    if( nf90_put_var(ncID, emis_RES_ID, emis, start=(/1,1,1,idt/), count=(/nx,ny,nc,1/) ) /= 0 ) &
        call wriwar('wriresu_res : error in nf90_put_var for variable emis')
    !
    !*** Put soil categories (no time dependent)
    !
    if( nf90_put_var(ncID, soil_RES_ID, soil) /= 0 ) &
        call wriwar('wriresu_res : error in nf90_put_var for variable soil')
    !
    if( nf90_inq_varid(ncID,vegfra_RES_name,vegfra_RES_ID) /= 0) &
        call runend('wriresu_res : Error getting vegfra_DBS_ID')
    if( nf90_put_var(ncID, vegfra_RES_ID, vegfra) /= 0 ) &
        call runend('wriresu_res : error in nf90_put_var for variable vegfra')
    !
    !*** Close the file
    !
    if( nf90_close(ncID) /= 0) call runend('wriresu_res : Error in nf90_close')
    !
    return
    end subroutine wriresu_res
