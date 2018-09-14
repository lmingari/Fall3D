subroutine readinp
  !************************************************************
  !*
  !*    Reads data from the input file
  !*
  !************************************************************
  use KindType
  use InpOut
  use Master
  use TimeFun
  use netcdf
  implicit none
  !
  character(len=s_mess) :: message,cvoid
  character(len=s_name) :: Name
  integer(ip) :: istat,ncID,nDim, nVar, nAttr,VarID,DimID
  integer(ip) :: i1yr,i1mo,i1dy,i1hr,i1mi,i1se,it,iVar
  integer(ip) :: varType,varDims,ivoid(10)
  !
  !***  Reads data input file
  !
  !
  !***  Plot characteristics
  !
  call get_input_rea(luinpname,'CROP_DOMAIN','LONMIN',lonmin,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea(luinpname,'CROP_DOMAIN','LONMAX',lonmax,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea(luinpname,'CROP_DOMAIN','LATMIN',latmin,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea(luinpname,'CROP_DOMAIN','LATMAX',latmax,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  !***  Topography (wirtten in m)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_TOPOGRAPHY',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_topog = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_topog = .true.
     call get_input_npar(luinpname,'MAP_TOPOGRAPHY','CONTOUR_LEVELS',nval_topog,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_TOPOGRAPHY','CONTOUR_LEVELS',cval_topog,nval_topog,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_topog = '=nf/1.0/'
  unit_topog = 'm'
  !
  !***  Total load (wirtten in kg/m2)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_TOTAL_LOAD',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_load0 = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_load0 = .true.
     call get_input_npar(luinpname,'MAP_TOTAL_LOAD','CONTOUR_LEVELS',nval_load0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_TOTAL_LOAD','CONTOUR_LEVELS',cval_load0,nval_load0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_load0 = '=nf/1.0/'
  unit_load0 = 'kg/m2'
  !
  !***  Wet deposition load (wirtten in kg/m2)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_WET_LOAD',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_loadw = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_loadw = .true.
     call get_input_npar(luinpname,'MAP_WET_LOAD','CONTOUR_LEVELS',nval_loadw,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_WET_LOAD','CONTOUR_LEVELS',cval_loadw,nval_loadw,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_loadw = '=nf/1.0/'
  unit_loadw = 'kg/m2'
  !
  !***  Class load (wirtten in kg/m2)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_CLASS_LOAD',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_loadc = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_loadc = .true.
     call get_input_npar(luinpname,'MAP_CLASS_LOAD','CONTOUR_LEVELS',nval_loadc,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_CLASS_LOAD','CONTOUR_LEVELS',cval_loadc,nval_loadc,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_loadc = '=nf/1.0/'
  unit_loadc = 'kg/m2'
  !
  !***  Class wet deposition (wirtten in kg/m2)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_CLASS_WET',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_wetc = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_wetc = .true.
     call get_input_npar(luinpname,'MAP_CLASS_WET','CONTOUR_LEVELS',nval_wetc,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_CLASS_WET','CONTOUR_LEVELS',cval_wetc,nval_wetc,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_wetc = '=nf/1.0/'
  unit_wetc = 'kg/m2'
  !
  !***  Deposit thickness (wirtten in cm)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_DEPOSIT_THICKNESS',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_thick = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_thick = .true.
     call get_input_npar(luinpname,'MAP_DEPOSIT_THICKNESS','CONTOUR_LEVELS',nval_thick,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_DEPOSIT_THICKNESS','CONTOUR_LEVELS',cval_thick,nval_thick,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_cha(luinpname,'MAP_DEPOSIT_THICKNESS','UNITS',cvoid,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(TRIM(cvoid).eq.'CM'.or.TRIM(cvoid).eq.'cm') then
        fact_thick = '=nf/1.0/'
        unit_thick = 'cm'
     else if(TRIM(cvoid).eq.'MM'.or.TRIM(cvoid).eq.'mm') then
        fact_thick = '=nf/10.0/'
        unit_thick = 'mm'
     else if(TRIM(cvoid).eq.'M'.or.TRIM(cvoid).eq.'m') then
        fact_thick = '=nf/0.01/'
        unit_thick = 'm'
     end if
  end if
  !
  !***  Ground concentration (wirtten in gr/m3)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_CONCE_GROUND',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_concg = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_concg = .true.
     call get_input_npar(luinpname,'MAP_CONCE_GROUND','CONTOUR_LEVELS',nval_concg,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_CONCE_GROUND','CONTOUR_LEVELS',cval_concg,nval_concg,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_concg = '=nf/1.0/'
  unit_concg = 'gr/m3'
  !
  !***  PMxx ground (wirtten in gr/m3)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_PMxx_GROUND',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_PMxxg = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_PMxxg = .true.
     call get_input_npar(luinpname,'MAP_PMxx_GROUND','CONTOUR_LEVELS',nval_PMxxg,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_PMxx_GROUND','CONTOUR_LEVELS',cval_PMxxg,nval_PMxxg,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_PMxxg = '=nf/1.0/'
  unit_PMxxg = 'gr/m3'
  !
  !***  Column mass (wirtten in gr/m2)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_COLUMN_MASS',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_cumul = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_cumul = .true.
     call get_input_npar(luinpname,'MAP_COLUMN_MASS','CONTOUR_LEVELS',nval_cumul,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_COLUMN_MASS','CONTOUR_LEVELS',cval_cumul,nval_cumul,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_cumul = '=nf/1.0/'
  unit_cumul = 'gr/m2'
  !
  !***  PMxx column mass (wirtten in gr/m2)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_COLUMN_PMxx',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_PMxxc = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_PMxxc = .true.
     call get_input_npar(luinpname,'MAP_COLUMN_PMxx','CONTOUR_LEVELS',nval_PMxxc,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_COLUMN_PMxx','CONTOUR_LEVELS',cval_PMxxc,nval_PMxxc,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_PMxxc = '=nf/1.0/'
  unit_PMxxc = 'gr/m2'
  !
  !***  FL (wirtten in gr/m3)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_FLIGHT_LEVEL',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_fl = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_fl = .true.
     call get_input_npar(luinpname,'MAP_FLIGHT_LEVEL','CONTOUR_LEVELS',nval_fl,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_FLIGHT_LEVEL','CONTOUR_LEVELS',cval_fl,nval_fl,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_fl = '=nf/1.0/'
  unit_fl = 'gr/m3'
  !
  !***  AOD (wirtten in -)
  !
  call get_input_cha(luinpname,'POSTPROCESS','MAP_AOD',cvoid,1,istat,message)
  if(istat.ne.0) then
     if(istat.gt.0) call wriwar(message)
     pp_aod05 = .false.
  else if(TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     pp_aod05 = .true.
     call get_input_npar(luinpname,'MAP_AOD','CONTOUR_LEVELS',nval_aod05,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea(luinpname,'MAP_AOD','CONTOUR_LEVELS',cval_aod05,nval_aod05,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
  end if
  fact_aod05 = '=nf/1.0/'
  unit_aod05 = '(-)'
  !
  !***  Reads time data from the netCDF file
  !
  if( nf90_open(TRIM(luncname),NF90_NOWRITE, ncID) /= 0 ) call runend('readinp : Error in nf90_open')
  !
  if( nf90_get_att(ncID, NF90_GLOBAL, nc_iyr_name, iyr) /= 0 ) &
       call runend('readinp : Error in f90_get_att for '//TRIM(nc_iyr_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, nc_imo_name, imo) /= 0 ) &
       call runend('readinp : Error in f90_get_att for '//TRIM(nc_imo_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, nc_idy_name, idy) /= 0 ) &
       call runend('readinp : Error in f90_get_att for '//TRIM(nc_idy_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, nc_irunb_name, irunb) /= 0 ) &
       call runend('readinp : Error in f90_get_att for '//TRIM(nc_irunb_name))
  if( nf90_get_att(ncID, NF90_GLOBAL, nc_irune_name, irune) /= 0 ) &
       call runend('readinp : Error in f90_get_att for '//TRIM(nc_irune_name))
  !
  ihr = irunb/3600
  imi = (irunb - (3600*ihr))/60
  !
  Name  = 'time'
  DimID = 4
  if( nf90_inquire_dimension(ncID,DimID,Name,nt) /= 0 ) &
       call runend('readinp : Error in nf90_inquire_dimension')
  allocate(times(nt))
  allocate(dates(nt))
  !
  Name = 'time'
  if( nf90_inq_varid(ncID,Name,VarID) /= 0) &
       call runend('readinp : Error getting VarID')
  if( nf90_get_var(ncID,VarID,times)  /= 0) &
       call runend('readinp : error reading variable times')
  times = times*3600.0_rp  ! convert to s
  !
  !*** Gets the current time instant date = hh:mmZddmmmyyyy adding times to 00:00UTC
  !
  do it = 1,nt
     call addtime (iyr,imo,idy,0,0,i1yr,i1mo,i1dy,i1hr,i1mi,i1se,times(it))
     call get_date(i1yr,i1mo,i1dy,i1hr,i1mi,dates(it))
  end do
  !
  !*** Inquires the number of classes (if any) and its name from the netcdf file
  !
  if( nf90_inquire(ncID, nDim, nVar, nAttr) /= 0 ) &  ! Global Inquire. Gets nDim nVar nAttr
       call runend('readinp : Error in nf90_inquire')
  !
  nclass = 0
  nclassw = 0
  do iVar = 1,nVar
     if( nf90_inquire_variable(ncID,iVar,Name,varType,varDims,ivoid) /= 0 ) &
          call runend('readinp : Error in nf90_inquire_variable')
     if(Name(1:5).eq.'LOAD_') then
        nclass = nclass + 1
        classname(nclass) = TRIM(name)
     end if
     if(Name(1:4).eq.'WET_') then
        nclassw = nclassw + 1
        classnamew(nclassw) = TRIM(name)
     end if
  end do
  !
  return
end subroutine readinp
