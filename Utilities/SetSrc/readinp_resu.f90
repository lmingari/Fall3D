  subroutine readinp_resu
  !********************************************************
  !*
  !*    Gets properties from the input file for the
  !*    resuspension mode
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  use Deposit
  use netcdf
  implicit none
  !
  logical               :: foundx,foundy,conve
  integer(ip)           :: istat,ic,ix,iy,iz,i,j,ix_dep,iy_dep
  real(rp)              :: s,t,shape(4),colat
  character(len=2)      :: ext
  character(len=s_work) :: cvoid
  character(len=s_mess) :: message
  real(rp)              :: rhoair0 = 1.2255, muair0=1.827d-5
  !
  !***  Initialization
  !
  lusncname = TRIM(lusrcname)//'.nc'
  write(lulog,1) TRIM(lusncname)
1 format(2x,'Resuspension   file   : ',a,/)
  !
  !*** Initialize names in deposit file
  !
  do ic = 1,nc_dep
     ext  = '00'
     if(ic.lt.10) then
        write(ext(2:2),'(i1.1)') ic
     else if(ic.le.100) then
        write(ext(1:2),'(i2.2)') ic
     end if
     clas_RES_name(ic) = TRIM(clas_RES_name(ic))//TRIM(ext)
  end do
  !
  !*** Read variables from the input file
  !
  call get_input_cha(luinpname,'RESUSPENSION_SOURCE','MOISTURE_CORRECTION',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  call upcase(cvoid)
  if(TRIM(cvoid).eq.'YES') then
     moisture_correction = .true.
     w_threshold = 3.0_rp           ! assumed threshold value of 0% for ash
  else
     moisture_correction = .false.
  end if
  !
  call get_input_cha(luinpname,'RESUSPENSION_SOURCE','VEGETATION_CORRECTION',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  call upcase(cvoid)
  if(TRIM(cvoid).eq.'YES') then
     vegetation_correction = .true.
  else
     vegetation_correction = .false.
  end if
  !
  call get_input_cha(luinpname,'RESUSPENSION_SOURCE','RECALCULATE_UST',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  call upcase(cvoid)
  if(TRIM(cvoid).eq.'YES') then
    recalculate_ust = .true.
    roughness_ust   = 0.01_rp
  else
    recalculate_ust = .false.
  endif
  !
  call get_input_cha(luinpname,'RESUSPENSION_SOURCE','EMISSION_SCHEME',emission_scheme,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  call upcase(emission_scheme)
  !
  call get_input_rea(luinpname,'RESUSPENSION_SOURCE','EMISSION_FACTOR',emission_factor,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  if(TRIM(emission_scheme).eq.'WESTPHAL') then
        call get_input_rea (luinpname,'RESUSPENSION_SOURCE','THRESHOLD_UST',tust_cte,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
  else if(TRIM(emission_scheme).eq.'MARTICORENA'.or.TRIM(emission_scheme).eq.'SHAO') then
        tust_cte = 0.0d0
  else
       call runend('Incorrect emission scheme type')
  end if
  !
  call get_input_rea(luinpname,'RESUSPENSION_SOURCE','MAX_INJECTION_HEIGHT_(M)',z_emission,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea(luinpname,'RESUSPENSION_SOURCE','DEPOSIT_THRESHOLD_(KGM2)',mindep,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  if(mindep.lt.mindep0) mindep=mindep0
  !
  !*** Reads MODEL block
  !
  call get_input_cha(luinpname,TRIM(model_name),'TERMINAL_VELOCITY_MODEL',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  if(TRIM(cvoid).eq.'ARASTOOPOUR') then
     modv = 1
  else if(TRIM(cvoid).eq.'GANSER') then
     modv = 2
  else if(TRIM(cvoid).eq.'WILSON') then
     modv = 3
  else if(TRIM(cvoid).eq.'DELLINO') then
     modv = 4
  else
     modv = 2
     call wriwar('Invalid settling velocity model. Ganser model assumed by default')
  end if
  !
  !*** Calculates psi depending on the velocity model
  !
  call setpsi(psi_dep,sphe_dep,diam_dep,modv,nc_dep)   ! Calculates psi(nc) depending on the model
  !
  do ic = 1,nc_dep
     call vsettl(diam_dep(ic),rhop_dep(ic),rhoair0,muair0,vset_dep(ic),modv,psi_dep(ic),conve)
     if(.not.conve) call runend('No convergence in vsettl routine')
  end do
  !
  !*** Reads the grid from the deposit file (netCDF). Note that the grids of both files
  !*** (deposit and dbs used in the simulation of resuspension) can be different
  !
  call get_dbs_property_cha(ludepname,'COORDINATES',coord_sys_dep,istat,message)
  if(TRIM(coord_sys_dep).ne.'LON-LAT') call runend('readinp_resu : UTM coordinates not implemented')
  !
  call get_dbs_dimension(ludepname,'lon',nx_dep,istat,message)
  call get_dbs_dimension(ludepname,'lat',ny_dep,istat,message)
  !
  allocate(lon_dep(nx_dep))
  allocate(lat_dep(ny_dep))
  !
  call get_dbs_value_dimension(ludepname,'lon',lon_dep,nx_dep,istat,message)
  call get_dbs_value_dimension(ludepname,'lat',lat_dep,ny_dep,istat,message)
  !
  xorigr_dep = lon_dep(1)
  yorigr_dep = lat_dep(1)
  dx_dep     = (lon_dep(nx_dep)-lon_dep(1))/(nx_dep-1)
  dy_dep     = (lat_dep(ny_dep)-lat_dep(1))/(ny_dep-1)
  !
  !*** Reads the deposit load (total and for each class)
  !
  allocate(work_dep(nx_dep,ny_dep))
  allocate(load_dep(nx_dep,ny_dep))
  allocate(clas_dep(nx_dep,ny_dep,nc_dep))
  !
  call get_dbs_dimension(ludepname,'time',nt_dep,istat,message)
  !
  if( nf90_open(TRIM(ludepname),NF90_NOWRITE, ncID) /= 0 ) &
      call runend('readinp_resu : Error in nf90_open for deposit file ')
  !
  !                                  Total load
  !
  if( nf90_inq_varid(ncID,load_RES_name,load_RES_ID) /= 0) &
      call runend('readinp_resu : Error getting netCDF ID for dimension '//TRIM(load_RES_name))
  if( nf90_get_var(ncID,load_RES_ID,load_dep,start=(/1,1,nt_dep/)) /= 0) &
      call runend('readinp_resu : Error reading variable '//TRIM(load_RES_name))
  !
  !                                  Class load
  !
  do ic = 1,nc_dep
     if( nf90_inq_varid(ncID,clas_RES_name(ic),clas_RES_ID(ic)) /= 0) &
        call runend('readinp_resu : Error getting netCDF ID for dimension '//TRIM(clas_RES_name(ic)))
     if( nf90_get_var(ncID,clas_RES_ID(ic),work_dep,start=(/1,1,nt_dep/)) /= 0) &
        call runend('readinp_resu : Error reading variable '//TRIM(clas_RES_name(ic)))
     clas_dep(1:nx_dep,1:ny_dep,ic) = work_dep(1:nx_dep,1:ny_dep)
  end do
  !
  if( nf90_close(ncID) /= 0) call runend('readinp_resu :  Error in nf90_close ')
  !
  !*** Reads the mesh (from the dbs meteo file)
  !
  call get_dbs_property_cha(ludbsname,'COORDINATES',coord_sys,istat,message)
  if(TRIM(coord_sys).ne.'LON-LAT') call runend('readinp_resu : UTM coordinates not implemented')
  !
  call get_dbs_dimension(ludbsname,'lon',nx,istat,message)
  call get_dbs_dimension(ludbsname,'lat',ny,istat,message)
  call get_dbs_dimension(ludbsname,'alt',nz,istat,message)
  !
  allocate(lon(nx))
  allocate(lat(ny))
  allocate(zmodel(nz))
  !
  !       Memory allocation
  !
  allocate(Hm1  (nx,ny))
  allocate(load (nx,ny))
  allocate(ust  (nx,ny))
  allocate(smoi (nx,ny))
  allocate(soil (nx,ny))
  allocate(vegfra(nx,ny))
  allocate(rmol (nx,ny))
  allocate(densi(nx,ny))
  !
  allocate(emis(nx,ny,nc))
  allocate(emis_mas(nx,ny,nc))
  allocate(clas(nx,ny,nc_dep))
  allocate(tust(nx,ny,nc_dep))
  !
  load(:,:)       = 0.0d0
  clas(:,:,:)     = 0.0d0
  emis_mas(:,:,:) = 0.0d0
  !
  call get_dbs_value_dimension(ludbsname,'lon',lon,nx,istat,message)
  call get_dbs_value_dimension(ludbsname,'lat',lat,ny,istat,message)
  call get_dbs_value_dimension(ludbsname,'alt',zmodel,nz,istat,message)
  !
  lonmin = lon(1)
  lonmax = lon(nx)
  latmin = lat(1)
  latmax = lat(ny)
  xorigr = lonmin
  yorigr = latmin
  dx     = (lonmax-lonmin)/(nx-1)
  dy     = (latmax-latmin)/(ny-1)
  do ix = 1,nx
    lon(ix) = lonmin + (ix-1)*dx
  end do
  do iy = 1,ny
    lat(iy) = latmin + (iy-1)*dy
  end do
  !
  !*** Compute the scaling factors
  !
  dX1 = Rearth*(dx*pi/180.0_rp)  ! Scaled (dx,dy)
  dX2 = Rearth*(dy*pi/180.0_rp)
  do iy = 1,ny
     colat = (90.0_rp-lat(iy))*pi/180.0_rp   ! colatitude in rad
     Hm1(1:nx,iy) = sin(colat)
  end do
  !
  !*** Interpolates deposit data
  !
  do iy = 1,ny
     !
     !
     foundy = .false.
     do j = 1,ny_dep-1
        if(lat(iy).ge.(yorigr_dep+(j-1)*dy_dep).and.lat(iy).le.(yorigr_dep+(j  )*dy_dep)) then
           foundy = .true.
           iy_dep = j           ! index iy position (bottom point)
        end if
     end do
     if(foundy) then
        t = (lat(iy)-(yorigr_dep+(iy_dep-1)*dy_dep))/dy_dep      ! parameter t in (0,1)
        t = 2.0*t-1.0                                            ! parameter s in (-1,1)
     end if
     !
     !
     do ix = 1,nx
        !
        foundx = .false.
        do i = 1,nx_dep-1
           if(lon(ix).ge.(xorigr_dep+(i-1)*dx_dep).and.lon(ix).le.(xorigr_dep+(i  )*dx_dep)) then
              foundx = .true.
              ix_dep = i         ! index ix position (left point)
           end if
        end do
        if(foundx) then
          s = (lon(ix)-(xorigr_dep+(ix_dep-1)*dx_dep))/dx_dep    ! parameter s in (0,1)
          s = 2.0*s-1.0                                          ! parameter s in (-1,1)
        end if
        !
        !  Point found
        !
        if(foundx.and.foundy) then
           call getshape(s,t,shape)
           !
           load(ix,iy) = load(ix,iy) + &
                         load_dep(ix_dep  ,iy_dep  )*shape(1) + &
                         load_dep(ix_dep+1,iy_dep  )*shape(2) + &
                         load_dep(ix_dep+1,iy_dep+1)*shape(3) + &
                         load_dep(ix_dep  ,iy_dep+1)*shape(4)
           !
           do ic = 1,nc_dep
              clas(ix,iy,ic) = clas(ix,iy,ic) + &
                               clas_dep(ix_dep  ,iy_dep  ,ic)*shape(1) + &
                               clas_dep(ix_dep+1,iy_dep  ,ic)*shape(2) + &
                               clas_dep(ix_dep+1,iy_dep+1,ic)*shape(3) + &
                               clas_dep(ix_dep  ,iy_dep+1,ic)*shape(4)
           end do
           !
        end if
        !
     end do
  end do
  deallocate(load_dep)
  deallocate(clas_dep)
  deallocate(work_dep)
  !
  !*** Reads time related variables
  !
  call get_dbs_dimension(ludbsname,'time',ndt,istat,message)         ! number of time points
  ndt = ndt -1                                                       ! number of time intervals ndt
  !
  call get_dbs_property_int(ludbsname,'START_TIME'    ,idbsbeg ,istat,message)
  call get_dbs_property_int(ludbsname,'END_TIME'      ,idbsend ,istat,message)
  call get_dbs_property_int(ludbsname,'TIME_INCREMENT',idbsincr,istat,message)
  call get_dbs_property_int(ludbsname,'YEAR'          ,idbsyr  ,istat,message)
  call get_dbs_property_int(ludbsname,'MONTH'         ,idbsmo  ,istat,message)
  call get_dbs_property_int(ludbsname,'DAY'           ,idbsday ,istat,message)
  iebeg = idbsbeg
  ieend = idbsend
  irend = idbsend
  !
  allocate(idt_src  (ndt))
  do idt = 1,ndt
     idt_src(idt) = idbsbeg + (idt-1)*idbsincr    ! values of time
  end do
  !
  !*** Number of emission layers
  !
  nz_emission = 1
  do iz = 1,nz
     if(zmodel(iz).le.z_emission) nz_emission = iz
  end do
  !
  !*** Writes the grid for the *.src.nc file
  !
  call wriresu_grid
  !
  !*** log file header
  !
   write(lulog,10)
10 format('  From time   To time    Interval       Total   ',/,&
          '     (s)        (s)      Mass (kg)    Mass (kg) ',/,&
          '  ----------------------------------------------')
  !
  return
  end subroutine readinp_resu
