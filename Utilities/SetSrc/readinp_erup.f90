  subroutine readinp_erup
  !********************************************************
  !*
  !*    Gets properties from the input file for the
  !*    eruption mode
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  use Coordinates
  use Plume
  implicit none
  !
  character(len=s_mess) :: message,cvoid
  integer(ip)           :: istat,ndt0,ivoid
  real(rp)              :: rvoid,vvoid(100)
  !
  !***  Initialization
  !
  lonmin = 0.0
  lonmax = 0.0
  latmin = 0.0
  latmax = 0.0
  X0 = 0.0
  Y0 = 0.0
  Z0 = 0.0
  !
  !***  Reads data input file
  !
  !***  Gets the number of eruption intervals ndt
  !
  call get_input_npar &
       (luinpname,'TIME_UTC','ERUPTION_START_(HOURS_AFTER_00)',ndt,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  !***  Allocates memory for ndt dependent variables
  !
  allocate(idt_src  (ndt))
  allocate(M0_dt    (ndt))
  allocate(HPlume_dt(ndt))
  allocate(Asuzu_dt (ndt))
  allocate(Lsuzu_dt (ndt))
  allocate(u0_dt    (ndt))
  allocate(T0_dt    (ndt))
  allocate(w0_dt    (ndt))
  !
  !***  TIME BLOCK
  !
  call get_input_rea &
       (luinpname,'TIME_UTC','ERUPTION_START_(HOURS_AFTER_00)',vvoid,min(ndt,100),istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  do idt = 1,ndt
     idt_src(idt) = INT(vvoid(idt)*3600.0)
  end do
  iebeg = idt_src(1)
  !
  call get_input_rea &
       (luinpname,'TIME_UTC','ERUPTION_END_(HOURS_AFTER_00)',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  ieend = INT(rvoid*3600.0)
  do idt = 1,ndt
     if(idt_src(idt).ge.ieend) call runend('Inconsiestency between eruption intervals and eruption end')
  end do
  !
  call get_input_rea &
       (luinpname,'TIME_UTC','RUN_END_(HOURS_AFTER_00)',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  irend = INT(rvoid*3600.0)
  !
  call get_dbs_property_int(ludbsname,'START_TIME',idbsbeg,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  if(iebeg.lt.idbsbeg) call runend('idbsbeg greater than iebeg')
  !
  call get_dbs_property_int(ludbsname,'END_TIME',idbsend,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  if(ieend.gt.idbsend) call runend('idbsend lower than ieend')
  !
  !***  Type of source
  !
  SELECT CASE(type_source)
  case('POINT')
     ns = 1            ! single source
     np = 0
     !
  case('SUZUKI')
     ns = 100           ! number of suzuki sources
     np = 0
     !
  case('PLUME')
     ns = 100          ! number of plume sources (plume+umbrella)
     np = 70           ! number of plume sources (with no umbrella)
     !
  case default
     call runend('Incorrect source type')
  END SELECT
  !
  !***  Reads SOURCE block
  !
  SELECT CASE(type_source)
  case('POINT')
     !
     !*** 1. POINT source case
     !
     call get_input_npar(luinpname,'POINT_SOURCE','HEIGHT_ABOVE_VENT_(M)',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'POINT_SOURCE','HEIGHT_ABOVE_VENT_(M)',Hplume_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) Hplume_dt(ndt0:ndt) = Hplume_dt(ndt0)
     !
     call get_input_cha &
          (luinpname,'POINT_SOURCE','MASS_FLOW_RATE_(KGS)',cvoid,1,istat,message)
     !
     if(TRIM(cvoid).eq.'ESTIMATE-MASTIN') then
        !
        MER_vs_h          = 'ESTIMATE-MASTIN'
        MER_wind_coupling = .false.
        !
        !*** Estimate MFR from height using  MFR = rho*(H/2)**(1/0.241)
        !
        do idt = 1,ndt
           M0_dt(idt) = rhomean*((0.5_rp*Hplume_dt(idt)/1d3)**(1.0_rp/0.241_rp))
        end do
        !
     else if(TRIM(cvoid).eq.'ESTIMATE-DEGRUYTER') then
        !
        MER_vs_h          = 'ESTIMATE-DEGRUYTER'
        MER_wind_coupling = .true.
        !
     else if(TRIM(cvoid).eq.'ESTIMATE-WOODHOUSE') then
        !
        MER_vs_h          = 'ESTIMATE-WOODHOUSE'
        MER_wind_coupling = .true.
        !
     else
        !
        MER_vs_h          = 'NONE'
        MER_wind_coupling = .false.
        !
        call get_input_npar(luinpname,'POINT_SOURCE','MASS_FLOW_RATE_(KGS)',ndt0,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        call get_input_rea (luinpname,'POINT_SOURCE','MASS_FLOW_RATE_(KGS)',M0_dt,ndt0,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        if(ndt0.lt.ndt) M0_dt(ndt0:ndt) = M0_dt(ndt0)
     end if
     !
  case('SUZUKI')
     !
     !*** 2. SUZUKI source case
     !
     call get_input_npar(luinpname,'SUZUKI_SOURCE','HEIGHT_ABOVE_VENT_(M)',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'SUZUKI_SOURCE','HEIGHT_ABOVE_VENT_(M)',Hplume_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) Hplume_dt(ndt0:ndt) = Hplume_dt(ndt0)
     !
     call get_input_npar(luinpname,'SUZUKI_SOURCE','A',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'SUZUKI_SOURCE','A',Asuzu_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) Asuzu_dt(ndt0:ndt) = Asuzu_dt(ndt0)
     !
     call get_input_npar(luinpname,'SUZUKI_SOURCE','L',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea  (luinpname,'SUZUKI_SOURCE','L',Lsuzu_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) Lsuzu_dt(ndt0:ndt) = Lsuzu_dt(ndt0)
     !
     call get_input_cha &
          (luinpname,'SUZUKI_SOURCE','MASS_FLOW_RATE_(KGS)',cvoid,1,istat,message)
     !
     if(TRIM(cvoid).eq.'ESTIMATE-MASTIN') then
        !
        MER_vs_h          = 'ESTIMATE-MASTIN'
        MER_wind_coupling = .false.
        !
        !*** Estimate MFR from height using  MFR = rho*(H/2)**(1/0.241)
        !
        do idt = 1,ndt
           M0_dt(idt) = rhomean*((0.5_rp*Hplume_dt(idt)/1d3)**(1.0_rp/0.241_rp))
        end do
        !
     else if(TRIM(cvoid).eq.'ESTIMATE-DEGRUYTER') then
        !
        MER_vs_h          = 'ESTIMATE-DEGRUYTER'
        MER_wind_coupling = .true.
        !
     else if(TRIM(cvoid).eq.'ESTIMATE-WOODHOUSE') then
        !
        MER_vs_h          = 'ESTIMATE-WOODHOUSE'
        MER_wind_coupling = .true.
        !
     else
        !
        MER_vs_h          = 'NONE'
        MER_wind_coupling = .false.
        !
        call get_input_npar(luinpname,'SUZUKI_SOURCE','MASS_FLOW_RATE_(KGS)',ndt0,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        call get_input_rea (luinpname,'SUZUKI_SOURCE','MASS_FLOW_RATE_(KGS)',M0_dt,ndt0,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        if(ndt0.lt.ndt) M0_dt(ndt0:ndt) = M0_dt(ndt0)
     end if
     !
  case('PLUME')
     !
     !*** 3. PLUME source case
     !
     call get_input_cha &
          (luinpname,'PLUME_SOURCE','SOLVE_PLUME_FOR',solve_plume_for,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     if(TRIM(solve_plume_for).eq.'MFR'.or.TRIM(solve_plume_for).eq.'mfr') then
        call get_input_rea (luinpname,'PLUME_SOURCE','MFR_SEARCH_RANGE',n_MFR,2,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
     else if(TRIM(solve_plume_for).eq.'HEIGHT'.or.TRIM(solve_plume_for).eq.'height') then
        continue
     else
        call runend('Specify the kind of plume solving strategy')
     end if
     !
     call get_input_npar(luinpname,'PLUME_SOURCE','HEIGHT_ABOVE_VENT_(M)',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'PLUME_SOURCE','HEIGHT_ABOVE_VENT_(M)',Hplume_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) Hplume_dt(ndt0:ndt) = Hplume_dt(ndt0)
     !
     call get_input_npar(luinpname,'PLUME_SOURCE','MASS_FLOW_RATE_(KGS)',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'PLUME_SOURCE','MASS_FLOW_RATE_(KGS)',M0_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) M0_dt(ndt0:ndt) = M0_dt(ndt0)
     !
     call get_input_npar(luinpname,'PLUME_SOURCE','EXIT_VELOCIY_(MS)',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'PLUME_SOURCE','EXIT_VELOCIY_(MS)',u0_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) u0_dt(ndt0:ndt) = u0_dt(ndt0)
     !
     call get_input_npar(luinpname,'PLUME_SOURCE','EXIT_TEMPERATURE_(K)',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'PLUME_SOURCE','EXIT_TEMPERATURE_(K)',T0_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) T0_dt(ndt0:ndt) = T0_dt(ndt0)
     !
     call get_input_npar(luinpname,'PLUME_SOURCE','EXIT_WATER_FRACTION_(%)',ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     call get_input_rea (luinpname,'PLUME_SOURCE','EXIT_WATER_FRACTION_(%)',w0_dt,ndt0,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if(ndt0.lt.ndt) w0_dt(ndt0:ndt) = w0_dt(ndt0)
     w0_dt(:) = w0_dt(:)/1d2                        ! Convert from %
     !
  END SELECT
  !
  !***  Reads vent coordinates
  !
  call get_input_cha(luinpname,'GRID','COORDINATES',coord_sys,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  if(TRIM(coord_sys).ne.'UTM'    .and. &
     TRIM(coord_sys).ne.'LON-LAT') call runend('Inocrrect system of coordinates')
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     call get_input_rea(luinpname,'GRID','LON_VENT',X0_LL,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     X0 = X0_LL
     !
     call get_input_rea(luinpname,'GRID','LAT_VENT',Y0_LL,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     Y0 = Y0_LL
     !
     call ll2utm(Y0_LL,X0_LL,UTM_zone,X0_UTM,Y0_UTM,WGS_84_DATUM,istat)
     if(istat.ne.0) call wriwar('ll2utm : unable to obtain UTM vent coordinates')
     !
     !***     Reads the vent height from the dbs or (if not found) from the input file
     !
     rvoid = 0.0
     call get_dbs_value_point(ludbsname,iebeg,'TOPG',X0_LL,Y0_LL,rvoid,Z0,ivoid,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) then
        call wriwar(TRIM(message)//' : Reading vent height from input file')
        call get_input_rea(luinpname,'GRID','VENT_HEIGHT_(M)',Z0,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
     end if
     !
  case('UTM')
     !
     call get_input_rea(luinpname,'GRID','X_VENT',X0_UTM,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     X0 = X0_UTM
     !
     call get_input_rea(luinpname,'GRID','Y_VENT',Y0_UTM,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     Y0 = Y0_UTM
     !
     call get_input_cha(luinpname,'GRID','UTMZONE',UTM_zone,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call utm2ll(X0_UTM,Y0_UTM,UTM_zone,X0_LL,Y0_LL,WGS_84_DATUM,istat)
     if(istat.ne.0) call wriwar('utm2ll : unable to obtain LON-LAT vent coordinates')
     !
     !***     Reads the vent height from the dbs or (if not found) from the input file
     !
     rvoid = 0.0
     call get_dbs_value_point(ludbsname,iebeg,'TOPG',X0_UTM,Y0_UTM,rvoid,Z0,ivoid,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) then
        call wriwar(TRIM(message)//' : Reading vent height from input file')
        call get_input_rea(luinpname,'GRID','VENT_HEIGHT_(M)',Z0,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
     end if
     !
  END SELECT
  !
  !***  Reads the real vent height from the input file. This is done to prevent
  !***  lack of accuracy in case the dbs gives coarse topography
  !
  dZ0 = Z0
  call get_input_rea(luinpname,'GRID','VENT_HEIGHT_(M)',dZ0,1,istat,message)
  if(istat.ne.0) call wriwar('Vent height not specifyed; extracted from topography')
  !
  dZ0 = dZ0 - Z0           ! vent correction
  if(correct_vent) then
     write(lulog,1) Z0,Z0+dZ0,dZ0
1    format(/,&
            2x,'Vent height from topography : ',f8.1,/,&
            2x,'Vent height from input      : ',f8.1,/,&
            2x,'Difference (corrected)      : ',f8.1,/)
  else
     write(lulog,2) Z0,Z0+dZ0,dZ0
2    format(/,&
            2x,'Vent height from topography : ',f8.1,/,&
            2x,'Vent height from input      : ',f8.1,/,&
            2x,'Difference (not used)       : ',f8.1,/)
     dZ0 = 0.0
  end if
  !
  !***  Reads GRID block
  !
  call get_input_int(luinpname,'GRID','NX',nx,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_int(luinpname,'GRID','NY',ny,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     call get_input_rea(luinpname,'GRID','LONMIN',lonmin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea(luinpname,'GRID','LONMAX',lonmax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea(luinpname,'GRID','LATMIN',latmin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea(luinpname,'GRID','LATMAX',latmax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     xorigr = lonmin
     yorigr = latmin
     dx     = (lonmax-lonmin)/(nx-1)
     dy     = (latmax-latmin)/(ny-1)
     !
  case('UTM')
     !
     call get_input_rea(luinpname,'GRID','XMIN',xmin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea(luinpname,'GRID','XMAX',xmax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea(luinpname,'GRID','YMIN',ymin,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     call get_input_rea(luinpname,'GRID','YMAX',ymax,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     !
     xorigr = xmin
     yorigr = ymin
     dx     = (xmax-xmin)/(nx-1)
     dy     = (ymax-ymin)/(ny-1)
     !
  END SELECT
  !
  !***  Gets the number of layers and allocates memory
  !
  call get_input_npar(luinpname,'GRID','ZLAYER_(M)',nz,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  allocate(topg(nx,ny))
  allocate(zmodel(nz))
  !
  !**   Reads topography and zmodel
  !
  call get_dbs_value_plane(ludbsname,iebeg,'TOPG',0.0,nx,ny,topg,ivoid,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) then
     call wriwar(TRIM(message)//' : Assuming zero topography')
     topg = 0.0
  end if
  !
  call get_input_rea(luinpname,'GRID','ZLAYER_(M)',zmodel,nz,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
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
  call setpsi(psi,sphe,diam,modv,nc)   ! Calculates psi(nc) depending on the model
  !
  !*** Allocates memory for the rest of arrays (depending on ns)
  !
  allocate(Xplum(ns))       ! (ns)
  allocate(Yplum(ns))
  allocate(Zplum(ns))
  allocate(Lplum(ns))
  allocate(Hplum(ns))
  allocate(Uplum(ns))
  allocate(Tplum(ns))
  allocate(Dplum(ns))
  allocate(Rplum(ns))
  allocate(Mair (ns))
  allocate(Mplum(nc,ns))
  Mplum(:,:)     = 0.0_rp
  !
  allocate(src(nc,nx,ny,nz))
  allocate(Vair(nz))
  allocate(Tair(nz))
  allocate(Aair(nz))
  allocate(Rair(nz))
  allocate(Nair(nz))
  !
  !*** log file header
  !
  SELECT CASE(type_source)
  case('POINT','SUZUKI')
     write(lulog,10)
10   format('  From time   To time      MER     Column height  Column height     Mass   ',/,&
            '     (s)        (s)      (kg/s)      (m a.vent)     (m a.s.l)       (kg)   ',/,&
            '  -------------------------------------------------------------------------')
  case('PLUME')
     write(lulog,20)
20   format('  From time   To time      MER     Column height  Column height     Mass   ',/,&
            '     (s)        (s)      (kg/s)    NBL (m a.s.l.) Total (m a.s.l)   (kg)   ',/,&
            '  -------------------------------------------------------------------------')
  END SELECT
  !
  return
  end subroutine readinp_erup
