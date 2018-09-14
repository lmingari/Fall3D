subroutine readat
  !****************************************************************
  !*
  !*    Reads data from the input file
  !*
  !****************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  use Coordinates
  implicit none
  !
  character(len=s_mess) :: message,cvoid
  integer  (ip)         :: istat
  real     (rp)         :: rvoid
  !
  !***  Reads scalar data from input file
  !
  !
  !***  TIME_UTC block
  !
  if( mpime == root ) then
     call get_input_int(finp,'TIME_UTC','YEAR',ibyr,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( ibyr , 1, root )
  !
  if( mpime == root ) then
     call get_input_int(finp,'TIME_UTC','MONTH',ibmo,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( ibmo , 1, root )
  !
  if( mpime == root ) then
     call get_input_int(finp,'TIME_UTC','DAY',ibdy,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( ibdy , 1, root )
  !
  if( mpime == root ) then
     call get_input_rea(finp,'TIME_UTC','ERUPTION_START_(HOURS_AFTER_00)',rvoid,1,istat,message)
     beg_time = rvoid*3600.0_rp        ! Convert to sec
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( beg_time, 1, root )
  !
  if( mpime == root ) then
     call get_input_rea(finp,'TIME_UTC','RUN_END_(HOURS_AFTER_00)',rvoid,1,istat,message)
     end_time = rvoid*3600.0_rp       ! Convert to sec
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( end_time, 1, root )
  !
  if( mpime == root ) then
     call get_input_cha(finp,'TIME_UTC','RESTART',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  call bcast( cvoid, LEN( cvoid ), root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) then
     restart = .false.
  else
     if(TRIM(cvoid) == 'YES') then
        restart = .true.
     else
        restart = .false.
     end if
  end if
  !
  !*** GRID block
  !
  coord_sys(:) = ' '
  if( mpime == root ) then
     call get_input_cha (finp,'GRID','COORDINATES',coord_sys,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( coord_sys, LEN( coord_sys ), root )
  if(TRIM(coord_sys).ne.'UTM'    .and. &
       TRIM(coord_sys).ne.'LON-LAT') call runend('Inocrrect system of coordinates')
  !
  SELECT CASE(coord_sys)
  case('LON-LAT')
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRID','LON_VENT',xvent_ll,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( xvent_ll, 1, root )
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRID','LAT_VENT',yvent_ll,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( yvent_ll, 1, root )
     !
     call ll2utm(yvent_ll,xvent_ll,UTM_zone,xvent_utm,yvent_utm,WGS_84_DATUM,istat)
     if(istat /= 0) call wriwar('ll2utm : unable to obtain UTM vent coordinates')
     !
  case('UTM')
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRID','X_VENT',xvent_utm,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( xvent_utm, 1, root )
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRID','Y_VENT',yvent_utm,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( yvent_utm, 1, root )
     !
     UTM_zone = ' '
     if( mpime == root ) then
        call get_input_cha(finp,'GRID','UTMZONE',UTM_zone,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( UTM_zone, LEN( UTM_zone ), root )
     !
     call utm2ll(xvent_utm,yvent_utm,UTM_zone,xvent_ll,yvent_ll,WGS_84_DATUM,istat)
     if(istat.ne.0) call wriwar('utm2ll : unable to obtain LON-LAT vent coordinates')
     !
  END SELECT
  !
  !*** SOURCE Block (gravity current)
  !
  cvoid(:) = ''
  if( mpime == root ) then
     call get_input_cha (finp,'GRAVITY_CURRENT','GRAVITY_CURRENT',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat.ne.0) then
     gravity_current = .false.
  else
     call bcast( cvoid, LEN( cvoid ), root )
     if (TRIM(cvoid) == 'NO') then
        gravity_current = .false.
     else if(TRIM(cvoid) == 'YES') then
        gravity_current = .true.
     else
        call runend('Invalid flag GRAVITY_CURRENT')
     end if
  end if
  !
  if(gravity_current) then
     !
     if( mpime == root ) then
        call wriwar('Gravity current model assumes a Plume source. Check consistence with SerSrc')
     end if
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRAVITY_CURRENT','C_FLOW_RATE',c_flow_rate,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( c_flow_rate, 1, root )
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRAVITY_CURRENT','LAMBDA_GRAV',lambda_grav,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( lambda_grav, 1, root )
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRAVITY_CURRENT','K_ENTRAIN',k_entrain,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( k_entrain, 1, root )
     !
     if( mpime == root ) then
        call get_input_rea(finp,'GRAVITY_CURRENT','BRUNT_VAISALA',brunt_vaisala,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( brunt_vaisala, 1, root )
     !
     !*** Set radial_wind_lon/lat equal to LON_VENT/LAT_VENT
     !
     radial_wind_lon = xvent_ll
     radial_wind_lat = yvent_ll
     radial_wind_x   = xvent_utm
     radial_wind_y   = yvent_utm
     !
  end if
  !
  !*** FALL3D block
  !
  cvoid(:) = ' '
  if( mpime == root ) then
     call get_input_cha (finp,'FALL3D','TERMINAL_VELOCITY_MODEL',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( cvoid, LEN( cvoid ), root )
  !
  if(TRIM(cvoid) == 'ARASTOOPOUR') then
     modv = 1
  else if(TRIM(cvoid) == 'GANSER') then
     modv = 2
  else if(TRIM(cvoid) == 'WILSON') then
     modv = 3
  else if(TRIM(cvoid) == 'DELLINO') then
     modv = 4
  else
     if( mpime == root ) then
        write(nlst,*) '                                         '
        write(nlst,*) 'Error: invalid velocity settling model   '
        write(nlst,*) '   Valid word choices are:               '
        write(nlst,*) '   1 = Arastoopour                       '
        write(nlst,*) '   2 = Ganser                            '
        write(nlst,*) '   3 = Wilson                            '
        write(nlst,*) '   4 = Dellino                           '
     end if
     call runend('Invalid velocity settling model')
  end if
  !
  call setpsi(psi,sphe,diam,modv,np)   ! Calculates psi(np) depending on the model
  !
  cvoid(:) = ' '
  if( mpime == root ) then
     call get_input_cha(finp,'FALL3D','VERTICAL_TURBULENCE_MODEL',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( cvoid, LEN( cvoid ), root )
  !
  if(TRIM(cvoid) == 'CONSTANT') then
     modkv = 0
  else if(TRIM(cvoid) == 'SIMILARITY') then
     modkv = 1
  else if(TRIM(cvoid) == 'SURFACE_LAYER') then
     modkv = 2
  else
     if( mpime == root ) then
        write(nlst,*) '                                         '
        write(nlst,*) 'Error: invalid vertical turbulence model '
        write(nlst,*) '   Valid word choices are:               '
        write(nlst,*) '   Constant                              '
        write(nlst,*) '   Similarity                            '
        write(nlst,*) '   Surface_layer                         '
     end if
     call runend('Invalid turbulence model')
  end if
  !
  if(modkv == 0) then
     if( mpime == root ) then
        call get_input_rea(finp,'FALL3D','VERTICAL_DIFFUSION_COEFFICIENT_(M2/S)',rkv0,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( rkv0, 1, root )
  end if
  !
  cvoid(:) = ' '
  if( mpime == root ) then
     call get_input_cha (finp,'FALL3D','HORIZONTAL_TURBULENCE_MODEL',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( cvoid, LEN( cvoid ), root )
  !
  if(TRIM(cvoid) == 'CONSTANT') then
     modkh = 0
  else if(TRIM(cvoid) == 'RAMS') then
     modkh = 1
     call get_input_rea(finp,'FALL3D','RAMS_CS',Csh,1,istat,message)
     if(istat /= 0) call wriwar('RAMS parameter not specified. A value of 0.2275 is assumed')
  else if(TRIM(cvoid) == 'CMAQ') then
     modkh = 2
  else
     if( mpime == root ) then
        write(nlst,*) '                                         '
        write(nlst,*) 'Error: invalid horizontal turbulence model'
        write(nlst,*) '   Valid word choices are:               '
        write(nlst,*) '   Constant                              '
        write(nlst,*) '   RAMS                                  '
        write(nlst,*) '   CMAQ                                  '
     end if
     call runend('Invalid turbulence model')
  end if
  !
  if(modkh == 0) then
     if( mpime == root ) then
        call get_input_rea(finp,'FALL3D','HORIZONTAL_DIFFUSION_COEFFICIENT_(M2/S)',rkh0,1,istat,message)
     end if
     call bcast( istat, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     call bcast( rkh0, 1, root )
  end if
  !
  cvoid(:) = ' '
  if( mpime == root ) then
     call get_input_cha (finp,'FALL3D','WET_DEPOSITION',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat.ne.0) then
     wetdeposition = .false.
  else
     call bcast( cvoid, LEN( cvoid ), root )
     if(TRIM(cvoid) == 'YES') wetdeposition = .true.
  end if
  !
  !*** AGGREGATION block 
  !
  call readat_aggr
  !
  return
end subroutine readat
  !
  !
  !
subroutine readat_aggr
  !****************************************************************
  !*
  !*    Reads aggregation block from the input file
  !*
  !****************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  implicit none
  !
  character(len=s_mess) :: message,cvoid
  integer  (ip)         :: istat,ic
  real     (rp)         :: rvoid,diam_aggr
  !
  !*** Checks if aggregation during transport must be switched on
  !
  if( mpime == root ) then
     call get_input_cha(finp,'AGGREGATION','AGGREGATION_IN_TRANSPORT',cvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  call bcast( cvoid, 1, root )
  if(istat > 0) call wriwar(message)
  if( TRIM(cvoid).eq.'yes'.or.TRIM(cvoid).eq.'YES') aggregation = .true.
  !
  !*** Check if an aggregate class exists and, if yes, read Vset_factor
  !*** to modify the settling velocity of aggregates
  !
  if(classname(np).eq.'aggregate') then
     !
     if( mpime == root ) then
        call get_input_rea(finp,'AGGREGATION','VSET_FACTOR',rvoid,1,istat,message)
     end if
     call bcast( istat, 1, root )
     call bcast( rvoid, 1, root )
     if(istat > 0) call wriwar(message)
     if(istat < 0) then
        vset_factor = 1.0     ! Default value for vset factor
     else
        vset_factor = rvoid
     end if
     !
  end if
  ! 
  !*** aggregation during transport?
  !
  if(.not.aggregation) then 
     return
  else
     if(classname(np).ne.'aggregate') call runend('Aggregated class does not exist in the grn file')
  end if
  !
  !*** Index of the first aggregating class
  !*** Aggregation occurrs at classes from iaggr:np-1
  !
  diam_aggr = diam(np)
  iaggr = 0
  do ic = np-1,1,-1
     if(diam_aggr.gt.diam(ic)) then
        iaggr = ic
     end if
  end do
  if(iaggr.eq.0) iaggr = np-1  ! ensure at least 1 aggregating class
  !
  if( mpime == root ) then
     call get_input_rea(finp,'AGGREGATION','FRACTAL_EXPONENT',rvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  call bcast( rvoid, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) then
     Df = 3.0             ! Default value for fractal exponent
  else
     Df = rvoid
  end if
  !
  if( mpime == root ) then
     call get_input_rea(finp,'AGGREGATION','N_DT_AGGREGATION',rvoid,1,istat,message)
  end if
  call bcast( istat, 1, root )
  call bcast( rvoid, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) then
     ndtagr = 1           ! Default value for aggregation time step
  else
     ndtagr = int(rvoid)
  end if
  !
  return
end subroutine readat_aggr
