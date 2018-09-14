  subroutine readgrn_erup
  !********************************************************
  !*
  !*    Gets granulometry for the eruption mode
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  character(len=s_mess) :: message,cvoid
  character(len=2     ) :: ext
  integer(ip)           :: istat,ic
  real   (rp)           :: work(2),fca
  !
  !*** Get the number of particles in the TGSD file
  !
  call get_granulometry_nclass(lutgsname,npart,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  !*** SO2 transport
  !
  call get_input_cha &
       (luinpname,'AEROSOLS','SO2',cvoid,1,istat,message)
  if(istat.ne.0) then
      SO2 = .false.
  else if(TRIM(cvoid).eq.'YES') then
      call get_input_rea &
           (luinpname,'AEROSOLS','SO2_PERCENTAGE_(%)',work,1,istat,message)
      if(istat.gt.0) call wriwar(message)
      if(istat.lt.0) call runend(message)
      SO2 = .true.
      SO2_percentage = 1d-2*work(1)        ! convert to [0,1]
  end if
  !
  !*** Determinates the total number of particles
  !
  if(aggregation) npart = npart + 1
  !
  !*** Determinates the number of gas species
  !
  ngas = 0
  if(H2O) ngas = ngas + 1
  if(SO2) ngas = ngas + 1
  !
  !*** Determinates the total number of classes and allocates
  !
  nc = npart + ngas
  !
  allocate(fc  (nc))
  allocate(rhop(nc))
  allocate(diam(nc))
  allocate(sphe(nc))
  allocate(psi (nc))
  allocate(classname(nc))
  fc   = 0.0_rp
  rhop = 0.0_rp
  diam = 0.0_rp
  sphe = 0.0_rp
  psi  = 0.0_rp
  classname = ' '
  !
  !*** Reads particle properties
  !
  call get_granulometry_value(lutgsname,'DIAMETER',diam,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  diam(1:npart)= diam(1:npart)/1d3    ! convert to m
  call get_granulometry_value(lutgsname,'DENSITY',rhop,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  call get_granulometry_value(lutgsname,'SPHERICITY',sphe,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  call get_granulometry_value(lutgsname,'FRACTION',fc,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  !*** Aggregate class properties
  !
  if(aggregation) then
    diam(npart) = diam_aggr
    rhop(npart) = rho_aggr
    sphe(npart) = 1.0_rp
    fc  (npart) = 0.0_rp
  end if
  !
  !*** Aerosol properties. For the gas aerosol components values for density, rho and
  !***  shape are assumed. Approximation can be justifyed by the low
  !***  mass fraction. This is used to generate a column with gas
  !***  distribuion (e.g. H2O for aggregation or SO2)
  !
  if(ngas.gt.0) then
    diam(npart+1:npart+ngas) = 1d-6      ! 1 mic assumed
    rhop(npart+1:npart+ngas) = 1d3
    sphe(npart+1:npart+ngas) = 1.0_rp
    if(ngas.eq.1) then
       if(H2O) fc(npart+1) = H2O_percentage
       if(SO2) fc(npart+1) = SO2_percentage
    else if(ngas.eq.2) then
       fc(npart+1) = H2O_percentage
       fc(npart+2) = SO2_percentage
    end if
  end if
  !
  !*** Determinates class names
  !
  do ic = 1,npart
     ext = '00'
     if(ic.lt.10) then
        write(ext(2:2),'(i1)') ic
     else
        write(ext(1:2),'(i2)') ic
     end if
     classname(ic) = 'class-'//ext
  end do
  if(aggregation) classname(npart) = 'aggregate'
  !
  if(ngas.eq.1) then
     if(H2O) classname(npart+1) = 'H2O'
     if(SO2) classname(npart+1) = 'SO2'
  else if(ngas.eq.2) then
     classname(npart+1) = 'H2O'
     classname(npart+2) = 'SO2'
  end if
  !
  !*** Computes the averaged particle density (without taking into account aggregates)
  !
  rhomean = 0.0_rp
  if(aggregation) then
     do ic = 1,npart-1
        rhomean = rhomean + rhop(ic)*fc(ic)
     end do
  else
     do ic = 1,npart
        rhomean = rhomean + rhop(ic)*fc(ic)
     end do
  end if
  !
  !*** Aggregation
  !
  SELECT CASE(type_aggr)
  case('NONE')
    continue
    !
  case('PERCENTAGE')
    !
    !   Computes aggregation according to a percentage
    !   All classes below diam_aggr are reduced with
    !   a fixed user-defined percentage
    !
    fca = 0.0_rp            ! fraction of mass in the aggregate class
    do ic = 1,npart-1
       if(diam_aggr.gt.diam(ic)) then
          fca    = fca + frac_aggr*fc(ic)
          fc(ic) = (1.0_rp-frac_aggr)*fc(ic)
       end if
    end do
    fc(npart) = fca
    !
  case('CORNELL')
    !
    !   Computes aggregation according to the Cornell model
    !   The aggregate class is made of
    !       50% of particles with 4<phi<3
    !       75% of particles with 4<phi<5
    !       90% of particles with phi>5
    !
    do ic = 1,npart-1
      if(diam(ic).le.0.000125_rp.and.diam(ic).gt.0.0000625_rp) then           ! 4<phi<3
         fc(npart) = fc(npart) + 0.5_rp*fc(ic)
         fc(ic   ) = 0.5_rp*fc(ic)
       else if(diam(ic).le.0.0000625_rp.and.diam(ic).ge.0.00003125_rp) then    ! 4<phi<5
         fc(npart) = fc(npart) + 0.75_rp*fc(ic)
         fc(ic   ) = 0.25_rp*fc(ic)
       else if(diam(ic).lt.0.00003125_rp) then                              ! phi>5
         fc(npart) = fc(npart) + 0.9_rp*fc(ic)
         fc(ic   ) = 0.1_rp*fc(ic)
       end if
    end do
    !
  case('COSTA')
    !
    !  Computed later depending on the plume model (time dependent)
    !
    continue
    !
  END SELECT
  !
  !*** Writes the file with the modifyed granulometry (aggregation + volatiles)
  !
  open(lugrn,file=TRIM(lugrnname) ,status='unknown',err=100)
  write(lugrn,'(i5)') nc
  do ic = 1,nc
     write(lugrn,10) 1d3*diam(ic),rhop(ic),sphe(ic),fc(ic),TRIM(classname(ic))
10  format(f10.6,1x,f8.1,1x,f7.3,1x,e16.9,2x,a)
  end do
  close(lugrn)
  return
100 call runend('Error opening Granulometry file '//TRIM(lugrnname))
  !
  return
  end subroutine readgrn_erup
