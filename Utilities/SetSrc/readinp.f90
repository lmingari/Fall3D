   subroutine readinp
  !********************************************************
  !*
  !*    Gets properties from the input file. Two different
  !*    modes can exist:
  !*     1- eruption
  !*     2- resuspension
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  character(len=s_mess) :: message,cvoid
  integer(ip)           :: istat
  real   (rp)           :: work(2)
  !
  !*** Determine the mode (source type)
  !
  call get_input_cha &
       (luinpname,'SOURCE','SOURCE_TYPE',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  if(TRIM(cvoid).eq.'POINT') then
     type_source = 'POINT'
  else if(TRIM(cvoid).eq.'SUZUKI') then
     type_source = 'SUZUKI'
  else if(TRIM(cvoid).eq.'PLUME') then
     type_source = 'PLUME'
     call openplumef                            !  Opens files for plume output
  else if(TRIM(cvoid).eq.'RESUSPENSION') then
     type_source = 'RESUSPENSION'
  else
     call runend('Incorrect source type')
  end if
  !
  !*** Determine the aggregation model
  !
  call get_input_cha &
       (luinpname,'AGGREGATION','AGGREGATION_MODEL',cvoid,1,istat,message)
  !
  if(istat.ne.0) then                                  ! AGGREGATION BLOCK not found
      aggregation = .false.
      type_aggr   = 'NONE'
  else if(TRIM(type_source).eq.'RESUSPENSION') then    ! Switch-off aggregation in resuspension mode
      aggregation = .false.
      type_aggr   = 'NONE'
  else if(TRIM(cvoid).eq.'NONE') then
      aggregation = .false.
      type_aggr   = 'NONE'
  else if(TRIM(cvoid).eq.'CORNELL') then
      aggregation = .true.
      type_aggr   = 'CORNELL'
  else if(TRIM(cvoid).eq.'PERCENTAGE') then
      aggregation = .true.
      type_aggr   = 'PERCENTAGE'
  else if(TRIM(cvoid).eq.'COSTA') then
      aggregation = .true.
      type_aggr   = 'COSTA'
  else
      call runend('Incorrect aggregation model')
  end if
  !
  !*** Checks for incompatibilities
  !
  if((type_aggr.eq.'COSTA').and.(type_source.ne.'PLUME')) &
     call runend('COSTA aggregation model is compatible only with the PLUME source option')
  !
  !*** Reads the aggregation block
  !
  if (aggregation) then
     call get_input_rea &
          (luinpname,'AGGREGATION','FI_AGGREGATES',work,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     phi_aggr  = work(1)
     diam_aggr = 2.0_rp**(-phi_aggr)  ! diameter in mm
     diam_aggr = 1d-3*diam_aggr       ! diameter in m
!
     call get_input_rea &
          (luinpname,'AGGREGATION','DENSITY_AGGREGATES',work,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     rho_aggr = work(1)
!
    if(type_aggr.eq.'PERCENTAGE') then
       call get_input_rea &
            (luinpname,'AGGREGATION','PERCENTAGE_(%)',work,1,istat,message)
       if(istat.gt.0) call wriwar(message)
       if(istat.lt.0) call runend(message)
       frac_aggr = 1d-2*work(1)                 ! convert to [0,1]
    end if
!
    if(type_aggr.eq.'COSTA') then
       call get_input_rea &
            (luinpname,'AGGREGATION','FRACTAL_EXPONENT',work,1,istat,message)
       if(istat.ne.0) then
          Df = 2.99
       else
          Df = work(1)
       end if
!
       call get_input_rea &
            (luinpname,'PLUME_SOURCE','EXIT_WATER_FRACTION_(%)',work,1,istat,message)
       if(istat.gt.0) call wriwar(message)
       if(istat.lt.0) call runend(message)
       H2O = .true.
       H2O_percentage = 1d-2*work(1)        ! convert to [0,1]
    end if
  end if
  !
  !*** Reads the granulometry and, if necessary, modifies TGSD to account
  !*** for aggregation
  !
  SELECT CASE(type_source)
  case('POINT','SUZUKI','PLUME')
     call readgrn_erup
  case('RESUSPENSION')
     call readgrn_resu
  case default
     call runend('Incorrect source type')
  END SELECT
  !
  !*** Reads the rest of input file
  !
  SELECT CASE(type_source)
  case('POINT','SUZUKI','PLUME')
     call readinp_erup
  case('RESUSPENSION')
     call readinp_resu
  case default
     call runend('Incorrect source type')
  END SELECT
  !
  return
  end subroutine readinp
