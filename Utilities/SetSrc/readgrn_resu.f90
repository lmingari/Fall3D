  subroutine readgrn_resu
  !********************************************************
  !*
  !*    Gets granulometry for the resuspension mode
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  use Deposit
  implicit none
  !
  character(len=s_mess) :: message
  character(len=2     ) :: ext
  integer(ip)           :: istat,ic
  !
  !**** Reads properties of the nc_dep classes in the deposit file (netCDF)
  !
  call get_dbs_dimension(ludepname,'np',nc_dep,istat,message)    ! number of classes (in the deposit)
  !
  allocate(rhop_dep(nc_dep))
  allocate(diam_dep(nc_dep))
  allocate(sphe_dep(nc_dep))
  allocate(psi_dep (nc_dep))
  allocate(vset_dep(nc_dep))
  !
  call get_dbs_value_dimension(ludepname,'diameter'  ,diam_dep,nc_dep,istat,message)  ! in mm
  diam_dep = 1d-3*diam_dep                                                            ! in m
  call get_dbs_value_dimension(ludepname,'density'   ,rhop_dep,nc_dep,istat,message)
  call get_dbs_value_dimension(ludepname,'sphericity',sphe_dep,nc_dep,istat,message)
  !
  !***  Maximum resuspended size (user-defined)
  !
  call get_input_rea(luinpname,'RESUSPENSION_SOURCE','MAX_RESUSPENSION_SIZE_(MIC)',diam_max,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  diam_max = diam_max*1d-6                      ! max diameter in m (read in microns)
  diam_max = max(diam_max,diam_dep(nc_dep))     ! there must be at least one saltating particle class
  !
  !*** Determines the number of classes that can be transported. Note that this can be lower
  !*** than the total number of classes in the deposit depending on diam_max
  !
  ngas  = 0
  npart = 0
  do ic = 1,nc_dep
     if(diam_dep(ic).le.diam_max) npart = npart + 1
  end do
  if(npart.eq.0) call runend('The deposit has no particles lower or equal than max size')
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
  !*** Class properties
  !
  do ic = 1,nc
     diam(ic) = diam_dep(nc_dep-nc+ic)
     rhop(ic) = rhop_dep(nc_dep-nc+ic)
     sphe(ic) = sphe_dep(nc_dep-nc+ic)
     fc  (ic) = 1.0_rp/nc
     !
     ext = '00'
     if(ic.lt.10) then
        write(ext(2:2),'(i1)') ic
     else
        write(ext(1:2),'(i2)') ic
     end if
     classname(ic) = 'class-'//ext
  end do
  !
  !*** Computes the averaged particle density (without taking into account aggregates)
  !
  rhomean = 0.0_rp
  do ic = 1,npart
     rhomean = rhomean + rhop(ic)*fc(ic)
  end do
  !
  !*** Writes the file with the effective granulometry (only resuspended particles)
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
  end subroutine readgrn_resu
