   subroutine reagrn
  !****************************************************************
  !*
  !*    Reads granulometry from the granulometry file
  !*
  !****************************************************************
  use KindType
  use InpOut
  use Master
  use Parallel
  implicit none
  !
  integer(ip)           :: istat,ic
  character(len=s_mess) :: message
  !
  !***  Reads data
  !
  !***  First number of particles and gas species
  !
  if( mpime == root ) then
     call get_granulometry_nclass(fgrn,nc,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( nc, 1, root )
  !
  if(nc < 1    ) call runend('At least 1 particle class must exist   ')
  if(nc > ncmax) call runend('Too many particle classes, change ncmax')
  !
  !*** Reads particles
  !
  if( mpime == root ) then
     call get_granulometry_value(fgrn,'DIAMETER',diam,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( diam, SIZE( diam ), root )
  diam(1:nc) = diam(1:nc)/1d3           ! Convert diameter to m
  !
  if( mpime == root ) then
     call get_granulometry_value(fgrn,'DENSITY',rhop,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( rhop, SIZE( rhop ), root )
  !
  if( mpime == root ) then
     call get_granulometry_value(fgrn,'SPHERICITY',sphe,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( sphe, SIZE( sphe ), root )
  sphe(1:nc) = min(sphe(1:nc),1.0_rp)
  !
  if( mpime == root ) then
      call get_granulometry_value(fgrn,'FRACTION',fc,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast( fc, SIZE( fc ), root )
  !
  if( mpime == root ) then
     call get_granulometry_name(fgrn,'NAME',classname,nc,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  do ic = 1,nc
     call bcast(classname(ic), LEN( classname(ic) ), root )
  end do
  !
  !*** Gets np (number of particle classes) and ng (number of gas classes) from names
  !
  np = nc
  ng = 0
  do ic = 1,nc
     if( (TRIM(classname(ic)).eq.'H2O').or. &
         (TRIM(classname(ic)).eq.'SO2') ) then
       ng = ng + 1
       np = np - 1
     end if
  end do
  !
  return
  end subroutine reagrn
