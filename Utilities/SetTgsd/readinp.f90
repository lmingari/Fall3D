subroutine readinp
  !************************************************************
  !*
  !*    Reads data from the input file
  !*
  !************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  character(len=s_mess) :: message
  integer(ip) :: istat
  real   (rp) :: work(2)
  !
  !***  Reads data input file
  !
  call get_input_cha &
       (luinpname,'GRANULOMETRY','DISTRIBUTION',TYPE_DIST,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  if(TYPE_DIST.eq.'GAUSSIAN') then
     mdist = 0
     ng    = 1
  else if(TYPE_DIST.eq.'BIGAUSSIAN') then
     mdist = 1
     ng    = 2
  else
     call runend('Type of distribution not contemplated')
  end if
  !
  call get_input_int &
       (luinpname,'GRANULOMETRY','NUMBER_OF_CLASSES',nc,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea &
       (luinpname,'GRANULOMETRY','FI_MEAN',fimean,ng,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea &
       (luinpname,'GRANULOMETRY','FI_DISP',fidisp,ng,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  !
  call get_input_rea &
       (luinpname,'GRANULOMETRY','FI_RANGE',work,2,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  fimin = work(1)
  fimax = work(2)
  if(fimin.gt.fimax) then
     fimin = work(2)
     fimax = work(1)
  end if
  !
  call get_input_rea &
       (luinpname,'GRANULOMETRY','DENSITY_RANGE',work,2,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  rhomin = work(1)
  rhomax = work(2)
  !
  call get_input_rea &
       (luinpname,'GRANULOMETRY','SPHERICITY_RANGE',work,2,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  sphemin  = work(1)
  sphemax  = work(2)
  !
  !***  Allocates memory
  !
  call getmem
  !
  return
end subroutine readinp
