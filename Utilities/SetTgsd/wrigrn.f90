subroutine wrigrn
  !****************************************************
  !*
  !*    Writes the granulometry file
  !*
  !****************************************************
  use InpOut
  use Master
  implicit none
  !
  integer(ip) :: ic
  !
  !***  Opens the file
  !
  open(lugrn,file=TRIM(lugrnname) ,status='unknown',err=100)
  !
  write(lugrn,'(i5)') nc
  do ic = 1,nc
     write(lugrn,10) diam(ic),rhop(ic),sphe(ic),fc(ic)
  end do
  !
10 format(f10.6,1x,f8.1,1x,f7.3,1x,e16.9)
  !
  close(lugrn)
  return
  !
100 call runend('Error opening Granulometry file '//TRIM(lugrnname))
  return
end subroutine wrigrn
