  subroutine openplumef
  !********************************************************************
  !*
  !*    Gets the names and opens the plume files
  !*
  !********************************************************************
  use KindType
  use InpOut
  implicit none
  !
  logical     :: found
  integer(ip) :: ilen
  !
  !***  Sets file names
  !
  ilen = LEN_TRIM(lusrcname)
  found = .false.
  do while(.not.found)     ! Remove file extension
     !
     if(lusrcname(ilen:ilen).eq.'.') then
        found = .true.
     else
        ilen = ilen - 1
        if(ilen.eq.0) call runend('openplumef: unable to find filenames')
     end if
  end do
  !
  luresname = lusrcname(1:ilen)//'Plume.res'
  lumasname = lusrcname(1:ilen)//'Plume.mass'
  luheiname = lusrcname(1:ilen)//'Plume.height'
  lutemname = lusrcname(1:ilen)//'Plume.tempe'
  !
  !***  Opens files
  !
  open(lures,file=TRIM(luresname),status='unknown',err=103)
  open(lumas,file=TRIM(lumasname),status='unknown',err=104)
  open(luhei,file=TRIM(luheiname),status='unknown',err=105)
  open(lutem,file=TRIM(lutemname),status='unknown',err=106)
  !
  write(lulog,3) TRIM(luresname)
3 format(2x,'Plume results  file   : ',a)
  write(lulog,4) TRIM(lumasname)
4 format(2x,'Plume mass     file   : ',a)
  write(lulog,5) TRIM(luheiname)
5 format(2x,'Plume height   file   : ',a)
  write(lulog,6) TRIM(lutemname)
6 format(2x,'Plume temper   file   : ',a)
  !
  return
  !
103 call runend('Error opening file '//TRIM(luresname))
104 call runend('Error opening file '//TRIM(lumasname))
105 call runend('Error opening file '//TRIM(luheiname))
106 call runend('Error opening file '//TRIM(lutemname))
  !
  end subroutine openplumef
