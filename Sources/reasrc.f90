  subroutine reasrc
  !***********************************************************************
  !*
  !*    Reads sources from the source file
  !*
  !***********************************************************************
  use KindType
  use InpOut
  use Master
  use Numeric
  use Parallel
  implicit none
  !
  character(len=s_mess) :: message
  integer(ip)           :: istat,iout,isou
  integer(rp), save     :: ipass = 0
  !
  if(ipass == 0) then
     ipass = 1
  else
     call outmem2         ! releases memory for vectors related to nsou
  end if
  !
  !***  Reads nsou at the current time instant
  !
  if( mpime == root ) then
     call get_source_nsrc(fsrc,INT(time),nsou,iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast(  nsou, 1, root )
  !
  !***  Allocates memory for vectors related to nsou (time dependent)
  !
  call getmem2
  !
  !***  Loads vectors
  !
  if( mpime == root ) then
     call get_source_coordinates(fsrc,INT(time),xs,ys,zs,iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast(xs , SIZE(xs), root )
  call bcast(ys , SIZE(ys), root )
  call bcast(zs , SIZE(zs), root )
  !
  if( mpime == root ) then
     call get_source_mesh(fsrc,INT(time),isrc,jsrc,ksrc,iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat > 0) call wriwar(message)
  if(istat < 0) call runend(message)
  call bcast(isrc , SIZE(isrc), root )
  call bcast(jsrc , SIZE(jsrc), root )
  call bcast(ksrc , SIZE(ksrc), root )
  !
  if( mpime == root ) then
      call get_source_value(fsrc,INT(time),tmrat,iout,istat,message)
  end if
  call bcast( istat, 1, root )
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  call bcast( tmrat, SIZE(tmrat), root )
  !
  !***  Set source_time
  !
  call bcast(iout, 1, root )
  source_time = DBLE(iout)
  !
  !*** Checks for inconsistence between mesh and source index. This may occurr if
  !*** a user defines a coarser mesh in the dbs file and forgets to run SetSrc 
  !
  do isou = 1,nsou
     if(isrc(isou).gt.nx) call runend('Inconsistence between dbs and src files. isrc > nx')
     if(jsrc(isou).gt.ny) call runend('Inconsistence between dbs and src files. jsrc > ny')
     if(ksrc(isou).gt.nz) call runend('Inconsistence between dbs and src files. ksrc > nz')
  end do
  !
  return
  end subroutine reasrc
