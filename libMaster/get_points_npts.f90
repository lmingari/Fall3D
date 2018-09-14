subroutine get_points_npts(fname,npts,istat,message)
  !**************************************************************************
  !*
  !*    Gets the number of tracking points
  !*
  !*
  !*    OUTPUT:
  !*    integer         npts     Number of points
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !***************************************************************************
  use kindtype
  implicit none
  !
  character(len=*)  :: fname,message
  integer(ip)       :: npts,istat
  !
  logical               :: go_on
  character(len=s_mess) :: mymessage
  character(len=1     ) :: cvoid
  integer(ip)           :: ilen,info
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  npts  = 0
  !
  !***  Opens the file
  !
  open(90,FILE=fname(1:LEN_TRIM(fname)),STATUS='unknown',ERR=100)
  !
  !***  Reads until the EOF
  !
  go_on = .true.
  do while(go_on)
     read(90,*,iostat=info) cvoid
     if(info /= 0) then
        go_on = .false.
     else
        npts = npts + 1
     end if
  end do
  !
  !***  End
  !
  close(90)
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_points_npts: error opening the points file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_points_npts
