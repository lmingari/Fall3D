subroutine get_source_nsrc(fname,timesec,nsrc,endsec,istat,message)
  !**************************************************************************
  !*
  !*    Gets the number of sources
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the source file
  !*    integer         timesec  Time (in seconds) at which data is extracted
  !*                             Time origin is date at 0000UTC
  !*
  !*    OUTPUT:
  !*    integer         nsrc     Number of sources
  !*    integer         endsec   Time (in seconds) at until which data is valid
  !*                             Time origin is date at 0000UTC
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !**************************************************************************
  use kindtype
  implicit none
  !
  character(len=*)  :: fname,message
  integer(ip)       :: nsrc,istat,timesec,endsec
  !
  logical               :: go_on
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen,itime1,itime2,ns,nc,is,ic
  real   (rp)           :: rvoid
  integer(rp)           :: ivoid
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Opens the file
  !
  open(90,FILE=fname(1:LEN_TRIM(fname)),STATUS='unknown',ERR=100)
  !
  !***  Reads until the required time
  !
  go_on = .true.
  do while(go_on)
     read(90,*,ERR=101) itime1,itime2
     read(90,*,ERR=101) ns    ,nc
     if((timesec.ge.itime1).and.(timesec.lt.itime2)) then
        go_on = .false.
     else
        read(90,*,ERR=101) rvoid            ! MFR
        do is = 1,ns
           read(90,*,ERR=101) rvoid,rvoid,rvoid,ivoid,ivoid,ivoid,(rvoid,ic=1,nc)
        end do
     end if
  end do
  !
  !***  Normal end
  !
  endsec = itime2
  nsrc   = ns
  close(90)
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_source_nsrc: error opening the source file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_source_nsrc: error reading the source file at time'
  ilen = LEN_TRIM(mymessage)
  write(mymessage(ilen+1:ilen+6),'(i6)') timesec
  ilen = LEN_TRIM(mymessage)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_source_nsrc
