subroutine get_granulometry_nclass(fname,nc,istat,message)
  !**************************************************************************
  !*
  !*    Gets the number of granulometric classes
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*
  !*    OUTPUT:
  !*    integer         nc       Number of granulometric classes
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !**************************************************************************
  use kindtype
  implicit none
  !
  character(len=*)  :: fname,message
  integer(ip)       :: nc,istat
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Opens the file
  !
  open(90,FILE=fname(1:LEN_TRIM(fname)),STATUS='old',ERR=100)
  !
  !***  Reads
  !
  read(90,*,ERR=101) nc
  close(90)
  !
  !***  Successful end
  !
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_granulometry_nclass: error opening the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_granulometry_nclass: error reading the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_granulometry_nclass
