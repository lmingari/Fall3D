subroutine get_areas_coordinates(fname,name,x,y,dx,dy,nsrc,istat,message)
  !**************************************************************************
  !*
  !*    Gets the points names and coordinatess
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the file
  !*    integer         nsrc     Number of sources
  !*
  !*    OUTPUT:
  !*    character    name(nsrc)   Names
  !*    real            x(nsrc)   X-coordinates
  !*    real            y(nsrc)   Y-coordinates
  !*    real            dx(nsrc)  DX-coordinates
  !*    real            dy(nsrc)  DY-coordinates
  !*    integer         istat         -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message       Exit message
  !*
  !**************************************************************************
  use kindtype
  implicit none
  !
  integer(ip)       :: nsrc,istat
  character(len=*)  :: fname,message
  character(len=*)  :: name(nsrc)
  real   (rp)       :: x(nsrc),y(nsrc),dx(nsrc),dy(nsrc)
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: isrc,ilen
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Opens the file
  !
  open(90,FILE=TRIM(fname),STATUS='unknown',ERR=100)
  !
  !***  Reads
  !
  do isrc = 1,nsrc
     read(90,*,ERR=101) name(isrc),x(isrc),y(isrc),dx(isrc),dy(isrc)
  end do
  !
  !***  Normal end
  !
  close(90)
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_areas_coordinates: error opening the source file '// &
       TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_areas_coordinates: error reading the points file '
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  close(90)
  !
  return
end subroutine get_areas_coordinates
