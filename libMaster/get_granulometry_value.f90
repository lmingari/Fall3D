subroutine get_granulometry_value(fname,word,val,istat,message)
  !**************************************************************************
  !*
  !*    Gets the a granulometric property
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*    character*(*)   word     Property to extract. Possible values are
  !*                             'DIAMETER' (in mm)
  !*                             'DENSITY'
  !*                             'FRACTION'
  !*                             'SPHERICITY'
  !*
  !*    OUTPUT:
  !*    real            val(nc)  Values of the property
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !*************************************************************************
  use kindtype
  implicit none
  !
  character(len=*)  :: fname,message,word
  integer(ip)       :: nc,istat
  real   (rp)       :: val(*)
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen,ipos,ic
  real   (rp)           :: rvoid(4)
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Checks that word exists in the dictionary
  !
  if(word(1:LEN_TRIM(word)).eq.'DIAMETER') then
     ipos = 1
  else if(word(1:LEN_TRIM(word)).eq.'DENSITY') then
     ipos = 2
  else if(word(1:LEN_TRIM(word)).eq.'SPHERICITY') then
     ipos = 3
  else if(word(1:LEN_TRIM(word)).eq.'FRACTION') then
     ipos = 4
  else
     istat = -1
     mymessage = 'get_granulometry_value: word not found in the dictionary'
     message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
     return
  end if
  !
  !***  Opens the file
  !
  open(90,FILE=TRIM(fname),STATUS='old',ERR=100)
  !
  !***  Reads
  !
  read(90,*,ERR=101) nc
  do ic = 1,nc
     read(90,*,ERR=101) rvoid(1),rvoid(2),rvoid(3),rvoid(4)
     val(ic) = rvoid(ipos)
  end do
  close(90)
  !
  !***  Successful end
  !
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_granulometry_value: error opening the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_granulometry_value: error reading the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_granulometry_value
