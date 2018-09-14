subroutine get_granulometry_name(fname,word,val,nc,istat,message)
  !**************************************************************************
  !*
  !*    Gets a string associated with a class
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*    character*(*)   word     Property to extract. Possible values are
  !*                             'NAME'
  !*
  !*    OUTPUT:
  !*    character*(*)   val(nc)  Values of the property
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message 
  !*
  !*****************************************************************************
  use kindtype
  implicit none
  !
  integer(ip)       :: nc,istat 
  character(len=*)  :: fname,message,word,val(nc)
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen,ipos,ic,ncc
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
  if(word(1:LEN_TRIM(word)).eq.'NAME') then
     ipos = 1
  else
     istat = -1
     mymessage = 'get_granulometry_name: word not found in the dictionary'
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
  read(90,*,ERR=101) ncc
  do ic = 1,nc
     read(90,*,ERR=101) rvoid(1),rvoid(2),rvoid(3),rvoid(4),val(ic)
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
  mymessage = 'get_granulometry_name: error opening the granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_granulometry_name: error reading the granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_granulometry_name
