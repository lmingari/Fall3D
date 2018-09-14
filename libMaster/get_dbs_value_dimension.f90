subroutine get_dbs_value_dimension(fname,word,x,nx,istat,message)
  !**************************************************************************
  !*
  !*    Gets a dimension from the netCDF database
  !*
  !*    INPUT:
  !*    character*(*) fname  Name of the database file
  !*    character*(*) word   Variable to select. Possibilities are: lon,lat,
  !*                         alt,diameter,density,sphericity
  !*
  !*    OUTPUT:
  !*    real            x (nx)     Dimension
  !*    integer         istat      -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message    Exit message
  !*
  !*************************************************************************
  use kindtype
  use netcdf
  implicit none
  character(len=*)      :: fname,message,word
  integer(ip)           :: nx,istat
  real(rp)              :: x(nx)
  !
  character(len=s_mess)   :: mymessage
  integer(ip)             :: ilen,ncID,wordID
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Checks that word exists in the dictionary
  !
  if(  word(1:LEN_TRIM(word)) /= 'lon'.and. &
       word(1:LEN_TRIM(word)) /= 'lat'.and. &
       word(1:LEN_TRIM(word)) /= 'alt'.and. &
       word(1:LEN_TRIM(word)) /= 'diameter'.and. &
       word(1:LEN_TRIM(word)) /= 'density'.and. &
       word(1:LEN_TRIM(word)) /= 'sphericity') then
     istat = -1
     mymessage = 'get_dbs_value_dimension: word '//TRIM(word)// &
          ' not found in the dictionary'
     message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
     return
  end if
  !
  !*** Reads data
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) goto 101
  if( nf90_inq_varid(ncID,word,wordID) /= 0) goto 102
  if( nf90_get_var(ncID,wordID,x,start=(/1/)) /= 0) goto 103
  !
  !***  Successful end
  !
  if( nf90_close(ncID) /= 0) goto 104
  return
  !
  !***  List of errors
  !
101 istat = -1
  mymessage = 'get_dbs_value_dimension : Error in nf90_open for file '// &
       TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
102 istat = -1
  mymessage = 'get_dbs_value_dimension : Error getting netCDF ID for '// &
       'variable '//TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
103 istat = -1
  mymessage = 'get_dbs_value_dimension : Error reading variable '//TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
104 istat = -1
  mymessage = 'get_dbs_value_dimension : Error in nf90_close '
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_dbs_value_dimension
