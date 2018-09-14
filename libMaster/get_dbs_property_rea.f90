subroutine get_dbs_property_rea(fname,attr,value,istat,message)
  !**************************************************************************
  !*
  !*    Gets a real global attribute from the database netCDF file
  !*
  !*    INPUT:
  !*    character(*)   fname    Name of the database file
  !*    character(*)   attr     word name to select. Possible values are
  !*
  !*
  !*    OUTPUT:
  !*    real           value
  !*    integer        istat    -1 ERROR  0 OK  1 WARNING
  !*    character(*)   message  Exit message
  !*
  !***************************************************************************
  use kindtype
  use netcdf
  implicit none
  !
  character(len=*)      :: fname,attr,message
  integer(ip)           :: istat
  real(rp)              :: value
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen
  integer(ip)           :: ncID
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) goto 101
  !
  !*** Reads the attribute
  !
  if( nf90_get_att(ncID, NF90_GLOBAL, attr, value) /= 0 ) goto 102
  !
  !***  Closes
  !
  if( nf90_close(ncID) /= 0) goto 103
  !
  !***  Successful end
  !
  return
  !
  !***  List of errors
  !
101 istat = -1
  mymessage = 'get_dbs_property_rea : Error in nf90_open for file '// &
       TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
102 istat = -1
  mymessage = 'get_dbs_property_rea : Error in f90_get_att for attribute '// &
       TRIM(attr)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
103 istat = -1
  mymessage = 'get_dbs_property_rea : Error in nf90_close '
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_dbs_property_rea
