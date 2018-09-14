subroutine get_dbs_value_volume(fname,timesec,word,nx,ny,nz,var,endsec, &
     istat,message)
  !**************************************************************************
  !*
  !*    Gets a variable from the netCDF database on a volume
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the database file
  !*    integer         timesec  Time (in seconds) at which data is extracted
  !*                             Time origin is date at 0000UTC
  !*    character*(*)   word     Variable to select. Possibilities are:
  !*            Time-dependent   3D var:    U,    V,     W,     T,       TV, P
  !*    OUTPUT:
  !*    real     var (nx,ny,nz)  Variable
  !*    integer  endsec          Time (in seconds) at until which data is valid
  !*                             Time origin is date at 0000UTC
  !*    integer         istat         -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message       Exit message
  !*
  !**************************************************************************
  use kindtype
  use netcdf
  implicit none
  !
  character(len=*)      :: fname,message,word
  integer(ip)           :: timesec,nx,ny,nz,endsec,istat
  real(rp)              :: var(nx,ny,nz)
  !
  character(len=s_mess)   :: mymessage
  integer(ip)             :: ilen,nt,it,itime
  integer(ip)             :: ncID,wordID,timeID
  integer(ip),allocatable :: time(:)
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Checks that word exists in the dictionary
  !
  if(word(1:LEN_TRIM(word)).ne.'U'.and. &
       word(1:LEN_TRIM(word)).ne.'V'.and. &
       word(1:LEN_TRIM(word)).ne.'W'.and. &
       word(1:LEN_TRIM(word)).ne.'T'.and. &
       word(1:LEN_TRIM(word)).eq.'TP'.or. &
       word(1:LEN_TRIM(word)).eq.'P'.or. &
       word(1:LEN_TRIM(word)).eq.'QV'.or. &
       word(1:LEN_TRIM(word)).eq.'RHO'.or. &
       word(1:LEN_TRIM(word)).ne.'TV') goto 101
  !
  !*** Gets database values
  !
  call get_dbs_property_int(fname,'NT',nt,istat,message)
  if(istat /= 0) return
  !
  !*** Reads data
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) goto 102
  !
  !*** Read time intervals
  !
  allocate(time(nt))
  if( nf90_inq_varid(ncID,'time',timeID) /= 0) goto 103
  if( nf90_get_var(ncID,timeID,time) /= 0) goto 104
  if((timesec.lt.time(1)).or.(timesec.gt.time(nt))) goto 105
  !
  itime  = nt
  endsec = time(nt)
  do it=1,nt-1
     if( (timesec.ge.time(it)).and.(timesec.lt.time(it+1)) ) then
        itime  = it
        endsec = time(it+1)
     end if
  end do
  !
  !***  Read data
  !
  if( nf90_inq_varid(ncID,word,wordID) /= 0) goto 103
  if( nf90_get_var(ncID,wordID,var,start=(/1,1,1,itime/)) /= 0) goto 104
  !
  !***  Successful end
  !
  if( nf90_close(ncID) /= 0) goto 105
  return
  !
  !***  List of errors
  !
101 istat = -1
  mymessage = 'get_dbs_value_volume: word '//TRIM(word)// &
       ' not found in the dictionary'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
102 istat = -1
  mymessage = 'get_dbs_value_volume : Error in nf90_open for file '// &
       TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
103 istat = -1
  mymessage = 'get_dbs_value_volume : Error getting netCDF ID for variable '//&
       TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
104 istat = -1
  mymessage = 'get_dbs_value_volume : Error reading variable '//TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
105 istat = -1
  mymessage = 'get_dbs_value_volume: time not in the database interval'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_dbs_value_volume
