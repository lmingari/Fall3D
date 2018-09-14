subroutine get_dbs_value_plane(fname,timesec,word,z,nx,ny,var,endsec,istat, &
     message)
  !**************************************************************************
  !*
  !*    Gets a variable from the netCDF database on a z-plane
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the database file
  !*    integer         timesec  Time (in seconds) at which data is extracted
  !*                             Time origin is date at 0000UTC
  !*    character*(*)   word     Variable to select. Possibilities are:
  !*
  !*        Time-independent 2D var: XLON, XLAT, UTM_X, UTM_Y, UTM_ZONE, ZONE
  !*        Time-dependent   2D var: PBLH,  UST,  TOPG,  PRATE
  !*        Time-dependent   3D var:    U,    V,     W,     T,       TV, P
  !*
  !*    real  z       Height at which variable is evaluated (terrain following)
  !*
  !*    OUTPUT:
  !*    real      var (nx,ny)   Variable at height z
  !*    integer   endsec        Time (in seconds) at until which data is valid
  !*                            Time origin is date at 0000UTC
  !*    integer         istat   -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !**************************************************************************
  use kindtype
  use netcdf
  implicit none
  !
  character(len=*)      :: fname,message,word
  integer(ip)           :: timesec,nx,ny,endsec,istat
  real(rp)              :: var(nx,ny),z
  !
  character(len=s_mess)   :: mymessage
  integer(ip)             :: ilen,varcase,nt,nz,it,itime,ialt1,ialt2,iz
  integer(ip)             :: ncID,wordID,altID,timeID
  real   (rp)             :: s
  integer(ip),allocatable :: time(:)
  real   (rp),allocatable :: alt(:),work1(:,:),work2(:,:)
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Checks that word exists in the dictionary and selects the case
  !
  if(word(1:LEN_TRIM(word)).eq.'XLON'.or. &
       word(1:LEN_TRIM(word)).eq.'XLAT'.or. &
       word(1:LEN_TRIM(word)).eq.'UTM_X'.or. &
       word(1:LEN_TRIM(word)).eq.'UTM_Y'.or. &
       word(1:LEN_TRIM(word)).eq.'UTM_ZONE'.or. &
       word(1:LEN_TRIM(word)).eq.'LDU'.or. &
       word(1:LEN_TRIM(word)).eq.'SOIL'.or. &
       word(1:LEN_TRIM(word)).eq.'VEGFRA'.or. &
       word(1:LEN_TRIM(word)).eq.'TOPG') then
     varcase = 0
  else if( &
       word(1:LEN_TRIM(word)).eq.'PBLH'.or. &
       word(1:LEN_TRIM(word)).eq.'UST'.or. &
       word(1:LEN_TRIM(word)).eq.'SMOI'.or. &
       word(1:LEN_TRIM(word)).eq.'RMOL'.or. &
       word(1:LEN_TRIM(word)).eq.'ZNT'.or. &
       word(1:LEN_TRIM(word)).eq.'SPD10'.or. &
       word(1:LEN_TRIM(word)).eq.'PRATE'.or. &
       word(1:LEN_TRIM(word)).eq.'L') then
     varcase = 1
  else if( &
       word(1:LEN_TRIM(word)).eq.'U'.or. &
       word(1:LEN_TRIM(word)).eq.'V'.or. &
       word(1:LEN_TRIM(word)).eq.'W'.or. &
       word(1:LEN_TRIM(word)).eq.'T'.or. &
       word(1:LEN_TRIM(word)).eq.'TP'.or. &
       word(1:LEN_TRIM(word)).eq.'P'.or. &
       word(1:LEN_TRIM(word)).eq.'QV'.or. &
       word(1:LEN_TRIM(word)).eq.'RHO'.or. &
       word(1:LEN_TRIM(word)).eq.'TV') then
     varcase = 2
  else
     istat = -1
     mymessage = 'get_dbs_value_plane: word '//TRIM(word)// &
          ' not found in the dictionary'
     message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
     return
  end if
  !
  !*** Gets database values
  !
  call get_dbs_property_int(fname,'NT',nt,istat,message)
  if(istat /= 0) return
  call get_dbs_property_int(fname,'NZ',nz,istat,message)
  if(istat /= 0) return
  !
  !*** Reads data
  !
  if( nf90_open(TRIM(fname),NF90_NOWRITE, ncID) /= 0 ) goto 101
  !
  !*** Read time intervals
  !
  allocate(time(nt))
  if( nf90_inq_varid(ncID,'time',timeID) /= 0) goto 102
  if( nf90_get_var(ncID,timeID,time) /= 0) goto 103
  if((timesec.lt.time(1)).or.(timesec.gt.time(nt))) goto 104
  !
  itime  = nt
  endsec = time(nt)
  do it=1,nt-1
     if( (timesec.ge.time(it)).and.(timesec.lt.time(it+1)) ) then
        itime  = it
        endsec = time(it+1)
     end if
  end do
  deallocate(time)
  !
  !*** CASE selection
  !
  select case(varcase)
     !
  case(0)
     !
     !***    Reads time-independent 2D data (z=0)
     !
     if( nf90_inq_varid(ncID,word,wordID) /= 0) goto 102
     if( nf90_get_var(ncID,wordID,var,start=(/1,1/)) /= 0) goto 103
     !
  case(1)
     !
     !***    Reads time-dependent 2D data (z=0)
     !
     if( nf90_inq_varid(ncID,word,wordID) /= 0) goto 102
     if( nf90_get_var(ncID,wordID,var,start=(/1,1,itime/)) /= 0) goto 103
     !
  case(2)
     !
     !***    Reads time-dependent 3D data.
     !***    It is necessary to have the zlayers to interpolate
     !
     allocate(alt(nz))
     if( nf90_inq_varid(ncID,'alt',altID) /= 0) goto 102
     if( nf90_get_var(ncID,altID,alt) /= 0) goto 103
     !
     if(z.le.alt(1)) then
        ialt1 = 1
        ialt2 = 1
        s     = 1.0_rp
     else if(z.ge.alt(nz)) then
        ialt1 = nz
        ialt2 = nz
        s     = 1.0_rp
     else
        do iz = 1,nz-1
           if( (z.ge.alt(iz)).and.(z.le.alt(iz+1)) ) then
              ialt1 = iz
              ialt2 = iz+1
              s     = (z-alt(iz))/(alt(iz+1)-alt(iz))
           end if
        end do
     end if
     !
     allocate(work1(nx,ny))
     allocate(work2(nx,ny))
     if(nf90_inq_varid(ncID,word,wordID) /= 0) goto 102
     if(nf90_get_var(ncID,wordID,work1,start=(/1,1,ialt1,itime/))/= 0) goto 103
     if(nf90_get_var(ncID,wordID,work2,start=(/1,1,ialt2,itime/))/= 0) goto 103
     var = (1.0_rp-s)*work1 + s*work2
     !
     deallocate(alt)
     deallocate(work1)
     deallocate(work2)
     !
  end select
  !
  !
  !***  Successful end
  !
  if( nf90_close(ncID) /= 0) goto 105
  return
  !
  !***  List of errors
  !
101 istat = -1
  mymessage = 'get_dbs_value_plane : Error in nf90_open for file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
102 istat = -1
  mymessage = 'get_dbs_value_plane : Error getting netCDF ID for variable '// &
       TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
103 istat = -1
  mymessage = 'get_dbs_value_plane : Error reading variable '//TRIM(word)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
104 istat = -1
  mymessage = 'get_dbs_value_plane: time not in the database interval'
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
105 istat = -1
  mymessage = 'get_dbs_value_plane : Error in nf90_close '
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_dbs_value_plane
