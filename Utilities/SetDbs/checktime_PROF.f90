subroutine checktime_PROF
  !*************************************************************************
  !*
  !*    Outputs the PROF and dbs time coverages and gets the time lag
  !*    (if any)
  !*
  !*************************************************************************
  use kindType
  use MathFun
  use TimeFun
  use Master
  use PROF_nc
  use InpOut
  implicit none
  !
  logical         :: found
  integer(ip)     :: iyr,imo,idy,ihr,imi,ise
  character(len=8):: date_PROF
  !
  !***  Profile file
  !***    (x,o,yo)
  !***    YYYYMMDD
  !***    resto fo file
  !
  open(ludat,FILE=TRIM(ludatname),status='old',ERR=100)
  !
  !***  Reads profile coordinates and checks that are inside the domain
  !
  read(ludat,*,ERR=101) xo_PROF, yo_PROF
  SELECT CASE(coord_sys)
  CASE('LON-LAT')
     if( (xo_PROF.lt.lonp(1,1)).or.(xo_PROF.gt.lonp(nx,1)) ) &
          call runend('checktime_PROF: profile not inside the lon domain')
     if( (yo_PROF.lt.latp(1,1)).or.(yo_PROF.gt.latp(1,ny)) ) &
          call runend('checktime_PROF: profile not inside the lat domain')
  CASE('UTM')
     if( (xo_PROF.lt.xutm(1,1)).or.(xo_PROF.gt.xutm(nx,1)) ) &
          call runend('checktime_PROF: profile not inside the UTM domain in x')
     if( (yo_PROF.lt.yutm(1,1)).or.(yo_PROF.gt.yutm(1,ny)) ) &
          call runend('checktime_PROF: profile not inside the UTM domain in Y')

  END SELECT
  !
  !*** Reads beginning date YYYYMMDD
  !
  read(ludat,*,ERR=101) date_PROF
  !
  ibyr_PROF = stoi1(date_PROF(1:1))*1d3 + &
       stoi1(date_PROF(2:2))*1d2 + &
       stoi1(date_PROF(3:3))*1d1 + &
       stoi1(date_PROF(4:4))*1d0
  ibmo_PROF = stoi1(date_PROF(5:5))*1d1 + &
       stoi1(date_PROF(6:6))*1d0
  ibdy_PROF = stoi1(date_PROF(7:7))*1d1 + &
       stoi1(date_PROF(8:8))*1d0
  !
  !***  Averiguates nt_PROF and sets timesec_PROF(nt_PROF) in sec after 0000UTC
  !
  nt_PROF = 1
  do while(0.ne.1)
     read(ludat,*,ERR=101,END=50) timesec_PROF(nt_PROF),timesec_PROF(nt_PROF+1)
     read(ludat,*,ERR=101,END=50) nz_PROF
     do iz_PROF = 1,nz_PROF
        read(ludat,*,ERR=101,END=50) rvoid
     end do
     nt_PROF = nt_PROF + 1
     if(nt_PROF.eq.(nt_PROF_max-1)) call runend('checktime_PROF: too many profile timesteps')
     !
  end do
50 close(ludat)
  !
  !***  Computes time_PROF in the format YYYYMMDDHHMMSS
  !
  do it_PROF = 1,nt_PROF
     call addtime(ibyr_PROF,ibmo_PROF,ibdy_PROF,0,iyr,imo,idy,ihr,imi,ise,timesec_PROF(it_PROF))
     time_PROF(it_PROF) = 1d10*iyr + 1d8 *imo + 1d6 *idy + 1d4 *ihr + 1d2*imi + ise
  end do
  !
  !***  Writes
  !
  write(lulog,10) time_PROF(1      ),timesec_PROF(1)     , &
       time_PROF(nt_PROF),timesec_PROF(nt_PROF)

10 format('Profile time coverage : ',/, &
       ' From         : ',f15.0,' ( ',f8.0,' sec after 0000UTC)',/, &
       ' To           : ',f15.0,' ( ',f8.0,' sec after 0000UTC)')
  write(lulog,11) time(1 ),timesec(1 ), &
       time(nt),timesec(nt), &
       INT(timesec(2)-timesec(1))
11 format('DBS time coverage : ',/, &
       ' From         : ',f15.0,' ( ',f8.0,' sec after 0000UTC)',/, &
       ' To           : ',f15.0,' ( ',f8.0,' sec after 0000UTC)',/, &
       ' Increment (s): ',i6)    !

  !
  !***  Calculates time lag by iteration (that is, the time in seconds between the PROF origin
  !***  and the DBS origin). Both origins are referred to 0000UTC, but may belong to different
  !***  days
  !
  found = .false.
  time_lag = 0.0_rp
  do while(.not.found)
     call addtime(ibyr_PROF,ibmo_PROF,ibdy_PROF,0, &
          iyr,imo,idy,ihr,imi,ise,time_lag)
     if( (iyr.eq.ibyr).and.(imo.eq.ibmo).and.(idy.eq.ibdy) ) then
        found = .true.
     else
        time_lag = time_lag + 86400.
     end if
     if(time_lag.gt.timesec_PROF(nt_PROF)) call runend('Time lag not found. Check consistency in dates')
  end do
  write(lulog,20) INT(time_lag/86400.0_rp),time_lag
20 format(' Time lag     : ',i3,' days (',f7.0,' sec)')
  !
  !*** Checks that the required time is within the bounds
  !
  if((time_lag+timesec(1 )).lt.timesec_PROF(1     )) &
       call runend('checktime_PROF : Minimum profile time is after required time')
  if((time_lag+timesec(nt)).gt.timesec_PROF(nt_PROF)) &
       call runend('checktime_PROF : Maximum profile time is before required time')
  !
  return
  !
  !*** List of errors
  !
100 call runend('checktime_PROF: error opening profile file '//TRIM(ludatname))
101 call runend('checktime_PROF: error reading profile file'//TRIM(ludatname))
  !
  return
end subroutine checktime_PROF
