subroutine checktime_ECMWF
  !*************************************************************************
  !*
  !*    Outputs the ECMWF and dbs time coverages and gets the time lag
  !*    (if any)
  !*
  !*************************************************************************
  use kindType
  use timeFun
  use Master
  use ECMWF_nc
  use InpOut
  implicit none
  !
  logical     :: found
  integer(ip) :: iyr ,imo ,idy ,ihr ,imi ,ise
  !
  !***  Writes
  !
  write(lulog,10) time_ECMWF(1     ),timesec_ECMWF(1)     , &
       time_ECMWF(nt_ECMWF),timesec_ECMWF(nt_ECMWF), &
       INT(timesec_ECMWF(2)-timesec_ECMWF(1))
10 format('ECMWF time coverage : ',/, &
       ' From         : ',f15.0,' ( ',f8.0,' sec after 0000UTC)',/, &
       ' To           : ',f15.0,' ( ',f8.0,' sec after 0000UTC)',/, &
       ' Increment (s): ',i6)
  write(lulog,11) time(1 ),timesec(1 ), &
       time(nt),timesec(nt), &
       INT(timesec(2)-timesec(1))
11 format('DBS time coverage : ',/, &
       ' From         : ',f15.0,' ( ',f8.0,' sec after 0000UTC)',/, &
       ' To           : ',f15.0,' ( ',f8.0,' sec after 0000UTC)',/, &
       ' Increment (s): ',i6)    !
  !
  !***  Calculates time lag by iteration (that is, the time in seconds between the ECMWF origin
  !***  and the DBS origin). Both origins are referred to 0000UTC, but may belong to different
  !***  days
  !
  found = .false.
  time_lag = 0.0_rp
  do while(.not.found)
     call addtime(ibyr_ECMWF,ibmo_ECMWF,ibdy_ECMWF,0, &
          iyr,imo,idy,ihr,imi,ise,time_lag)
     !
     if( (iyr.eq.ibyr).and.(imo.eq.ibmo).and.(idy.eq.ibdy) ) then
        found = .true.
     else
        time_lag = time_lag + 86400.0_rp
     end if
     if(time_lag.gt.timesec_ECMWF(nt_ECMWF)) call runend('Time lag not found. Check interval consistency')
  end do
  write(lulog,20) INT(time_lag/86400.0_rp),time_lag
20 format(' Time lag     : ',i3,' days (',f7.0,' sec)')
  !
  !*** Checks that the required time is within the bounds
  !
  if((time_lag+timesec(1 )).lt.timesec_ECMWF(1     )) &
       call runend('read_ECMWF : Minimum datatime is after  required time')
  if((time_lag+timesec(nt)).gt.timesec_ECMWF(nt_ECMWF)) &
       call runend('read_ECMWF : Maximum datatime is before required time')
  !
  return
end subroutine checktime_ECMWF
