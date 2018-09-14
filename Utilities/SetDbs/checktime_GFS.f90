subroutine checktime_GFS
  !*************************************************************************
  !*
  !*    Outputs the GFS and dbs time coverages and gets the time lag
  !*    (if any)
  !*
  !*************************************************************************
  use kindType
  use timeFun
  use Master
  use GFS_nc
  use InpOut
  implicit none
  !
  logical     :: found
  integer(ip) :: iyr ,imo ,idy ,ihr ,imi ,ise
  !
  !***  Writes
  !
  write(lulog,10) time_GFS(1     ),timesec_GFS(1)     , &
       time_GFS(nt_GFS),timesec_GFS(nt_GFS), &
       INT(timesec_GFS(2)-timesec_GFS(1))
10 format('GFS time coverage : ',/, &
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
  !***  Calculates time lag by iteration (that is, the time in seconds between the GFS origin
  !***  and the DBS origin). Both origins are referred to 0000UTC, but may belong to different
  !***  days
  !
  found = .false.
  time_lag = 0.0_rp
  do while(.not.found)
     call addtime(ibyr_GFS,ibmo_GFS,ibdy_GFS,0, &
          iyr,imo,idy,ihr,imi,ise,time_lag)
     !
     if( (iyr.eq.ibyr).and.(imo.eq.ibmo).and.(idy.eq.ibdy) ) then
        found = .true.
     else
        time_lag = time_lag + 86400.0_rp
     end if
     if(time_lag.gt.timesec_GFS(nt_GFS)) call runend('Time lag not found. Check interval consistency')
  end do
  write(lulog,20) INT(time_lag/86400.0_rp),time_lag
20 format(' Time lag     : ',i3,' days (',f7.0,' sec)')
  !
  !*** Checks that the required time is within the bounds
  !
  if((time_lag+timesec(1 )).lt.timesec_GFS(1     )) &
       call runend('read_GFS : Minimum datatime is after  required time')
  if((time_lag+timesec(nt)).gt.timesec_GFS(nt_GFS)) &
       call runend('read_GFS : Maximum datatime is before required time')
  !
  return
end subroutine checktime_GFS
